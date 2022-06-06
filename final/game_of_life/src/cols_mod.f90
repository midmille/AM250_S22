! AM250_S22/final/game_of_life/src/cols_mod.f90
! Author : Miles D. Miller, University of California Santa Cruz
! Created : 9:30 AM, June 4, 2022
! About : This is the module for the parallelization of the game of life
!         using a column based partitioning.


MODULE cols_mod
    
    USE mpi 

    USE funcs_mod, ONLY: Life, Write_Life_Step, Init_Partition, Partition_Cols 
    USE funcs_mod, ONLY: Share_BC

    implicit none

    ! Parameters
    ! ----------
    ! [Double.]
    integer, parameter                 :: dp = kind(0.d0)


CONTAINS


    !-----------------------------------------------------------------
    SUBROUTINE Share_Life_Cols(pid,                                  &
    &                     master,                                    &
    &                     Nx,                                        &
    &                     Ny,                                        &
    &                     A,                                         &
    &                     Nxs,                                       &
    &                     Nys,                                       & 
    &                     Nxsp1,                                     &
    &                     Nysp1,                                     &
    &                     Asg)
    !-----------------------------------------------------------------
        ! This algorithm implements the scatter to share the parent matrix
        ! from the master processor to the child sub matrices of the 
        ! rest of the processors.  

        ! Parameters
        ! ----------
        ! [Pid.]
        integer, intent(in)                      :: pid
        ! [The master procs id.]
        integer, intent(in)                      :: master
        ! [The number of cols, rows of the parent matrix.]
        integer, intent(in)                      :: Nx, Ny
        ! [The parent matrix, only allocated and initialized on master.]
        integer, intent(in)                      :: A(Ny, Nx)
        ! [The number of cols, rows of child matrix partition.]
        integer, intent(in)                      :: Nxs, Nys
        ! [Their p1 counterparts.]
        integer, intent(in)                      :: Nxsp1, Nysp1

        ! Returns
        ! -------
        ! [The child matrix with ghost nodes.]
        integer, intent(inout)                   :: Asg(0:Nysp1, 0:Nxsp1)

        ! Variables
        ! ----------
        ! [The send and recv count.]
        integer                                  :: send_cnt, recv_cnt 
        ! [A flattened]
        integer                                  :: A_flt(Nx*Ny)
        ! [Asg flattened.]
        integer                                  :: As_flt(Nxs*Nys)
        ! [The error returned error code.]
        integer                                  :: ierr


        ! ROUTINE START
        ! =============

        ! [The send and recv counts.]
        send_cnt = Nys*Nxs
        recv_cnt = Nys*Nxs

        ! [Scatter sub matrix partitions as contiguously found in memory.]
        CALL MPI_SCATTER(A,                                          &
        &                send_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                Asg(1:Nys,1:Nxs),                           &
        &                recv_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                master,                                     &
        &                MPI_COMM_WORLD,                             &
        &                ierr)
    
    END SUBROUTINE Share_Life_Cols


    !-----------------------------------------------------------------
    SUBROUTINE Gather_Life_Cols(pid,                                 &
    &                      master,                                   &
    &                      Nxs,                                      & 
    &                      Nys,                                      &
    &                      Nxsp1,                                    &
    &                      Nysp1,                                    &
    &                      Asg,                                      &
    &                      Nx,                                       &
    &                      Ny,                                       &
    &                      A)
    !-----------------------------------------------------------------
        ! This algorithm implements the mpi gether routine to gather the 
        ! sub matrices back into the parent matrix.

        ! Parameters
        ! ----------
        ! [Pid.]
        integer, intent(in)                      :: pid
        ! [The master procs id.]
        integer, intent(in)                      :: master
        ! [The number of cols, rows of child matrix partition.]
        integer, intent(in)                      :: Nxs, Nys
        ! [Their p1 counterparts.]
        integer, intent(in)                      :: Nxsp1, Nysp1
        ! [The child matrix with ghost nodes.]
        integer, intent(in)                      :: Asg(0:Nysp1, 0:Nxsp1)
        ! [The number of cols, rows of the parent matrix.]
        integer, intent(in)                      :: Nx, Ny

        ! Returns
        ! -------
        ! [The parent matrix, only allocated and initialized on master.]
        integer, intent(inout)                   :: A(Ny, Nx)

        ! Variables
        ! ----------
        ! [The send and recv count.]
        integer                                  :: send_cnt, recv_cnt 
        ! [A flattened]
        integer                                  :: A_flt(Nx*Ny)
        ! [Asg flattened.]
        integer                                  :: As_flt(Nxs*Nys)
        ! [The error returned error code.]
        integer                                  :: ierr

        ! ROUTINE START
        ! =============

        ! [The send and recv counts.]
        send_cnt = Nys*Nxs
        recv_cnt = Nys*Nxs

        ! [Gather partition sub arrays contigously in memory.]
        CALL MPI_GATHER(Asg(1:Nys,1:Nxs),                            &
        &               send_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               A,                                           &
        &               recv_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               master,                                      &
        &               MPI_COMM_WORLD,                              &
        &               ierr)
        
    END SUBROUTINE Gather_Life_Cols


    !-----------------------------------------------------------------
    SUBROUTINE UGN_Cols(pid,                                         & 
    &                   Np,                                          &
    &                   Nxs,                                         &
    &                   Nys,                                         &
    &                   Nxsp1,                                       &
    &                   Nysp1,                                       &
    &                   Asg)
    !-----------------------------------------------------------------
        ! This routine updates the ghots nodes of the sub partitioned 
        ! matrices. It implements a column send and receive of the BCs
        ! with the neighibour to the left and right of the ring at the 
        ! same time. It then updates the ghost nodes of each sub matrix 
        ! with the BC from neighbours.

        ! Parameters
        ! ----------
        ! [The pid.]
        integer, intent(in)                      :: pid
        ! [numprocs.]
        integer, intent(in)                      :: Np
        ! [Nxs, Nys]
        integer, intent(in)                      :: Nxs, Nys
        ! [Their p1 countrerparts.]
        integer, intent(in)                      :: Nxsp1, Nysp1


        ! Returns
        ! -------
        ! [The sub matrix on each procs, whose ghost nodes are updated.]
        integer, intent(inout)                   :: Asg(0:Nysp1, 0:Nxsp1)

        ! Variables
        ! ---------
        ! [The right and left pids.] 
        integer                                  :: right, left    
        ! [The count of send vector.]
        integer                                  :: cnt

        ! ROUTINE START
        ! =============

        ! [Create the left and right neighbour ids.]
        ! [For the pid on the left most side (total domain has periodic BC).]
        IF (pid .EQ. 0) THEN
            ! [left is the id of the last processor in the group.]
            left = Np - 1
            right = pid + 1
        ! [For the processor on the right most side.]
        ELSE IF (pid .EQ. Np-1) THEN 
            left = pid - 1
            ! [The right processor is just the first processor.]
            right = 0
        ELSE 
            left = pid - 1
            right = pid + 1 
        END IF

        ! [First implement row exchange.]
        Asg(0,1:Nxs) = Asg(Nys, 1:Nxs)
        Asg(Nysp1,1:Nxs) = Asg(1,1:Nxs)

        ! [The vecotr being sent is an entire row including ghots cells.]
        cnt = Nysp1+1

        ! [Send and receive the columns for the
        CALL Share_BC(pid,                                           &
        &             right,                                         &
        &             left,                                          &
        &             cnt,                                           &
!        &             sendR_bc,                                      &
        &             Asg(0:Nysp1,Nxs),                              &
!        &             sendL_bc,                                      &
        &             Asg(0:Nysp1,1),                                &
!        &             recvL_bc,                                      &
        &             Asg(0:Nysp1, 0),                               &
!        &             recvR_bc)
        &             Asg(0:Nysp1, Nxsp1))                

    END SUBROUTINE UGN_Cols


    !-----------------------------------------------------------------
    SUBROUTINE Run_Column_Life(pid,                                  & 
    &                          master,                               &
    &                          Np,                                   &
    &                          Nx,                                   &
    &                          Ny,                                   &
    &                          Nt,                                   & 
    &                          Nw,                                   &
    &                          A,                                    &
    &                          savefile_head,                        &
    &                          cols_outdir,                          &
    &                          woflag,                               &
    !                          [OUTPUT]                              ! 
    &                          t)
    !-----------------------------------------------------------------
        ! This routine solves the game of life by partitioning the
        ! parallelization over the columns of the domain.

        ! Parameters
        ! ----------
        ! [The processor id.]
        integer, intent(in)                  :: pid
        ! [The master id.]
        integer, intent(in)                  :: master
        ! [The total number of processors.]
        integer, intent(in)                  :: Np
        ! [The x,y,t dimensions.]
        integer, intent(in)                  :: Nx, Ny, Nt
        ! [The write interval.]
        integer, intent(in)                  :: Nw
        ! [The parent life matrix.]
        integer, intent(in)                  :: A(Ny, Nx)
        ! [The save file header.]
        character(len=*), intent(in)         :: savefile_head
        ! [The serial output directory.]
        character(len=*), intent(in)         :: cols_outdir
        ! [The write output flag.]
        logical, intent(in)                  :: woflag

        ! Returns
        ! -------
        ! [The elapsed time for the algorithm.]
        real (dp), intent(out)               :: t

        ! Variables
        ! ---------
        ! [The a array at a given step.]
        integer                              :: A_k(Ny, Nx)
        ! [The size of the ghosted sub domain matrix for each.]
        integer                              :: Nxs, Nys
        integer                              :: Nxsp1, Nysp1
        ! [The paritioned sub domain for each processor. With ghots cells.]
        integer, allocatable                 :: Asg(:,:) 
        ! [Time indexing.]
        integer                              :: k
        ! [Temp other indexing.]
        integer::i,j

        ! ROUTINE START
        ! =============

        ! [Temp place holder for time output.]
        t = 0
        ! [First step is to partition the problem according to the domain size, 
        !  the number of processors, and the designmated decomposition.]
        CALL Partition_Cols(Nx, Np, Nxs)
        ! [Setting the size of each sub domain given Ns from above.]
        Nys = Ny
        Nxsp1 = Nxs + 1
        Nysp1 = Nys + 1

        ! [Allocate the partitions of domain to the partition size.]
        CALL Init_Partition(Nxsp1, Nysp1, Asg)

        ! [Broadcat designated partition to each processor from A on master.]
        CALL Share_Life_Cols(pid,                                    &
        &               master,                                      &
        &               Nx,                                          &
        &               Ny,                                          &
        &               A,                                           &
        &               Nxs,                                         &
        &               Nys,                                         & 
        &               Nxsp1,                                       &
        &               Nysp1,                                       &
        &               Asg)
        
        ! [Loop time.]
        DO k=1,Nt 

            ! [Updated the ghost nodes on each partition. This includes the
            !  necessary ringed send and receives of the columns.]
            CALL UGN_Cols(pid,                                       & 
            &             Np,                                        &
            &             Nxs,                                       &
            &             Nys,                                       &
            &             Nxsp1,                                     &
            &             Nysp1,                                     &
            &             Asg)
 
            ! [Run life on each partition.]
            CALL Life(Nxsp1, Nysp1, Asg)

            ! [If the write output flag is true.]
            IF (woflag .EQ. .TRUE.) THEN 
               ! [Write first step and the output every Nw time step.]
                IF ((k .EQ. 1) .OR. (MOD(k,Nw) .EQ. 0)) THEN 
                    ! [Gather everything to A_k, which is A at step k.]
                    CALL Gather_Life_Cols(pid,                                   &
                    &                master,                                     &
                    &                Nxs,                                        & 
                    &                Nys,                                        &
                    &                Nxsp1,                                      &
                    &                Nysp1,                                      &
                    &                Asg,                                        &
                    &                Nx,                                         &
                    &                Ny,                                         &
                    &                A_k)

                    IF (pid .EQ. master) THEN      
                        PRINT *, "A after scatter and gather:"
                        PRINT *, "---------------------------"
                        DO i=1,Ny
                            DO j=1,Nx
                                WRITE(*,'(I2)', advance='no') A_k(i,j)
                            ENDDO
                            WRITE(*,*) ''
                        ENDDO
                    ENDIF

                    ! [Write the result on master processor.]
                    IF (pid .EQ. master) THEN
                        ! [Write A_k.]
                        CALL Write_Life_Step(k,                      &
                        &                    Nx,                     &
                        &                    Ny,                     &
                        &                    A_k,                    &
                        &                    savefile_head,          &
                        &                    cols_outdir)
                    ENDIF
                ENDIF 
            ENDIF
        ENDDO 


    END SUBROUTINE Run_Column_Life

 
END MODULE cols_mod
