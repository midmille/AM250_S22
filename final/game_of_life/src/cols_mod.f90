! AM250_S22/final/game_of_life/src/cols_mod.f90
! Author : Miles D. Miller, University of California Santa Cruz
! Created : 9:30 AM, June 4, 2022
! About : This is the module for the parallelization of the game of life
!         using a column based partitioning.


MODULE cols_mod
    
    USE mpi 

    USE funcs_mod, ONLY: Life, Write_Life_Step

    implicit none

    ! Parameters
    ! ----------
    ! [Double.]
    integer, parameter                 :: dp = kind(0.d0)


CONTAINS
    

    !----------------------------------------------------------------
    SUBROUTINE Partition_Cols(Nx,                                   & 
    &                         Np,                                   &
    &                         Nxs)
    !----------------------------------------------------------------
        ! This routine is for the partitioning of the domain to different 
        ! processors. 
        ! It checks for the cohesion of the desired number of processors
        ! to the size of the domain. 
        ! Since the partition is for columns we only have to work with the
        ! number of columns.

        ! Parameters
        ! ----------
        ! [The number of columns.] 
        integer, intent(in)                      :: Nx
        ! [The number of processors.]
        integer, intent(in)                      :: Np

        ! Returns
        ! -------
        ! [The number of columns per processor.]  
        integer, intent(out)                     :: Nxs

        ! Variables
        ! ---------
        ! None 

        ! ROUTINE START
        ! =============
        
        ! [Require that processor count is less than or equal to number of columns.]
        IF (Np .GT. Nx) THEN 
            PRINT *, "The number of columns is less than number of processors"
            STOP "Error: Bad column number and processor count cohesion."
        ENDIF

        ! [Require that the number of columns be evenly divisible by the number of processor.]
        IF (MOD(Nx, Np) .NE. 0) THEN 
            PRINT *, "For column partition, please ensure MOD(Nx,Np) == 0"
            STOP "Error: Bad column number and processor count cohesion."
        ENDIF

        ! [This gives the number of columns per parttition.] 
        Nxs = Nx / Np

    END SUBROUTINE 


    !-----------------------------------------------------------------
    SUBROUTINE Init_Partition(Nxsp1,                                 &
    &                         Nysp1,                                 &
    &                         Asg)
    !-----------------------------------------------------------------
        ! This allocates the memory to the partition array.
        
        ! Parameters
        ! ----------
        ! [The number of columns in the partitioned matrix sub domain.]
        integer, intent(in)                      :: Nxsp1
        ! [The number of rows in the partitioned matrix sub domain.]
        integer, intent(in)                      :: Nysp1

        ! Returns
        ! -------
        ! [The array with to allocate with ghost nodes.]
        integer, allocatable, intent(inout)      :: Asg(:,:)

        ! ROUTINE START
        ! =============

        ! [Allocate the array with ghost cells on every processor.]
        ALLOCATE(Asg(0:Nysp1, 0:Nxsp1))
        Asg = 0.0

    END SUBROUTINE Init_Partition

    
    !-----------------------------------------------------------------
    SUBROUTINE Share_Life(pid,                                       &
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

        ! [Make A flat such that scatter is contiguous]
        A_flt = RESHAPE(A, (/Nx*Ny/))
        PRINT *, 'Shape A_flt', SHAPE(A_flt)
        PRINT *, 'Shape As_flt', SHAPE(As_flt)

        PRINT *, 'Send Count:', send_cnt
        PRINT *, 'Recv Count:', recv_cnt

        ! [Scatter]
        CALL MPI_SCATTER(A_flt,                                      &
        &                send_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                As_flt,                                     &
        &                recv_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                master,                                     &
        &                MPI_COMM_WORLD,                             &
        &                ierr)

        Asg(1:Nys,1:Nxs) = RESHAPE(As_flt, (/Nys, Nxs/))

    
    END SUBROUTINE Share_Life


    !-----------------------------------------------------------------
    SUBROUTINE Gather_Life(pid,                                      &
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

        ! [Make Ags flat such that gather is contiguous]
        As_flt = RESHAPE(Asg(1:Nys, 1:Nxs), (/Nys*Nxs/))

        ! [Scatter]
        CALL MPI_GATHER(As_flt,                                      &
        &               send_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               A_flt,                                       &
        &               recv_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               master,                                      &
        &               MPI_COMM_WORLD,                              &
        &               ierr)

        A = RESHAPE(A_flt, (/Ny, Nx/))

    END SUBROUTINE Gather_Life

    
    !-----------------------------------------------------------------
    SUBROUTINE Share_BC(pid,                                         &
    &                   right,                                       &
    &                   left,                                        &
    &                   cnt,                                         &
    &                   sendR_bc,                                    &
    &                   sendL_bc,                                    &
    &                   recvL_bc,                                    &
    &                   recvR_bc)
    !-----------------------------------------------------------------
        ! This uses MPI ISEND and IRECV to share the column BCs of the 
        ! sub matrices with their neighbours.
        ! It assumes the domain is deconstructed such that neighbouring
        ! processors have neighbouring partitions of the parent domain.
        ! For example somthing like:
        !     p1           p2          p3      ....
        !  a1 a2 a3     a4 a5 a6    a7 a8 a9   .... 
        ! This should be accouplished if the partitionin of the columns 
        ! is done using the scatter function.

        ! Parameters
        ! ----------
        ! [pid.]
        integer, intent(in)                      :: pid
        ! [The processor number of the left and right.]
        integer, intent(in)                      :: right, left
        ! [The count of the buffer to send and receive.]
        integer, intent(in)                      :: cnt
        ! [The BC vector buffer to send to the right.]
        integer, intent(in)                      :: sendR_bc(cnt)
        ! [The BC vector buffer to send to the left.]
        integer, intent(in)                      :: sendL_bc(cnt)

        ! Returns
        ! -------
        ! [The BC vector buffer to receive from the left.]
        integer, intent(inout)                   :: recvL_bc(cnt)
        ! [The BC vector buffer to receive from the right.]
        integer, intent(inout)                   :: recvR_bc(cnt)

        ! Variables
        ! ---------
        ! [The tags for the right send and the left send.]
        integer                                  :: tagR, tagL
        ! [The requests for the wait all, send/recv and right/left.]
        integer                                  :: request_sendR, request_recvL
        integer                                  :: request_sendL, request_recvR
        ! [The ierr returned from each send/recv/]
        integer                                  ::ierr
        ! [The array of requests.]
        integer                                  :: requests(4) 
        ! [The returned array of statuses the 4 is for the 4 requests.]
        integer                                  :: statuses(MPI_STATUS_SIZE, 4)


        ! ROUTINE START
        ! =============

        ! [Sending the data to the left and right in the ring.]
        ! [The tag for the right send, left recv pair.]
        tagR = 12
        ! [The tag for the left send, right recv pair.]
        tagL = 34

        ! [Send the data to the right.] 
        CALL MPI_ISEND(sendR_bc,                               &
        &              cnt,                                    &
        &              MPI_INTEGER,                            &
        &              right,                                  &
        &              tagR,                                   &
        &              MPI_COMM_WORLD,                         &
        !              [OUTPUT]                                ! 
        &              request_sendR,                          &
        &              ierr)
        ! [Receive data from the left.]
        CALL MPI_IRECV(recvL_bc,                               &
        &              cnt,                                    &
        &              MPI_INTEGER,                            &
        &              left,                                   &
        &              tagR,                                   &
        &              MPI_COMM_WORLD,                         &
        !              [OUTPUT]                                ! 
        &              request_recvL,                          &
        &              ierr)

        ! [Send the data to the left.] 
        CALL MPI_ISEND(sendL_bc,                               &
        &              cnt,                                    &
        &              MPI_INTEGER,                            &
        &              left,                                   &
        &              tagL,                                   &
        &              MPI_COMM_WORLD,                         &
        !              [OUTPUT]                                ! 
        &              request_sendL,                          &
        &              ierr)
        ! [Receive data from the right.]
        CALL MPI_IRECV(recvR_bc,                               &
        &              cnt,                                    &
        &              MPI_INTEGER,                            &
        &              right,                                  &
        &              tagL,                                   &
        &              MPI_COMM_WORLD,                         &
        !              [OUTPUT]                                ! 
        &              request_recvR,                          &
        &              ierr)


        ! [The request array.]
        requests = (/request_sendR, request_recvL, request_sendL, request_recvR/)
    
        ! [Wait for the recv and send to compelete.]
        CALL MPI_WAITALL(4, requests, statuses, ierr)
     
    END SUBROUTINE Share_BC


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
        CALL Share_Life(pid,                                         &
        &               master,                                      &
        &               Nx,                                          &
        &               Ny,                                          &
        &               A,                                           &
        &               Nxs,                                         &
        &               Nys,                                         & 
        &               Nxsp1,                                       &
        &               Nysp1,                                       &
        &               Asg)

        
        DO k=0,Np
            IF (pid .EQ. k) THEN      
                PRINT *, "pid:", pid
                PRINT *, "A after scatter and gather:"
                DO i=1,Ny
                    DO j=1,Nx
                        WRITE(*,'(I2)', advance='no') A(i,j)
                    ENDDO
                    WRITE(*,*) ''
                ENDDO
            ENDIF
        ENDDO
     
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
                    CALL Gather_Life(pid,                                        &
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
                        PRINT *, "pid:", pid
                        PRINT *, "A after scatter and gather:"
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
