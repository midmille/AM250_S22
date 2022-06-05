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
    &                         Ns)
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
        integer, intent(out)                     :: Ns

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
        Ns = Nx / Np

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
        integer, intent(inout)                   :: Asg(0:Nysp1, 0:Nxsp1)
        ! [The number of cols, rows of the parent matrix.]
        integer, intent(in)                      :: Nx, Ny

        ! Returns
        ! -------
        ! [The parent matrix, only allocated and initialized on master.]
        integer, intent(in)                      :: A(Ny, Nx)

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
        As_flt = 
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




    END SUBROUTINE Gather_Life


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
        ! [The number of columns per processor.]
        integer                              :: Ns
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
        CALL Partition_Cols(Nx, Np, Ns)
        ! [Setting the size of each sub domain given Ns from above.]
        Nxs = Ns
        Nys = Ny
        Nxsp1 = Ns + 1
        Nysp1 = Ny + 1

        ! [Allocate the partitions of domain to the partition size.]
        CALL Init_Partition(Nxsp1, Nysp1, Asg)

        PRINT *, "Asg before scatter:"
        DO i=0,Nysp1
            DO j=0,Nxsp1
                WRITE(*,'(I2)', advance='no') Asg(i,j)
            ENDDO
            WRITE(*,*) ''
        ENDDO


        ! [Broadcat designated partition to each processor from A on master.]
        CALL Share_Life(pid,                                             &
        &               master,                                          &
        &               Nx,                                              &
        &               Ny,                                              &
        &               A,                                               &
        &               Nxs,                                             &
        &               Nys,                                             & 
        &               Nxsp1,                                           &
        &               Nysp1,                                           &
        &               Asg)

        
        PRINT *, "Asg after scatter:", pid
        DO i=0,Nysp1
            DO j=0,Nxsp1
                WRITE(*,'(I2)', advance='no') Asg(i,j)
            ENDDO
            WRITE(*,*) ''
        ENDDO
     
        ! [Loop time.]
!        DO k=1,Nt 
!
!            ! [If the write output flag is true.]
!            IF (woflag .EQ. .TRUE.) THEN 
!               ! [Write first step and the output every Nw time step.]
!                IF ((k .EQ. 1) .OR. (MOD(k,Nw) .EQ. 0)) THEN 
!                    CALL Gather_Write_Life_Step()
!
!                ENDIF 
!            ENDIF
! 
!            ! [Updated the ghost nodes on each partition. This includes the
!            !  necessary ringed send and receives of the columns.]
!            CALL UGN()
!
!            ! [Run life on each partition.]
!            CALL Life()
!        
!        ENDDO 


    END SUBROUTINE Run_Column_Life

 
END MODULE cols_mod
