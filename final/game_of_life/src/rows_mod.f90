! Created : 9:30 AM, June 4, 2022
! About : This is the module for the parallelization of the game of life
!         using a rows based partitioning.
MODULE rows_mod
    
    USE mpi 

    USE funcs_mod, ONLY: Life, Write_Life_Step, Init_Partition, Partition_Rows
    USE funcs_mod, ONLY: Share_BC

    implicit none

    ! Parameters
    ! ----------
    ! [Double.]
    integer, parameter                 :: dp = kind(0.d0)


CONTAINS

    
    !----------------------------------------------------------------
    SUBROUTINE Create_Rowtype(Nx,                                    &
    &                         Ny,                                    &
    &                         Nxs,                                   &
    &                         Nys,                                   &
    &                         rowtype)
    !----------------------------------------------------------------
        ! This routine creates and commits the rows type data type for
        ! the provided matrix dimensions.

        ! Parameters
        ! ----------
        ! [The number of columns and rows for of the matrix this rowtype i
        !  will belong to.]
        integer, intent(in) :: Nx, Ny, Nxs, Nys

        ! Returns
        ! -------
        ! [The rowtype for the designated matrix dims.]
        integer, intent(out)                     :: rowtype

        ! Variables
        ! ---------
        ! [The temporary rowtype with large extent.]
        integer                                  :: rowtypetmp
        ! [The err] 
        integer                                  :: ierr
        ! [The extent and lower bound, necessarily kind=MPI_ADDRESS_KIND.] 
        integer(kind=MPI_ADDRESS_KIND)           :: lb, extent
        ! [The type size for MPI_INTEGER.]
        integer                                  :: typesize

        ! ROUTINE START
        ! =============
        
        ! [Creating the rowtype.]
        ! [The Count is the number of columns, so the entries in each row.
        !  The block is the number of rows in the sub matrix, Nys, we are grouping 
        !  the rows by partition.
        !  The stride is the number of columns, the stride includes the block.
        !  See professor Brummels MPI notes for the image of count, block and
        !  stride. 
        !  Also see the SCAI MPI Derived Datatypes in main.f90 references.
        CALL MPI_TYPE_VECTOR(Nx, Nys, Ny, MPI_INTEGER, rowtypetmp, ierr)
!        CALL MPI_TYPE_COMMIT(rowtypetmp, ierr)

        ! [Retreive the bit size of an MPI_INTEGER.]
        CALL MPI_TYPE_SIZE(MPI_INTEGER, typesize, ierr)

!        PRINT *, "MPI type size of rowtypetmp", typesize
!        CALL MPI_TYPE_GET_EXTENT(rowtypetmp, lb, extent, ierr)
!        PRINT *, "MPI type extent", lb, extent
!        CALL MPI_TYPE_SIZE(MPI_INTEGER, typesize, ierr)
!        PRINT *, "MPI type size of MPI_INTEGER", typesize

        ! [Set the lower bound for the extended data type to be one since A is 
        !  one based index.]
        lb = 1

        ! [Set the extend to be the 
        extent = Nys*typesize
        ! [Form the resized the rowtype. It will have size Nx*Nys*typesize and
        ! extent Nys*typesize. So it will have a the column size of a row 
        ! submatrix partition as the extent and the number of those column
        ! chunks as the size. I am not exatly sure why this works however....
        ! Must ask Brummel during office hours.] 
        CALL MPI_TYPE_CREATE_RESIZED(rowtypetmp, lb, extent, rowtype, ierr)

        ! [Commiting the rowtype.]
        CALL MPI_TYPE_COMMIT(rowtype, ierr)


!        CALL MPI_TYPE_GET_EXTENT(rowtype, lb, extent, ierr)
!        PRINT *, "MPI typ  extent rotype", lb, extent
!        CALL MPI_TYPE_SIZE(rowtype, typesize, ierr)
!        PRINT *, "MPI type size of rowtype", typesize


    END SUBROUTINE Create_Rowtype

    !-----------------------------------------------------------------
    SUBROUTINE Share_Life_Rows(pid,                                  &
    &                     master,                                    &
    &                     rowtype,                                   &
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
        ! [The MPI rowtype from MPI_TYPE_VECTOR.]
        integer, intent(in)                      :: rowtype
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
        integer                                  :: extent


        ! ROUTINE START
        ! =============

        ! [The send and recv counts.]
        send_cnt = 1
        recv_cnt = Nys*Nxs
!        PRINT *, 'Send count:', send_cnt
!        PRINT *, 'Recv count:', recv_cnt
!        PRINT *, "Ny", Ny
!        PRINT *, "Nys", Nys

!        A_flt = RESHAPE(A, (/Nx*Ny/))

        ! [Scatter sub matrix partitions according to row type.]
        CALL MPI_SCATTER(A(:,:),                                     &
        &                send_cnt,                                   &
        &                rowtype,                                    &
        &                Asg(1:Nys,1:Nxs),                           &
        &                recv_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                master,                                     &
        &                MPI_COMM_WORLD,                             &
        &                ierr)
 
!        ! [Make A flat such that scatter is contiguous]
!        A_flt = RESHAPE(A, (/Nx*Ny/))

        !Asg(1:Nys,1:Nxs) = RESHAPE(As_flt, (/Nys, Nxs/))

    
    END SUBROUTINE Share_Life_Rows


    !-----------------------------------------------------------------
    SUBROUTINE Gather_Life_Rows(pid,                                 &
    &                      master,                                   &
    &                      rowtype,                                  &
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
        ! [The row type for the handling of A in the scatter gather.]
        integer, intent(in)                      :: rowtype
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
        recv_cnt = 1

        ! [Make Ags flat such that gather is contiguous]
!        As_flt = RESHAPE(Asg(1:Nys, 1:Nxs), (/Nys*Nxs/))

        ! [Gather partition sub arrays contigously in memory.]
        CALL MPI_GATHER(Asg(1:Nys,1:Nxs),                            &
        &               send_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               A(:,:),                                      &
        &               recv_cnt,                                    &
        &               rowtype,                                     &
        &               master,                                      &
        &               MPI_COMM_WORLD,                              &
        &               ierr)
       
!        A = RESHAPE(A_flt, (/Ny, Nx/))

    END SUBROUTINE Gather_Life_Rows

   

    !-----------------------------------------------------------------
    SUBROUTINE UGN_Rows(pid,                                         & 
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

        ! [First implement col exchange.]
        Asg(1:Nys,0) = Asg(1:Nys, Nxs)
        Asg(1:Nys,Nxsp1) = Asg(1:Nys,1)

        ! [The vecotr being sent is an entire row including ghots cells.]
        cnt = Nxsp1+1

        ! [Send and receive the columns for the
        CALL Share_BC(pid,                                           &
        &             right,                                         &
        &             left,                                          &
        &             cnt,                                           &
!        &             sendR_bc,                                      &
        &             Asg(Nys,0:Nxsp1),                              &
!        &             sendL_bc,                                      &
        &             Asg(1,0:Nxsp1),                                &
!        &             recvL_bc,                                      &
        &             Asg(0,0:Nxsp1),                               &
!        &             recvR_bc)
        &             Asg(Nysp1, 0:Nxsp1))                

    END SUBROUTINE UGN_Rows


    !-----------------------------------------------------------------
    SUBROUTINE Run_Row_Life(pid,                                     & 
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
        ! [The rowtype from MPI_TYPE_VECTOR.]
        integer                              :: rowtype, rowtypetmp
        ! [Time indexing.]
        integer                              :: k
        ! [Temp other indexing.]
        integer::i,j
        integer::ierr

        ! ROUTINE START
        ! =============

        ! [Temp place holder for time output.]
        t = 0
        ! [First step is to partition the problem according to the domain size, 
        !  the number of processors, and the designmated decomposition.]
        CALL Partition_Rows(Nx, Np, Nys)
        ! [Setting the size of each sub domain given Ns from above.]
        Nxs = Nx
        Nxsp1 = Nxs + 1
        Nysp1 = Nys + 1

        ! [Allocate the partitions of domain to the partition size.]
        CALL Init_Partition(Nxsp1, Nysp1, Asg)
        PRINT *, 'SHAPE(Asg):', SHAPE(Asg)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        ! [Create the rowtype for the matrix A for its gather and scatter
        !  operations.]
        CALL Create_Rowtype(Nx, Ny, Nxs, Nys, rowtype)

        ! [Get the rowtype for the scatter and gather of A.]

        IF (pid .EQ. master) THEN      
            PRINT *, "A before scatter:"
            PRINT *, "---------------------------"
            DO i=1,Ny
                DO j=1,Nx
                    WRITE(*,'(I2)', advance='no') A(i,j)
                ENDDO
                WRITE(*,*) ''
            ENDDO
        ENDIF

        ! [Broadcat designated partition to each processor from A on master.]
        CALL Share_Life_Rows(pid,                                    &
        &               master,                                      &
        &               rowtype,                                     &
        &               Nx,                                          &
        &               Ny,                                          &
        &               A,                                           &
        &               Nxs,                                         &
        &               Nys,                                         & 
        &               Nxsp1,                                       &
        &               Nysp1,                                       &
        &               Asg)

!        IF (pid .EQ. master) THEN      
!            PRINT *, "Asg after scatter:"
!            PRINT *, "---------------------------"
!            DO i=1,Nys
!                DO j=1,Nxs
!                    WRITE(*,'(I2)', advance='no') Asg(i,j)
!                ENDDO
!                WRITE(*,*) ''
!            ENDDO
!        ENDIF

!        CALL Gather_Life_Rows(pid,                                   &
!        &                master,                                     &
!        &                rowtype,                                    &
!        &                Nxs,                                        & 
!        &                Nys,                                        &
!        &                Nxsp1,                                      &
!        &                Nysp1,                                      &
!        &                Asg,                                        &
!        &                Nx,                                         &
!        &                Ny,                                         &
!        &                A_k)
!
!        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!        IF (pid .EQ. master) THEN      
!            PRINT *, "A after gather:"
!            PRINT *, "---------------------------"
!            DO i=1,Ny
!                DO j=1,Nx
!                    WRITE(*,'(I2)', advance='no') A_k(i,j)
!                ENDDO
!                WRITE(*,*) ''
!            ENDDO
!        ENDIF


        ! [Loop time.]
        DO k=1,Nt 

            ! [Updated the ghost nodes on each partition. This includes the
            !  necessary ringed send and receives of the columns.]
            CALL UGN_Rows(pid,                                       & 
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

                    ! [Gether the life parent array for the write step.]
                    CALL Gather_Life_Rows(pid,                                   &
                    &                master,                                     &
                    &                rowtype,                                    &
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

        CALL MPI_TYPE_FREE(rowtype, ierr)
        DEALLOCATE(Asg)

    END SUBROUTINE Run_Row_Life

 
END MODULE rows_mod
