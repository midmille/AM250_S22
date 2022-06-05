! AM250_S22/final/game_of_life/src/funcs_mod.f90
! Author : Miles D. Miller, University of California Santa Cruz
! Created : 2:02 PM, June 3, 2022
! About : This is the general purpose module for the global routines used in the
!         game of life algorithm. 

MODULE funcs_mod

    USE mpi

    implicit none

    ! Paramters
    ! ---------
    ! [Double.]
    integer, parameter                 :: dp = kind(0.d0)

CONTAINS

    
    !-----------------------------------------------------------------
    SUBROUTINE Init_Rand_Life(Nx,                                    &
    &                         Ny,                                    &
    &                         A)
    !-----------------------------------------------------------------
        ! This routine creates a random ones or zero initialization for
        ! the provided matrix A. 

        ! Parameters
        ! ----------
        ! [The number of columns.]
        integer, intent(in)                      :: Nx
        ! [The number of rows.]
        integer, intent(in)                      :: Ny

        ! Returns
        ! -------
        integer, intent(inout)                   :: A(Ny,Nx)

        ! Variables
        ! ---------
        ! [The temporary random real array.]
        real (dp)                                :: R(Ny,Nx)
        ! [Indexing.]
        integer                                  :: i,j

        ! ROUTINE START
        ! =============

        ! [First initialize the temp array to be random reals 0<R(i,j)<1.
        CALL RANDOM_NUMBER(R)

        ! [Loop over array entries.]
        DO j=1,Nx
            DO i=1,Ny
                ! [If less than 0.5, the zero.]
                IF (R(i,j) < 0.5) THEN 
                    A(i,j) = 0
                ! [If greater than 0.5, thenm one.]
                ELSEIF (R(i,j) > 0.5) THEN 
                    A(i,j) = 1
                ENDIF
            ENDDO 
        ENDDO

    END SUBROUTINE Init_Rand_Life


    !-----------------------------------------------------------------
    SUBROUTINE Init_Glide_Life(Nx,                                   &
    &                          Ny,                                   &
    &                          A)
    !-----------------------------------------------------------------
        ! This routine sets the glider formation as the game of life
        ! inititalization for the matrix A. 
        ! It follows the glider form given by the final project report
        ! question 3.

        ! Parameters
        ! ----------
        ! [The number of columns.]
        integer, intent(in)                      :: Nx
        ! [The number of rows.]
        integer, intent(in)                      :: Ny

        ! Returns
        ! -------
        integer, intent(inout)                   :: A(Ny,Nx)

        ! ROUTINE START
        ! =============

        ! [Set A to zero again.]
        A = 0

        ! [Set the glider formation.]
        A(2,1) = 1
        A(3,2) = 1
        A(1,3) = 1
        A(2,3) = 1 
        A(3,3) = 1

    END SUBROUTINE Init_Glide_Life


    !-----------------------------------------------------------------
    SUBROUTINE Init_Life(pid,                                        &
    &                    master,                                     & 
    &                    iflag,                                      & 
    &                    Nx,                                         &
    &                    Ny,                                         &
    &                    A) 
    !-----------------------------------------------------------------
    
        ! This routine is for the initialization of the life matrix on 
        ! the master processor. 
        ! The routine allocates the array, initilizes its state accoriding to 
        ! The designation of iflag, and then returns the initialized A.
 
        ! Parameters
        ! ----------
        ! [The processdor id.]
        integer, intent(in)                      :: pid
        ! [The master id.]
        integer, intent(in)                      :: master
        ! [The initialization flag.]
        character(len=*)                         :: iflag
        ! [The x-dimension.]
        integer                                  :: Nx
        ! [The y-dimension.]
        integer                                  :: Ny

        ! Returns
        ! -------
        ! [The allocatable parent life domain matrix.]
        integer, allocatable, intent(inout)      :: A(:,:)

        ! Variables
        ! ---------
        ! None
    
        ! ROUTINE START
        ! =============
     
!        IF (pid .EQ. master) THEN 
            ! [Allocate A.]
            ALLOCATE(A(Ny, Nx))
            A = 0

            ! [Initialize A according to the iflag.]
            IF (iflag .EQ. 'rand') THEN 
                PRINT *, "Initilizing Game of Life, initialization: Random"
                CALL Init_Rand_Life(Nx, Ny, A)
            ELSE IF (iflag .EQ. 'glide') THEN 
                PRINT *, "Initilizing Game of Life, initialization: Glider"
                CALL Init_Glide_Life(Nx, Ny, A)
            ENDIF
!        ENDIF

    END SUBROUTINE Init_Life
    

    !-----------------------------------------------------------------
    SUBROUTINE End_Life(pid,                                         &
    &                   master,                                      & 
    &                   A) 
    !-----------------------------------------------------------------
    
        ! This routine is for the initialization of the life matrix on 
        ! the master processor. 
        ! The routine allocates the array, initilizes its state accoriding to 
        ! The designation of iflag, and then returns the initialized A.
 
        ! Parameters
        ! ----------
        ! [The processdor id.]
        integer, intent(in)                      :: pid
        ! [The master id.]
        integer, intent(in)                      :: master
        ! [The allocatable parent life domain matrix.]
        integer, allocatable, intent(inout)      :: A(:,:)

        ! Variables
        ! ---------
        ! None
    
        ! ROUTINE START
        ! =============
     
!        IF (pid .EQ. master) THEN
            DEALLOCATE(A)
!        ENDIF

    END SUBROUTINE End_Life

 
    !-----------------------------------------------------------------
    SUBROUTINE Life(Nxp1,                                            &
    &               Nyp1,                                            &
    &               A)                                            
    !-----------------------------------------------------------------
        ! This routine performs the game of life algorithm on the provided 
        ! initialized matrix. 
        ! This algorithm assumes that ghost noes are already in play.

        ! Parameters
        ! ----------
        ! [The interior y-dim plus one..]
        integer, intent(in)                      :: Nxp1
        ! [The interior y-dim plus one.]
        integer, intent(in)                      :: Nyp1
        
        ! Returns
        ! -------
        ! [The evolved game of life matrix.]
        integer, intent(inout)                   :: A(0:Nyp1, 0:Nxp1)
        
        ! Variables
        ! ---------
        ! [The temporary A, such that we don't overwrite cells.]
        integer                                  :: Atmp(0:Nyp1, 0:Nxp1)
        ! [The count of alive cells.]
        integer                                  :: cnt
        ! [Indexing.]
        integer                                  :: i,j
        ! [The start and stop of the  i and j loops.]
        integer                                  :: istr, istp
        integer                                  :: jstr, jstp

        ! ROUTINE START
        ! =============
         
        ! [Init Atmp as zeros.]
        Atmp = 0

        ! [Set the loop start and stop variables.]
        istr = 1
        istp = Nyp1 - 1
        jstr = 1
        jstp = Nxp1 - 1

        ! [Loop over the entries of the provided matrix.]
        DO j=jstr,jstp
            DO i=istr, istp
                ! [Init the cnt to 0.]
                cnt = 0
                ! [Count the alive cells by wrapping around the given point.]
                cnt = cnt + A(i+1, j-1)
                cnt = cnt + A(i+1, j)
                cnt = cnt + A(i+1, j+1)
                cnt = cnt + A(i, j+1)
                cnt = cnt + A(i-1, j+1)
                cnt = cnt + A(i-1, j)
                cnt = cnt + A(i-1, j-1)
                cnt = cnt + A(i, j-1)

                ! [Rules for alive or dead.]
                IF (cnt .EQ. 3) THEN 
                    Atmp(i,j) = 1
                ELSEIF (cnt .EQ. 2) THEN 
                    Atmp(i,j) = A(i,j)
                ELSE
                    Atmp(i,j) = 0
                ENDIF
            ENDDO
        ENDDO 

        ! [Overwrite A as the new evolved A.]
        A = Atmp
    
    END SUBROUTINE Life


    !----------------------------------------------------------------
    SUBROUTINE Write_Life_Mat(filename,                             &
    &                         Nx,                                   &
    &                         Ny,                                   & 
    &                         A)
    !----------------------------------------------------------------
        ! This routine write the provided matrix to the designated 
        ! file.

        ! Parameters
        ! ----------
        ! [The name of the file to write the matrix to.]
        character(len=*), intent(in)             :: filename
        ! [The number of columns.]
        integer, intent(in)                      :: Nx
        ! [The number of rows.]
        integer, intent(in)                      :: Ny
        ! [The matrix to write.]
        integer, intent(in)                      :: A(Ny, Nx)

        ! Returns
        ! -------
        ! None

        ! Variables
        ! ---------
        ! Indexing
        integer                                  :: i,j

        ! ROUTINE START
        ! =============
        
        ! [Open the file with tag=10.]
        OPEN(10, file=filename)

        ! [Loop over the rows.]
        DO i=1,Ny
            DO j=1,Nx
                WRITE(10,'(I2)', advance='no') A(i,j)
            ENDDO
            WRITE(10,*) ''
        ENDDO

        ! [Close file.]
        CLOSE(10)

    END SUBROUTINE Write_Life_Mat


    !------------------------------------------------------------
    SUBROUTINE Write_Life_Step(k,                               & 
    &                          Nx,                              & 
    &                          Ny,                              & 
    &                          A,                               & 
    &                          savefile_head,                   &
    &                          outdir)
    !------------------------------------------------------------
        ! This algorithm writes the matrix A to file for the current
        ! time step.

        ! Parameters
        ! ----------
        ! [The current time step.] 
        integer, intent(in)                      :: k 
        ! [The number of columns.]
        integer, intent(in)                      :: Nx
        ! [The number of rows.]
        integer, intent(in)                      :: Ny
        ! [The matrix at step k.]
        integer, intent(in)                      :: A(Ny, Nx)
        ! [The savefile header.] 
        character(len=*), intent(in)             :: savefile_head
        ! [The output directory for the parallelization method.]
        character(len=*), intent(in)             :: outdir

        ! Returns
        ! -------
        ! None

        ! Variables
        ! ---------
        ! [The character variable for k.]
        character(len=6)                         :: k_char
        ! [The savefile name.]
        character(len=100)                       :: savefile


        ! ROUTINE START
        ! =============

        ! [Form the savefile name for the given time step.]
        ! [Write the step index into a character.]
        WRITE(k_char, '(I6)'), 100000 + k
        ! [Concatinate the strings.]
        savefile = outdir // savefile_head // k_char // ".dat"
        
        PRINT *, savefile
        ! [Write the matrix to the above file.]
        CALL Write_Life_Mat(savefile, Nx, Ny, A)

    END SUBROUTINE Write_Life_Step


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

    END SUBROUTINE Partition_Cols


    !----------------------------------------------------------------
    SUBROUTINE Partition_Rows(Ny,                                   & 
    &                         Np,                                   &
    &                         Nys)
    !----------------------------------------------------------------
        ! This routine is for the partitioning of the domain to different 
        ! processors. 
        ! It checks for the cohesion of the desired number of processors
        ! to the size of the domain. 
        ! Since the partition is for columns we only have to work with the
        ! number of columns.

        ! Parameters
        ! ----------
        ! [The number of rows.] 
        integer, intent(in)                      :: Ny
        ! [The number of processors.]
        integer, intent(in)                      :: Np

        ! Returns
        ! -------
        ! [The number of rows per processor.]  
        integer, intent(out)                     :: Nys

        ! Variables
        ! ---------
        ! None 

        ! ROUTINE START
        ! =============
        
        ! [Require that processor count is less than or equal to number of rows.]
        IF (Np .GT. Ny) THEN 
            PRINT *, "The number of rows is less than number of processors"
            STOP "Error: Bad row number and processor count cohesion."
        ENDIF

        ! [Require that the number of rows be evenly divisible by the number of processor.]
        IF (MOD(Ny, Np) .NE. 0) THEN 
            PRINT *, "For row partition, please ensure MOD(Ny,Np) == 0"
            STOP "Error: Bad row number and processor count cohesion."
        ENDIF

        ! [This gives the number of rows per parttition.] 
        Nys = Ny / Np

    END SUBROUTINE 


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

!        ! [Make A flat such that scatter is contiguous]
!        A_flt = RESHAPE(A, (/Nx*Ny/))

        ! [Scatter]
        CALL MPI_SCATTER(A,                                          &
        &                send_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                Asg(1:Nys,1:Nxs),                           &
        &                recv_cnt,                                   &
        &                MPI_INTEGER,                                &
        &                master,                                     &
        &                MPI_COMM_WORLD,                             &
        &                ierr)

!        Asg(1:Nys,1:Nxs) = RESHAPE(As_flt, (/Nys, Nxs/))

    
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
!        As_flt = RESHAPE(Asg(1:Nys, 1:Nxs), (/Nys*Nxs/))

        ! [Scatter]
        CALL MPI_GATHER(Asg(1:Nys,1:Nxs),                            &
        &               send_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               A,                                           &
        &               recv_cnt,                                    &
        &               MPI_INTEGER,                                 &
        &               master,                                      &
        &               MPI_COMM_WORLD,                              &
        &               ierr)

!        A = RESHAPE(A_flt, (/Ny, Nx/))

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


END MODULE funcs_mod
