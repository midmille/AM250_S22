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
        A(2,3) = 1

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
        character(len=100)                       :: iflag
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
    
        ! ROUTINE START
        ! =============
     
        ! [Allocate A.]
        ALLOCATE(A(Ny, Nx))
        A = 0

        ! [Initialize A according to the iflag.]
        IF (iflag .EQ. 'rand') THEN 
            PRINT *, "Initilizing Game of Life, initialization: Random"
            CALL Init_Rand_Life()
        ELSE IF (iflag .EQ. 'glide') THEN 
            PRINT *, "Initilizing Game of Life, initialization: Glider"
            CALL Init_Glide_Life()

    END SUBROUTINE Init_Life
    

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
        integer, intent(inout)                   :: A(0:Nxp1, 0:Nyp1)
        
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
        istp = Np1 - 1
        jstr = 1
        jstp = Np1 - 1

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
                    Atmp = 1
                ELSEIF (cnt .EQ. 2) THEN 
                    Atmp = A(i,j)
                ELSE
                    Atmp = 0
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

        ! [Open the file with tag=10.]
        OPEN(10, file=filename)

        ! [Loop over the rows.]
        DO i=1,Ny
            WRITE(10,*) (A(i,j), j=1,Nx)
        ENDDO

        ! [Close file.]
        CLOSE(10)

    END SUBROUTINE Write_Life_Mat


END MODULE funcs_mod
