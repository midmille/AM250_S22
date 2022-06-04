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
    SUBROUTINE Life(Nx,                                              &
    &               Ny,                                              &
    &               A)                                            
    !-----------------------------------------------------------------

END MODULE funcs_mod
