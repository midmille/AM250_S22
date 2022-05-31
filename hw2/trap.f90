! This file provides the solution to AM250 HW2 Prob1. 
! Using the trapezoidal rule to solve for the integral of two functions

PROGRAM Trap

  implicit none

  ! Parameters
  ! ----------
  integer, parameter                   :: dp=kind(0.d0)

  ! Local Variables
  ! ---------------
  ! [The lower and upper bound of integration.]
  real (dp)                            :: lbnd
  real (dp)                            :: ubnd
  ! [Number of grid points.]
  integer                              :: N 
  ! [The resulting integral.]
  real (dp)                            :: integ

  ! [Request upper bound]
  PRINT *, "Please enter integral lower bound:"
  READ *, lbnd
  PRINT *, ''

  ! [Request the lower bound.] 
  PRINT *, "Please enter integral upper bound:"
  READ *, ubnd
  PRINT *, ''

  ! [Read in the resolution.]
  PRINT *, "Please enter the number of grid points (N):"
  READ *, N
  PRINT *, ''

  ! [Init the integral val.]
  integ = 0.0_dp
  ! [Perform the integration for f1, x**2]
  CALL Int_Trap(ubnd, lbnd, N, integ)
  ! [Print the integral result to screen.]
  PRINT *, 'The resulting integral is:', integ

CONTAINS

  !--------------------------------------------------------
  REAL(dp) FUNCTION f(x)
  !--------------------------------------------------------
    ! This calculates x^2. 
    implicit None
    real (dp), intent(in)              :: x

    f = x**2
    !f = SIN(x)

  END FUNCTION 


  !--------------------------------------------------------
  SUBROUTINE Int_Trap(ubnd,                               &
  &                   lbnd,                               &
  &                   N,                                  & 
  &                   integ)
  !--------------------------------------------------------
    ! This routine solves for the integral of a function 
    ! within the given bounds using the trapezoidal rule 
    ! on a constant grid. 
    ! 
    ! Parameters
    ! ----------
    ! ubnd: Float
    !     The upper bound of integration.
    ! lbnd: Float
    !     The lower boudn of integration. 
    ! N: Integer
    !     The number of grid points. 
    ! 
    ! Returns
    ! -------
    ! integral: Float
    !     The resulting integral. 
    implicit none

    ! Parameters
    ! ----------
    ! [Upper and lower bounds.]
    real (dp), intent(in)              :: ubnd, lbnd
    ! [Number of grid points.]
    integer, intent(in)                :: N

    ! Returns
    ! -------
    ! [The resulting integral.]
    real (dp), intent(out)             :: integ

    ! Local Vars
    ! ----------
    ! [Indexing.]
    integer                            :: k 
    ! [The grid spacing.]
    real (dp)                          :: dx
    ! [The grid location.]
    real (dp)                          :: x
    ! [The values of at at different grid locations]
    real (dp)                          :: fkm1, fk


    ! [Calculate the grid spacing.] 
    dx = (ubnd - lbnd) / N

    ! [Init x to be lbnd.]
    x = lbnd
    ! [Init the integral to be zero.]
    integ = 0.0_dp
    ! [Loop over the points.]
    DO k=1,N
      ! [The function value at f(x_{k-1}).]
      fkm1 = f(x) 
      ! [Update x.]
      x = x + dx
      ! [The function value at f(x_{k}).]
      fk = f(x)
      ! [The trapezoidal rule for the integral solution.]
      integ = integ + (fkm1 + fk)
    ENDDO

    ! [Dividing sum by two and multiplying by dx.]
    integ = (dx/2)*integ

  END SUBROUTINE Int_Trap


END PROGRAM Trap

