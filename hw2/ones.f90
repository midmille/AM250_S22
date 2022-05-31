! This file provides the solution to AM250 HW2 Prob2. 
! This involves the evolution of a ones array follwoing the rules provided by
! the homework.

PROGRAM Ones

  implicit none 

  ! Parameters
  ! ----------
  integer, parameter                   :: dp=kind(0.d0)

  ! Local Variables
  ! ---------------
  ! [The size of the square matrix.]
  integer                              :: N, Np1
  ! [The random ones array.]
  integer, allocatable               :: A(:,:)
  ! [The resulting array after the ones algorithm.] 
  integer, allocatable               :: B(:,:)
  ! [The number of ones around an entry.]
  integer                              :: cnt
  ! [Indexing.]
  integer                              :: i,j,k

  ! [Request the dimension of the array.]
  PRINT *, "Please enter matrix dimesion, N: "
  READ *, N
  PRINT *, ''
    
  ! [N plus one.]
  Np1 = N+1

  ! [Allocate the array.]
  ! [The indexing in A is constructed such that it includes the halo.]
  ALLOCATE(A(0:Np1,0:Np1))
  ALLOCATE(B(N,N))
  ! [Init the array.]
  A = 0.0_dp
  B = 0.0_dp

  ! [Init the random sequence.]
  CALL SRAND(2)
  ! [Fill the indside of A with random ones.]
  DO i=1,N
    DO j=1,N
      ! [RAND() produces real between 0,1]
      ! [INT() make the number 0 if < 1, and 1 if >1 but <2.]
      A(i,j) = INT(RAND() + 0.5_dp)
    ENDDO
  ENDDO

  ! [Now construct the halo.]
  ! [Top and bottom rows.]
  A(0, :) = A(N,:)
  A(Np1, :) = A(1,:)
  ! [Left and right columns.]
  A(:, 0) = A(:,N)
  A(:, Np1) = A(:,1)
  ! [The corners.]
  A(0,0) = A(N,N)
  A(Np1, 0) = A(0, N)
  A(Np1, Np1) = A(1,1)
  A(0, Np1) = A(N, 0)

  ! [Implement the ones algorithm on A but store solution to B.]
  DO i=1,N 
    DO j=1,N
      ! [Init the ones counter to zero.]
      cnt = 0 
      ! [Count the number of ones around the given point.]
      ! [Wrapping around clockwise from top left corner.]
      ! [Total of 8 entries being added to counter.]
      cnt = cnt + A(i+1, j-1)
      cnt = cnt + A(i+1, j)
      cnt = cnt + A(i+1, j+1)
      cnt = cnt + A(i, j+1)
      cnt = cnt + A(i-1, j+1)
      cnt = cnt + A(i-1, j)
      cnt = cnt + A(i-1, j-1)
      cnt = cnt + A(i, j-1)

      ! [Check if the counter is 3.]
      IF (cnt .EQ. 3) THEN 
        ! [Make the center entry a one.]
        B(i,j) = 1
      ENDIF
    ENDDO
  ENDDO

  ! [Now print the original matrix.]
  PRINT *, "The original random ones matrix:"
  DO i=1,N
    WRITE(*,*) (A(i,j), j=1, N)
  ENDDO

  ! [Print the resulting matrix.]
  PRINT *, "The resulting matrix after ones algorithm:"
  DO i=1,N
    WRITE(*,*) (B(i,j), j=1, N)
  ENDDO


END PROGRAM Ones
