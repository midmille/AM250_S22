! This files provides the solution to AM 250 HW4 Problem 6. 
! It implements the parallel calculation of pi.

PROGRAM pi

  USE mpi

  implicit none 

  ! Parameters
  ! ----------
  integer, parameter                   :: dp=kind(0.d0)
  ! [The count of what is being sent.]
  integer, parameter                   :: cnt=1
  ! [The number of shots.]
  integer, parameter                   :: N = 1000000

  ! Local Variables
  ! ---------------
  ! [The mpi error code variable.]
  integer                              :: ierr
  ! [The id of the individual processesor.]
  integer                              :: myid
  ! [The total number of processors.]
  integer                              :: numprocs
  ! [The root/main processor.]
  integer                              :: root
  ! [Indexing.]
  integer                              :: k
  ! [The shots per processor. It is necesary that the N_shots/numprocs == INT]
  integer                              :: N_pp
  ! [The radius, x, and y]
  real (dp)                            :: x, y, r
  ! [The number of shots in the circle total and per processor respectively.]
  integer                              :: N_intot, N_inpp
  ! [The resulting value for pi.]
  real (dp)                            :: pi_res

  ! [Init MPI.]
  CALL MPI_INIT(ierr)
  ! [ID of current processor.]
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  ! [The number of processors.]
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  ! [The root processor.]
  root = 0

  ! [Neccessary  that N is a multiple of numprocs.]
  N_pp = N / numprocs
  
  ! [Use the processor id for the seed of the random number s.t. each processor
  !  actually uses different numbers]
  CALL SRAND(myid+1)

  ! [Init the number of shots in circle to be 0.]
  N_inpp = 0

  ! [Loop over the N_pp.]
  DO k=1,N_pp
    ! [The x and y coordinates of the circle.]
    ! [Since the radius function is an even function, 
    !  it does not matter that x and y are exclusively between 0 and 1.]
    x = RAND()
    y = RAND() 
    
    ! [The radius of the circle.]
    r = SQRT(x**2 + y**2)

    ! [If the shots lands within the circle.]
    IF (r .LT. 1) THEN 
      N_inpp = N_inpp + 1
    END IF
  END DO 

  ! [Print the processor and their N_inpp.]
  PRINT *, 'Processor: ', myid, ' N_inpp: ', N_inpp
        
  ! [Init the N_intot to be 0.]
  N_intot = 0

  ! [Synchronize everything.]
!  CALL MPI_BARRIER(MPI_COMM_WORLD)

  ! [Gather and sum the N_inpp.]
  CALL MPI_REDUCE(N_inpp,                                  &
  &               N_intot,                                 &
  &               cnt,                                     &
  &               MPI_INTEGER,                             &
  &               MPI_SUM,                                 &
  &               root,                                    & 
  &               MPI_COMM_WORLD,                          &
  &               ierr)

  ! [Final computation and printing on root.]
  IF (myid .EQ. root) THEN 
    ! [Calculate pi as a ratio of shots in and out.]
    pi_res = 4 * (DBLE(N_intot)/DBLE(N))

    ! [Print result.]
    PRINT *, ' Final Ratio: ', N_intot, '/', N 
    PRINT *, ' pi is estimated to be: ', pi_res
  END IF 

  ! [Terminate MPI.]
  CALL MPI_FINALIZE(ierr)  

END PROGRAM 
