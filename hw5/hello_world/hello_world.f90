! This file is to implement a hello world algortihm for homework 5 AM
! 250. 

PROGRAM hello_world

    implicit none 

    ! Parameters
    ! ----------
    integer, parameter                 :: dp=kind(0.d0)

    ! Local Variables
    ! ---------------
    integer                            :: NTHREADS 
    integer                            :: TID
    integer                            :: OMP_GET_NUM_THREADS 
    integer                            :: OMP_GET_THREAD_NUM

! [Fork a team of threads by giving them their own copies of variables. ]
!$OMP PARALLEL PRIVATE(NTHREADS, TID)

    ! [Get the thread number.]
    TID = OMP_GET_THREAD_NUM()
    PRINT *, 'Hello World from thread = ', TID 

    ! [This is done by the master ndoe.]
    IF (TID .EQ. 0) THEN 
        NTHREADS = OMP_GET_NUM_THREADS()
        PRINT *, 'Number of threads = ', NTHREADS 
    END IF

! [Join with master and end parallel.]
!$OMP END PARALLEL

END PROGRAM



