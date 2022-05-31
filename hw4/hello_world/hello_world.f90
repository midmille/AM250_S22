! This files provides the solution to AM 250 HW4 Prob 1. 
! The Hello World Algorithm. 

PROGRAM hello_world
    
  USE mpi

  implicit none

  ! Parameters
  ! ----------
  integer, parameter                   :: dp=kind(0.d0)

  ! Local Variables
  ! ---------------
  ! [The mpi error code variable.]
  integer                              :: ierr
  ! [The id of the individual processesor.]
  integer                              :: myid
  ! [The total number of processors.]
  integer                              :: numprocs

  ! [Initialize  MPI.]
  CALL MPI_INIT(ierr)

  ! [Retreive the ID of the current processor.]
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

  ! [Retrieve the number of processors.]
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  ! [Output a hellow from each processor.]
  PRINT *, "Hello world from processor ", myid, " out of ", numprocs

  ! [Finalize the MPI program.]
  CALL MPI_FINALIZE(ierr)


END PROGRAM hello_world
