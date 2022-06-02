! This file provides the solution to AM 250 HW4 Problem 5. 
! This is the ring based rend and receive algorithm. 

PROGRAM ring 

  USE mpi

  implicit none

  ! Parameters
  ! ----------
  integer, parameter                   :: dp=kind(0.d0)
  ! [The count of what is being sent.]
  integer, parameter                   :: count=1

  ! Local Variables
  ! ---------------
  ! [The mpi error code variable.]
  integer                              :: ierr
  ! [The id of the individual processesor.]
  integer                              :: myid
  ! [The ids of the processors to the right and left.]
  integer                              :: left, right
  ! [The total number of processors.]
  integer                              :: numprocs
  ! [The tag for this send/receive process.]
  integer                              :: tag
  ! [The request number for each ISEND, IRECV]
  integer                              :: request_send, request_recv
  ! [The request array. ]
  integer                              :: requests(2)
  ! [The source and destination ids.]
  integer                              :: source, destination
  ! [The buffer being sent and received.]
  integer                              :: buffer
  ! [The returned status, with size given by mpi status size.]
  integer                              :: stat(MPI_STATUS_SIZE)
  ! [The returned statuses array, the 2 is for the 2 requests..]
  integer                              :: statuses(MPI_STATUS_SIZE, 2)
  ! [Indexing.]
  integer                              :: k

  ! [Init MPI.]
  CALL MPI_INIT(ierr)
  ! [ID of current processor.]
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  ! [The number of processors.]
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  tag = 1234

  ! [Init the buffer data is just the processor id for now.]
  buffer = myid

  IF (myid .EQ. 0) THEN
    ! [left is the id of the last processor in the group.]
    left = numprocs - 1
    right = myid + 1
  ELSE IF (myid .EQ. numprocs-1) THEN 
    left = myid - 1
    ! [The right processor is just the first processor.]
    right = 0
  ELSE 
    left = myid - 1
    right = myid + 1 
  END IF

  ! [Loop over the ring data transfers.]
  DO k=1,numprocs

    ! [Send the data to the right.] 
    CALL MPI_ISEND(buffer,                                 &
    &              count,                                  &
    &              MPI_INTEGER,                            &
    &              right,                                  &
    &              tag,                                    &
    &              MPI_COMM_WORLD,                         &
    !              [OUTPUT]                                ! 
    &              request_send,                           &
    &              ierr)
    ! [Receive data from the left.]
    CALL MPI_IRECV(buffer,                                 &
    &              count,                                  &
    &              MPI_INTEGER,                            &
    &              left,                                   &
    &              tag,                                    &
    &              MPI_COMM_WORLD,                         &
    !              [OUTPUT]                                ! 
    &              request_recv,                           &
    &              ierr)

    ! [The request array.]
    requests = (/request_send, request_recv/)

    ! [Wait for the recv and send to compelete.]
    CALL MPI_WAITALL(2, requests, statuses, ierr)
    
    PRINT *, 'Iteration: ', k, ' Processor: ', myid, ' Buffer: ', buffer 

  END DO 

  ! [Terminate MPI.]
  CALL MPI_FINALIZE(ierr)  

END PROGRAM 
  
