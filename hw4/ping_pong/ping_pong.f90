! This file provides the solution toi AM 250 HW4 Prob3. 
! This is a ping pong based send receive mpi program. 

PROGRAM ping_pong

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
  ! [The tag for this send/receive process.]
  integer                              :: tag
  ! [The request parameter for the immediate send receive. The array
  !  is of length 2 since there is the request for the send and receive 
  !  that both must be received.]
  integer                              :: request1, request2
  integer                              :: requests(2)
  ! [The source and destination ids.]
  integer                              :: source, destination
  ! [The count of what is being sent.]
  integer, parameter                   :: count=5
  ! [The buffer being sent and received.]
  integer                              :: buffer(count)
  ! [The returned status, with size given by mpi status size.]
  integer                              :: stat(MPI_STATUS_SIZE)
  integer                              :: statuses(2,MPI_STATUS_SIZE)

  ! [Init MPI.]
  CALL MPI_INIT(ierr)
  ! [ID of current processor.]
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  ! [The number of processors.]
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  ! [The tag for this ping pong process.]
  tag = 1234

  ! [The source and destination ids.]
  source = 0 
  destination = 1

  ! [The send compenent, ensure send processor is of id == 0.]
  IF (myid .EQ. source) THEN 
    ! [The source and destination for this processor.]
    ! [Set the thing we are sending.]
    buffer = (/1, 2, 3, 4, 5/)

    CALL MPI_SEND(buffer,                                  &
    &              count,                                  &
    &              MPI_INTEGER,                            &
    &              destination,                            &
    &              tag,                                    &
    &              MPI_COMM_WORLD,                         & 
    !              [OUTPUT]
    &              stat,                                   & 
    &              ierr)

    ! [Print the details of the send.]
    PRINT *, ' Processor: ', myid, ' sent ', buffer

    CALL MPI_RECV(buffer,                                  &
    &              count,                                  &
    &              MPI_INTEGER,                            &
    &              destination,                            &
    &              tag,                                    &
    &              MPI_COMM_WORLD,                         & 
    !              [OUTPUT]
    &              stat,                                   &
    &              ierr)

   ! [Print the details of the send.]
    PRINT *, ' Processor: ', myid, ' recv ', buffer

  ENDIF 

  ! [The receive component, ensure the destination processor is id==1.]
  IF (myid .EQ. destination) THEN 
    CALL MPI_RECV(buffer,                                  &
    &              count,                                  &
    &              MPI_INTEGER,                            &
    &              source,                                 &
    &              tag,                                    &
    &              MPI_COMM_WORLD,                         & 
    !              [OUTPUT]
    &              stat,                                   &
    &              ierr)

    ! [Print the details of the send.]
    PRINT *, ' Processor: ', myid, ' recv ', buffer

    !  [Edit the buffer slightly.]
    buffer = buffer + 1.d0

    CALL MPI_SEND(buffer,                                  &
    &              count,                                  &
    &              MPI_INTEGER,                            &
    &              source,                                 &
    &              tag,                                    &
    &              MPI_COMM_WORLD,                         & 
    !              [OUTPUT]
    &              stat,                                   &
    &              ierr)

    ! [Print the details of the send.]
    PRINT *, ' Processor: ', myid, ' sent ', buffer

  ENDIF 

  ! [Terminate MPI.]
  CALL MPI_FINALIZE(ierr)

END PROGRAM 


