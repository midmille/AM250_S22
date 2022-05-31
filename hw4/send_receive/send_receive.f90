! This file provides the solution toi AM 250 HW4 Prob2. 
! This is a simple send-receive mpi program. 

PROGRAM send_receive

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
  ! [The source and destination ids.]
  integer                              :: source, destination
  ! [The count of what is being sent.]
  integer, parameter                   :: count=5
  ! [The buffer being sent and received.]
  integer                              :: buffer(count)
  ! [The returned status, with size given by mpi status size.]
  integer                              :: stat(MPI_STATUS_SIZE)

  ! [Init MPI.]
  CALL MPI_INIT(ierr)
  ! [ID of current processor.]
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  ! [The number of processors.]
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  ! [The tag is the key for the send/receive process.]
  tag = 1234
  ! [This is the id of the source.]
  source = 0
  ! [Destination id.]
  destination = 1

  ! [The send compenent, ensure send processor is of id == 0.]
  IF (myid .EQ. source) THEN 
    ! [Set the thing we are sending.]
    buffer = (/1, 2, 3, 4, 5/)
    ! [Send the buffer using MPI_SEND.]
    CALL MPI_SEND(buffer,                                  &
    &             count,                                   &
    &             MPI_INTEGER,                             &
    &             destination,                             &
    &             tag,                                     &
    &             MPI_COMM_WORLD,                          & 
    &             ierr)
    ! [Print the details of the send.]
    PRINT *, ' Processor: ', myid, ' sent ', buffer
  ENDIF 

  ! [The receive component, ensure the destination processor is id==1.]
  IF (myid .EQ. destination) THEN 
    ! [Receive the sent buffer.]
    CALL MPI_RECV(buffer,                                  & 
    &             count,                                   & 
    &             MPI_INTEGER,                             & 
    &             source,                                  &
    &             tag,                                     &
    &             MPI_COMM_WORLD,                          &
    &             stat,                                    & 
    &             ierr)
    ! [Print the details of the receive. ]
    PRINT *, ' Processor: ', myid, ' sent ', buffer
  ENDIF 

  ! [Terminate MPI.]
  CALL MPI_FINALIZE(ierr)

END PROGRAM 
