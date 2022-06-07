! This file provides the solution toi AM 250 HW4 Prob3. 
! This is a ping pong based send receive mpi program. 

PROGRAM latency

    USE mpi
    
    implicit none 
  
    ! Parameters
    ! ----------
    integer, parameter                   :: dp=kind(0.d0)
    ! [The number of times to loop the ping pong.]
    integer, parameter                   :: N=10000
    ! [The name of the input file.]
    character(len=100), parameter        :: in_file="input.txt"
    ! [The name of the output file.]
    character(len=100), parameter        :: out_file="output.txt"

  
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
    integer                              :: cnt
!    integer, parameter                   :: cnt=5
    ! [The buffer being sent and received.]
    integer, allocatable                 :: buffer(:)
    ! [Indexing.]
    integer                              :: k
    ! [The timing variables.]
    real (dp)                            :: t1, t2, t
    ! [The returned status, with size given by mpi status size.]
    integer                              :: stat(MPI_STATUS_SIZE)
    integer                              :: statuses(2,MPI_STATUS_SIZE)


    ! [Read in the size of the buffer message to be sent.]
    OPEN(10, file=in_file)
    READ(10,*) cnt
    CLOSE(10)
  
    ! [Allocate the buffer and Init zeros on both processors.]
    ALLOCATE(buffer(cnt))
    buffer = 0.0

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

    ! [Init the pin pong buffer on the main processor before looping the
    !  messages.]
    IF (myid .EQ. source) THEN 
        ! [The source and destination for this processor.]
        ! [Set the thing we are sending.]
        buffer = 1
    ENDIF
  
    ! [Start the mpi timer.]
    t1 = MPI_WTIME()
    DO k=1,N
        ! [The send compenent, ensure send processor is of id == 0.]
        IF (myid .EQ. source) THEN 
     
            CALL MPI_SEND(buffer,                                  &
            &              cnt,                                    &
            &              MPI_INTEGER,                            &
            &              destination,                            &
            &              tag,                                    &
            &              MPI_COMM_WORLD,                         & 
            !              [OUTPUT]
            &              ierr)
     
            ! [Print the details of the send.]
!            PRINT *, ' Processor: ', myid, ' sent ', buffer
     
            CALL MPI_RECV(buffer,                                  &
            &              cnt,                                    &
            &              MPI_INTEGER,                            &
            &              destination,                            &
            &              tag,                                    &
            &              MPI_COMM_WORLD,                         & 
            !              [OUTPUT]
            &              stat,                                   &
            &              ierr)
     
        ! [Print the details of the send.]
!            PRINT *, ' Processor: ', myid, ' recv ', buffer
     
        ENDIF 
     
        ! [The receive component, ensure the destination processor is id==1.]
        IF (myid .EQ. destination) THEN 
            CALL MPI_RECV(buffer,                                  &
            &              cnt,                                    &
            &              MPI_INTEGER,                            &
            &              source,                                 &
            &              tag,                                    &
            &              MPI_COMM_WORLD,                         & 
            !              [OUTPUT]
            &              stat,                                   &
            &              ierr)
     
            ! [Print the details of the send.]
!            PRINT *, ' Processor: ', myid, ' recv ', buffer
     
            !  [Edit the buffer slightly.]
!            buffer = buffer + 1.d0
     
            CALL MPI_SEND(buffer,                                  &
            &              cnt,                                    &
            &              MPI_INTEGER,                            &
            &              source,                                 &
            &              tag,                                    &
            &              MPI_COMM_WORLD,                         & 
            !              [OUTPUT]
            &              ierr)
     
            ! [Print the details of the send.]
!            PRINT *, ' Processor: ', myid, ' sent ', buffer
  
        ENDIF 
    ENDDO

    ! [End the mpi timer.]
    t2 = MPI_WTIME()
    ! [The total elapsed time.]
    t = t2 - t1
    
    IF (myid .EQ. source) THEN 
        ! [Print the result to job output file.]
        PRINT *, "Elapsed time: ", t
        ! [Write the elapsed time to output file.]
        OPEN(10, file=out_file)
        WRITE(10,*) t
        CLOSE(10)
    ENDIF 
  
    ! [Terminate MPI.]
    CALL MPI_FINALIZE(ierr)
  
    ! [Deallocate the buffer.]
    DEALLOCATE(buffer)

END PROGRAM 
  
  
