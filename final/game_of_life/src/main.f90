! This program is the main driver routine for the game of life.

PROGRAM main

    USE mpi

    USE funcs_mod, ONLY: Print_Thing1
    USE parallel_cols_mod, ONLY: Print_Thing2

    implicit none
    
    ! Problem Parameters
    ! ------------------
    ! [Double Precision for real numbers.]
    integer, parameter                 :: dp=kind(0.d0)
    ! [The size of the x-dimension, the number of columns.]
    integer, parameter                 :: Nx = 100
    ! [The size of the y-dimension, the number of rows.]
    integer, parameter                 :: Ny = 100
    ! [The number of time steps to take.]
    integer, parameter                 :: Nt = 20
    ! [The number of time steps between each gather and write.]
    integer, parameter                 :: Nw = 5
    ! [The parallelization flag, the options are, 
    !  'serial' : Runs the problem in series on master--No Parallel,
    !  'cols'   : Parallelizes the problem into column partitions, 
    !  'rows'   : Parallelizes the problem into rows paritions,
    !  'tile'   : Parallelizes the problem into rectangular tile partitions.]
    character(len=100), parameter      :: pflag = "cols"
    ! [The Initialization flag, the options are, 
    !  'rand'  : Initializes the domain to be randomly alive or dead, 
    !  'glide' : Initializes the domain to have glider formation in top left.]
    character(len=100), parameter      :: iflag = 'rand'
    ! [The header of the save file, format will be '{savefile_head}_{k}.dat'
    !  where k is the given time step for the save.]
    character(len=100), parameter      :: savefile_head = "life_out"
    ! [The out directory for the serial implementation.]
    character(len=100), parameter      :: serial_outdir = "output/serial_out/"
    ! [The out directory for the column parallelization.]
    character(len=100), parameter      :: cols_outdir = "output/cols_out/"

    ! Problem Variables
    ! -----------------
    ! [The parent matrix, the entire domain.]
    integer, allocatable               :: A(:,:)


    ! MPI Parameters
    ! --------------
    ! [The master processor id is zero.]
    integer, parameter                  :: master=0

    ! MPI Parameters
    ! --------------
    ! [The MPI error return code.]
    integer                            :: ierr
    ! [The processor id.]
    integer                            :: pid
    ! [The total number of available processors.]
    integer                            :: numprocs



    ! PROGRAM START
    ! =============

    ! [Init MPI.]
    CALL MPI_INIT(ierr)
    ! [ID of current processor.]
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
    ! [The number of processors.]
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    ! [This creates and allocates the parent matrix on master processor.]
    CALL Init_Life()

    ! [BARRIER??????

    ! [Run the life model for given parallelization flag.]
    IF (pflag .EQ. 'none') THEN  
        CALL Run_Serial_Life()
    ELSE IF (pflag .EQ. 'cols') THEN 
        CALL Run_Column_Life()
    ELSE IF (pflag .EQ. 'rows') THEN 
        CALL Run_Row_Life()
    ELSE IF (pflag .EQ. 'tile') THEN
        CALL Run_Tile_Life()
    END IF

    ! [This deallocates the matrix A on master.]
    CALL End_Life()

    ! [Terminate MPI.]
    CALL MPI_FINALIZE(ierr)

END PROGRAM main
