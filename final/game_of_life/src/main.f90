! AM250_S22/final/game_of_life/src/main.f90
! Author : Miles D. Miller, University of California Santa Cruz
! Created : 12:30 PM, June 3, 2022
! About : This is the main driver program for the game of life final 
!         project for AM 250. 

PROGRAM main

    USE mpi

    USE funcs_mod, ONLY: Init_Life, End_Life
    USE serial_mod, ONLY: Run_Serial_Life
    USE cols_mod, ONLY: Run_Column_Life
    USE rows_mod, ONLY: Run_Row_Life

    implicit none
    
    ! Problem Parameters
    ! ------------------
    ! [Double Precision for real numbers.]
    integer, parameter                 :: dp=kind(0.d0)
    ! [The size of the x-dimension, the number of columns.]
    ! [The input file for Nx.] 
    character(len=15), parameter       :: Nx_infile ='input/Nx_in.dat'
!    integer, parameter                 :: Nx = 10
    ! [The size of the y-dimension, the number of rows.]
    ! [The input file for Ny.] 
    character(len=15), parameter       :: Ny_infile = 'input/Ny_in.dat'
!    integer, parameter                 :: Ny = 10
    ! [The number of time steps to take.]
    integer, parameter                 :: Nt = 80
    ! [The number of time steps between each gather and write.]
    integer, parameter                 :: Nw = 1
    ! [The parallelization flag, the options are, 
    !  'serial' : Runs the problem in series on master--No Parallel,
    !  'cols'   : Parallelizes the problem into column partitions, 
    !  'rows'   : Parallelizes the problem into rows paritions,
    !  'tile'   : Parallelizes the problem into rectangular tile partitions.]
!    character(len=100), parameter      :: pflag = "serial"
    character(len=100), parameter      :: pflag = "cols"
!    character(len=100), parameter      :: pflag = "rows"
    ! [The Initialization flag, the options are, 
    !  'rand'  : Initializes the domain to be randomly alive or dead, 
    !  'glide' : Initializes the domain to have glider formation in top left.]
    character(len=100), parameter      :: iflag = 'glide'
    ! [The write out put flag, boolean, options are, 
    !  .TRUE.  : Life Matrix Output is written to desginated file.
    !  .FALSE. : Life Mat no written, reasons might be for timing studies.]
    logical, parameter                 :: woflag = .TRUE.
!    logical, parameter                 :: woflag = .FALSE.
    ! [The header of the save file, format will be '{savefile_head}_{k}.dat'
    !  where k is the given time step for the save.]
    character(len=8), parameter        :: savefile_head = "life_out"
    ! [The out directory for the serial implementation.]
    character(len=18), parameter       :: serial_outdir = "output/serial_out/"
    ! [The out directory for the column parallelization.]
    character(len=16), parameter       :: cols_outdir = "output/cols_out/"
    ! [The out directory for the row parallelization.]
    character(len=16), parameter       :: rows_outdir = "output/rows_out/"
    ! [The time out dir written to the global output dir since this program 
    !  can only be executed for on eparallelization at a time.]
    character(len=15), parameter         :: time_outfile = "output/time.dat"

    ! Problem Variables
    ! -----------------
    ! [The number of columns, will be read in from Nx input file above.]
    integer                            :: Nx
    ! [The number of rows, will be read in from Ny input file above.]
    integer                            :: Ny
    ! [The parent matrix, the entire domain.]
    integer, allocatable               :: A(:,:)
    ! [The Nx plus one and Ny plus one variables.]
    integer                            :: Nxp1 
    integer                            :: Nyp1
    ! [The time elapsed for given implementation.]
    real (dp)                          :: t

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
    integer                            :: Np



    ! PROGRAM START
    ! =============

    ! START PARALLEL
    ! --------------

    ! FILE INPUTS
    ! -----------
    ! [Read in the Nx.]
    OPEN(10, file=Nx_infile)
    READ(10,*) Nx
    CLOSE(10)
    ! [Read in the Ny.]
    OPEN(11, file=Ny_infile)
    READ(11,*) Ny
    CLOSE(11)

    Nxp1 = Nx + 1
    Nyp1 = Ny + 1

    ! [Init MPI.]
    CALL MPI_INIT(ierr)
    ! [ID of current processor.]
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
    ! [The number of processors.]
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Np, ierr)

    ! [This creates and allocates the parent matrix.]
    CALL Init_Life(pid, master, iflag, Nx, Ny, A)

    ! [Run the life model for given parallelization flag.]
    IF (pflag .EQ. 'serial') THEN  

        PRINT *, "Running Game of Life in series..."
        PRINT *, "---------------------------------"

        CALL Run_Serial_Life(pid,                                    &
        &                    master,                                 & 
        &                    Nx,                                     & 
        &                    Ny,                                     &
        &                    Nt,                                     &
        &                    Nw,                                     & 
        &                    Nxp1,                                   &
        &                    Nyp1,                                   & 
        &                    A,                                      &
        &                    savefile_head,                          &
        &                    serial_outdir,                          & 
        &                    woflag,                                 &
        !                    [OUTPUT]                                !
        &                    t)
 
    ! [Run the life model for column paralleization.]
    ELSE IF (pflag .EQ. 'cols') THEN 

        PRINT *, "Running Game of Life in parallel: Column Partitions"
        PRINT *, "---------------------------------------------------"
        PRINT *, "Hello from processor:", pid
        PRINT *, "---------------------        ----"
        PRINT *, ""

        CALL Run_Column_Life(pid,                                    &
        &                    master,                                 &
        &                    Np,                                     &
        &                    Nx,                                     &
        &                    Ny,                                     &
        &                    Nt,                                     &   
        &                    Nw,                                     &
        &                    A,                                      &
        &                    savefile_head,                          &
        &                    cols_outdir,                            &
        &                    woflag,                                 &
        !                    [OUTPUT]                                ! 
        &                    t)
 
    ELSE IF (pflag .EQ. 'rows') THEN 

        PRINT *, "Running Game of Life in parallel: Row Partitions"
        PRINT *, "---------------------------------------------------"
        PRINT *, "Hello from processor:", pid
        PRINT *, "---------------------        ----"
        PRINT *, ""

        CALL Run_Row_Life(pid,                                       &
        &                    master,                                 &
        &                    Np,                                     &
        &                    Nx,                                     &
        &                    Ny,                                     &
        &                    Nt,                                     &   
        &                    Nw,                                     &
        &                    A,                                      &
        &                    savefile_head,                          &
        &                    rows_outdir,                            &
        &                    woflag,                                 &
        !                    [OUTPUT]                                ! 
        &                    t)
 
!    ELSE IF (pflag .EQ. 'tile') THEN
!        CALL Run_Tile_Life()
    END IF

    ! [Write the elapsed time to file.]
    ! [The time is only on master.]
    IF (pid .EQ. master) THEN 
        ! [Print the elapsed time.]
        PRINT *, "TOTAL ELAPSED COMPUTE TIME: ", t
        ! [Write to the time_outfile.]
        OPEN(12, file=time_outfile)
        WRITE(12,*) t
        CLOSE(12)
    ENDIF

    ! [This deallocates the matrix A on master.]
    CALL End_Life(pid, master, A)

    ! [Terminate MPI.]
    CALL MPI_FINALIZE(ierr)

END PROGRAM main
