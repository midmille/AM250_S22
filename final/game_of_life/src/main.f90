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
    integer, parameter                 :: Nx = 10
    ! [The size of the y-dimension, the number of rows.]
    integer, parameter                 :: Ny = 10
    ! [The number of time steps to take.]
    integer, parameter                 :: Nt = 40
    ! [The number of time steps between each gather and write.]
    integer, parameter                 :: Nw = 1
    ! [The parallelization flag, the options are, 
    !  'serial' : Runs the problem in series on master--No Parallel,
    !  'cols'   : Parallelizes the problem into column partitions, 
    !  'rows'   : Parallelizes the problem into rows paritions,
    !  'tile'   : Parallelizes the problem into rectangular tile partitions.]
    character(len=100), parameter      :: pflag = "rows"
    ! [The Initialization flag, the options are, 
    !  'rand'  : Initializes the domain to be randomly alive or dead, 
    !  'glide' : Initializes the domain to have glider formation in top left.]
    character(len=100), parameter      :: iflag = 'glide'
    ! [The write out put flag, boolean, options are, 
    !  .TRUE.  : Life Matrix Output is written to desginated file.
    !  .FALSE. : Life Mat no written, reasons might be for timing studies.]
    logical, parameter                 :: woflag = .TRUE.
    ! [The header of the save file, format will be '{savefile_head}_{k}.dat'
    !  where k is the given time step for the save.]
    character(len=8), parameter        :: savefile_head = "life_out"
    ! [The out directory for the serial implementation.]
    character(len=18), parameter       :: serial_outdir = "output/serial_out/"
    ! [The out directory for the column parallelization.]
    character(len=16), parameter       :: cols_outdir = "output/cols_out/"

    ! Problem Variables
    ! -----------------
    ! [The parent matrix, the entire domain.]
    integer, allocatable               :: A(:,:)
    ! [The Nx plus one and Ny plus one variables.]
    integer                            :: Nxp1 = Nx+1
    integer                            :: Nyp1 = Ny+1
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

    ! [Init MPI.]
    CALL MPI_INIT(ierr)
    ! [ID of current processor.]
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
    ! [The number of processors.]
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Np, ierr)

    ! [This creates and allocates the parent matrix on master processor.]
    CALL Init_Life(pid, master, iflag, Nx, Ny, A)

    ! [BARRIER??????

    ! [Run the life model for given parallelization flag.]
    IF (pflag .EQ. 'serial') THEN  
        PRINT *, "Running Game of Life in series..."
        PRINT *, "---------------------------------"
        PRINT *, "Hello from processor:", pid
        PRINT *, "---------------------        ----"
        PRINT *, ''

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
        &                    cols_outdir,                            &
        &                    woflag,                                 &
        !                    [OUTPUT]                                ! 
        &                    t)
 
!    ELSE IF (pflag .EQ. 'tile') THEN
!        CALL Run_Tile_Life()
    END IF

    ! [This deallocates the matrix A on master.]
    CALL End_Life(pid, master, A)

    ! [Terminate MPI.]
    CALL MPI_FINALIZE(ierr)

END PROGRAM main
