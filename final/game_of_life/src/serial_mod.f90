MODULE serial_mod

    USE mpi
    
    USE funcs_mod, ONLY: Init_Life
    
    implicit none

    ! Parameters
    ! ----------
    ! [Double.]
    integer, parameter                 :: dp = kind(0.d0)

    CONTAINS 

        !------------------------------------------------------------
        SUBROUTINE Run_Serial_Life(pid,                             &
        &                          master,                          & 
        &                          Nx,                              & 
        &                          Ny,                              &
        &                          Nt,                              &
        &                          Nw,                              & 
        &                          A,                               &
        &                          savefile_head,                   &
        &                          serial_outdir,                   & 
        !                          [OUTPUT]                         !
        &                          t)
        !------------------------------------------------------------
            ! This is th eusbroutine for running the game of life in 
            ! serial on the master processor.
            
            ! Parameters
            ! ----------
            ! [The processor id.]
            integer, intent(in)                  :: pid
            ! [The master id.]
            integer, intent(in)                  :: master
            ! [The x,y,t dimensions.]
            integer, intent(in)                  :: Nx, Ny, Nt
            ! [The write interval.]
            integer, intent(in)                  :: Nw
            ! [The parent life matrix.]
            integer, intent(in)                  :: A(Ny, Nx)
            ! [The save file header.]
            character(len=100), intent(in)       :: savefile_head
            ! [The serial output directory.]
            character(len=100), intent(in)       :: serial_outdir

            ! Returns
            ! -------
            ! [The elapsed time for the algorithm.]
            real (dp), intent(out)               :: t
            




        END SUBROUTINE Run_Serial_Life

END MODULE serial_mod

