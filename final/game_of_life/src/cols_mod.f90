MODULE parallel_cols_mod
    
    USE mpi 

    USE funcs_mod, ONLY: Print_Thing1

    implicit none

CONTAINS
    
    !-----------------------------------------------------------------
    SUBROUTINE Run_Column_Life(pid,                                  & 
    &                          master,                               &
    &                          Nx,                                   &
    &                          Ny,                                   &
    &                          Nt,                                   & 
    &                          Nw,                                   &
    &                          A,                                    &
    &                          savefile_head,                        &
    &                          cols_outdir,                          &
    &                          woflag,                               &
    !                          [OUTPUT]                              ! 
    &                          t)
    !-----------------------------------------------------------------
        ! This routine solves the game of life by partitioning the
        ! parallelization over the columns of the domain.

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
        character(len=*), intent(in)         :: savefile_head
        ! [The serial output directory.]
        character(len=*), intent(in)         :: serial_outdir
        ! [The write output flag.]
        logical, intent(in)                  :: woflag

        ! Returns
        ! -------
        ! [The elapsed time for the algorithm.]
        real (dp), intent(out)               :: t

        ! Variables
        ! ---------
        ! [For the serial application the ghost cells are implemented from
        !  begining. Thus this A with halo.]
        integer                              :: Ag(0:Nyp1, 0:Nxp1)
        ! [Time indexing.]
        integer                              :: k

    END SUBROUTINE Run_Column_Life

 
END MODULE parallel_cols_mod
