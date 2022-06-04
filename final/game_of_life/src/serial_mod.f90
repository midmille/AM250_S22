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
    SUBROUTINE UGN(Nx,                                          &
    &              Ny,                                          &
    &              Nxp1,                                        & 
    &              Nyp1,                                        &
    &              Ag)
    !------------------------------------------------------------
        ! UGN: Update Ghost Nodes.
        ! This routine updates the ghost nodes of the provided 
        ! matrix for periodic BC. 
        ! The provided matrix should already include a halo.

        ! Parameters
        ! ----------
        ! [The Nx and Ny dimensions.]
        integer, intent(in)                  :: Nx
        integer, intent(in)                  :: Ny
        ! [The Nxp1 and Nyp1 dimensions.]
        integer, intent(in)                  :: Nxp1
        integer, intent(in)                  :: Nyp1

        ! Returns 
        ! -------
        ! [The haloed matrix whose ghosty nodes will be updated.].
        integer, intent(inout)               :: Ag(0:Nyp1, 0:Nxp1)

        ! ROUTINE START
        ! =============

        ! [Implement row exchange.]
        Ag(0,1:Nx) = Ag(Ny,1:Nx)
        Ag(Nyp1,1:Nx) = Ag(1,1:Nx)

        ! [Implement column exchange (includes diagonals here).]
        Ag(0:Nyp1,0) = Ag(0:Nyp1, Nx) 
        Ag(0:Nyp1,Nxp1) = Ag(0:Nyp1, 1)

    END SUBROUTINE UGN


    !------------------------------------------------------------
    SUBROUTINE Write_Life_Step(k,                               & 
    &                          Nx,                              & 
    &                          Ny,                              & 
    &                          A,                               & 
    &                          savefile_head,                   &
    &                          outdir)
    !------------------------------------------------------------
        ! This algorithm writes the matrix A to file for the current
        ! time step.

        ! Parameters
        ! ----------
        ! [The current time step.] 
        integer, intent(in)                  :: k 
        ! [The number of columns.]
        integer, intent(in)                  :: Nx
        ! [The number of rows.]
        integer, intent(in)                  :: Ny
        ! [The matrix at step k.]
        integer, intent(in)                  :: A(Ny, Nx)
        ! [The savefile header.] 
        integer

        ! Returns
        ! -------
        ! 

        ! Variables
        ! ---------
        ! 




    !------------------------------------------------------------
    SUBROUTINE Run_Serial_Life(pid,                             &
    &                          master,                          & 
    &                          Nx,                              & 
    &                          Ny,                              &
    &                          Nt,                              &
    &                          Nw,                              & 
    &                          Nxp1,                            &
    &                          Nyp1,                            & 
    &                          A,                               &
    &                          savefile_head,                   &
    &                          serial_outdir,                   & 
    &                          wo,                              &
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
        ! [The Nxp1 and Nyp1 variables.]
        integer, intent(in)                  :: Nxp1, Nyp1
        ! [The parent life matrix.]
        integer, intent(in)                  :: A(Ny, Nx)
        ! [The save file header.]
        character(len=*), intent(in)         :: savefile_head
        ! [The serial output directory.]
        character(len=*), intent(in)         :: serial_outdir
        ! [The write output flag.]
        logical, intent(in)                  :: wo

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


        ! ROUTINE START
        ! =============

        ! [Run this algorithm only on master.]
        IF (pid .EQ. master) THEN 

            ! [Init Ag to zeros.]
            Ag = 0
            ! [Init center of Ag to be A.]
            Ag(1:Ny, 1:Nx) = A 
    
            ! [Loop time.]
            DO k=1,Nt

                ! [Update Ghost Cells in Ag.]
                CALL UGN(Nx, Ny, Nxp1, Nyp1, Ag)

                ! [Run Life.]
                CALL Life(Nxp1, Nyp1, Ag)  

                ! [If the write output flag is true.]
                IF (wo .EQ. .TRUE.) THEN 
                    ! [Write first step and the output every Nw time step.]
                    IF (k .EQ. 1) .OR. (MOD(k,Nw) .EQ. 0) THEN 
                        CALL Write_Step()
                    ENDIF 
                ENDIF

            ENDDO

        ENDIF 




    END SUBROUTINE Run_Serial_Life

END MODULE serial_mod

