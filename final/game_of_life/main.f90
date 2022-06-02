! This program is the main driver routine for the game of life.

PROGRAM main

    USE funcs_mod, ONLY: Print_Thing1
    USE parallel_cols_mod, ONLY: Print_Thing2

    implicit none

    CALL Print_Thing1()
    CALL Print_Thing2()

END PROGRAM main
