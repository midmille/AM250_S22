MODULE parallel_cols_mod
    
    USE funcs_mod, ONLY: Print_Thing1
    implicit none

CONTAINS
    SUBROUTINE Print_Thing2()
        PRINT *, 'Thing2'
        CALL Print_Thing1
    END SUBROUTINE Print_Thing2

END MODULE parallel_cols_mod
