! This file implements the solution to homework 6 question 2 for AM 250. 
! This is the matrix multiplication of two matrices. 
!
! one multiplication will be performed in parallel while the other in
! series. 

PROGRAM mat_mul

    implicit none

    ! Parameters
    ! ----------
    ! [The double precision scalar.]
    integer, parameter                 :: dp=kind(0.d0)
    ! [The size of the matrices.]
    integer, parameter                 :: N=2000

    ! Variables
    ! ---------
    ! Problem Variables:
    ! [The two matrices to be multiplied.]
    real (dp)                          :: A(N,N), B(N,N)
    ! [The product matrix, C]
    real (dp)                          :: C(N,N)
    ! [The resulting minloc array.]
    integer                            :: res(2)
    ! [The minimum value of C and the shared array version, and its thread temp.]
    real (dp)                          :: m, mt
    ! [The locations of the minimum value of C and their shared array versions.]
    integer                            :: mi, mj, mit, mjt
    ! [Indexing.]
    integer                            :: i,j,k
    ! [Timer variables.]
    real (dp)                          :: t1, t2
    ! OMP Varaibles: 
    ! [The thead number]
    integer                            :: TID, OMP_GET_THREAD_NUM
    ! [initrinsi omp wall time function.]
    real (dp)                          :: OMP_GET_WTIME

    ! [Initialize A and B to be random number arrays.]
    CALL RANDOM_NUMBER(A)
    CALL RANDOM_NUMBER(B)


    ! ----------------------------------------------------------------------
    ! [Parallel matrix multiplication and minimum value done using DO loops.]
    ! -----------------------------------------------------------------------

    ! [Start timer.]
    t1 = OMP_GET_WTIME()
! [Initialize this parallel block.]
!$OMP PARALLEL SHARED(A,B,C, m, mi, mj) PRIVATE(i, j, k, mt, mit, mjt, t1, t2)

    ! [Initialize the minimum and its temp to be the largest the minimum could
    !  possibly be, N, since the Random numbers will all be less than one.]
    m = N
    mt = N
! [Parallelize the DO loops.]
!$OMP DO
    ! [Fortran is column major so much faster to use the cached memory by 
    !  looping over the columns first.]
    ! [Loop columns.]
    DO j=1,N
        ! [Loop the rows.]
        DO i=1,N
            ! [Loop over the sumnation for the mat mul.]
            DO k=1,N
                ! [Perform the matrix multiplication.]
                C(i,j) = C(i,j) + A(i,k)*B(k,j)
            END DO
            ! [Find minimum after the sumnation for C(i,j) is complete.]
            ! [If the first row of C.]
            IF (C(i,j) .LT. mt) THEN 
                ! [Each thread should have its own version of these variables.]
                mt = C(i,j)
                mit = i
                mjt = j
            END IF 
        END DO 
    END DO
! [End the paralleziation of the DO loops.]
!$OMP END DO

! [This critical section is to determine which threads has the minimum mt value
!  and then save its location saved as well.]
!$OMP CRITICAL
    IF (mt < m) THEN 
        m = mt
        mi = mit
        mj = mjt
    END IF 
!$OMP END CRITICAL

!$OMP MASTER
    ! [Print the result from master only.]
    PRINT *, 'Parallel DO Loop Implementation:'
    PRINT *, '--------------------------------'
    PRINT *, 'MIN: ', m, 'ROW: ', mi, 'COL: ', mj
    PRINT *, ''
!$OMP END MASTER

! [Wrap up parallelization.]
!$OMP END PARALLEL
    ! [End timer.]
    t2 = OMP_GET_WTIME()
    ! [Elapsed time.]
    PRINT *, 'Time Do Loop: ', t2 - t1
    

    ! ----------------------------
    ! [Workshare Parallelization.]
    ! ----------------------------

    ! [Start timer.]
    t1 = OMP_GET_WTIME()
! [Initialize Parallelization.]
!$OMP PARALLEL SHARED(A, B, C, res)

!$OMP WORKSHARE
    ! [These twon processes are parallelized sequentially.]
    C = MATMUL(A,B)
    res = MINLOC(C)
!$OMP END WORKSHARE

!$OMP MASTER
    ! [Print the result from master.]
    PRINT *, 'Work Share Implementation:'
    PRINT *, '---------------------------'
    PRINT *, 'MIN: ', C(res(1), res(2)), 'ROW: ', res(1), 'COL: ', res(2)
    PRINT *, ''
!$OMP END MASTER
!$OMP END PARALLEL
    ! [End timer.]
    t2 = OMP_GET_WTIME()
    ! [Elapsed time.]
    PRINT *, 'Time Work Share: ', t2 - t1

END PROGRAM 
