MODULE LinAl
  implicit none

  integer, parameter       :: dp=kind(0.d0)
  !integer, save :: msize, nsize
  !real, dimension(:,:), allocatable, save :: mat

CONTAINS

  !********************************************************

  !----------------------------------------------
  SUBROUTINE read_mat(filename,                 &
  &                   msize,                    &
  &                   nsize,                    &
  &                   dims_in_file,             & 
  &                   matrix)
  !----------------------------------------------
    ! This file reads in a matrix from file. 
    ! 
    ! Parameters
    ! ----------
    ! filename: String
    !     The name of the file. 
    ! msize: Integer
    !     The number of rows in the file. 
    ! nsize: Integer
    !     The number of columns in the file.
    ! dims_in_file: Boolean
    !     This flag is used to determine if the first line of the file contains
    !     the dimensions of the matrix. If this flag is set to FALSE then the 
    !     routine will assume the matrix starts on th efirst line of the column. 
    !     If this flag is set to TRUE then the routine assumes that the
    !     dimensions are written on the first line of the file and that the
    !     matrix starts on the second line. 
    ! matrix: Float, (msize, nsize)
    !     The input zero array to be filled with data from file. 
    ! 
    ! Returns
    ! -------
    ! matrix: Float, (msize, nsize)
    !     The array that is read from file. 
    implicit none
    
    ! The name of the file, a str.
    character(len=*)          :: filename
    ! The number of rows.  
    integer, intent(in)       :: msize 
    ! The number of columns.
    integer, intent(in)       :: nsize 
    ! dims_in_file flag. 
    logical, intent(in)       :: dims_in_file
    ! The matrix to read values for.
    real (dp), intent(inout)  :: matrix(msize, nsize)


    integer :: i,j

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
    ! Note that entries must be separated by a tab.


    open(10,file=filename)

    ! If the dims of the matrix are included in file. 
    IF (dims_in_file .EQV. .TRUE.) THEN
      ! Read the matrix dimensions
      READ(10,*) i,j
      ! Read matrix
      DO i=1,msize
         READ(10,*) ( matrix(i,j), j=1,nsize )
      ENDDO

    ELSEIF (dims_in_file .EQV. .FALSE.) THEN
      ! Read thre matrix starting at first line.
      DO i=1,msize
        READ(10,*) (matrix(i,j), j=1,nsize)
      ENDDO

   ENDIF
   
   CLOSE(10) 
     
  END SUBROUTINE read_mat


  !--------------------------------------------------------
  !-----------------Write-Matrix-To-File-------------------
  SUBROUTINE write_mat(filename,                          &
  &                    ma,                                &
  &                    na,                                &
  &                    dims_in_file,                      &
  &                    A)
  !--------------------------------------------------------
    ! This file writes a given matrix into a file in acsii form. 
    !
    ! Parameters
    ! ----------
    ! filename: String 
    !     The name of file. 
    ! ma: Integer
    !     The number of rows in matrix.
    ! na: Integer
    !     The number of columns in the matrix. 
    ! dims_in_file: Boolean
    !     This flag is used if it is desired for the dimensions of 
    !     the matrix A are to be written to first line of the file. 
    !     if the flag is .FALSE. then the dimensions are not written and 
    !     the first row of the matrix is written on the first line. If the 
    !     flag is set to .TRUE. then the dimensions of the given matrix are 
    !     written to the first line of the file.
    ! A: Float, (ma,na)
    !     The matrix who will be written to file. 
    ! 
    ! Returns
    ! -------
    ! None, [Just creates the file]

    ! Global Vars
    ! -----------
    ! The name of the file. 
    character(len=*)                   :: filename
    ! The number of rows and cols. 
    integer, intent(in)                :: ma, na
    ! The write dimensions flag. 
    logical, intent(in)                :: dims_in_file
    ! Th ematrix to write to file. 
    real (dp), intent(in)              :: A(ma,na)
    
    ! Local Vars
    ! ----------
    ! Indexing 
    integer                            :: i,j

    OPEN(10, file=filename)

    ! If the dims are to be written. 
    IF (dims_in_file .EQV. .TRUE.) THEN 
      ! Write the dimensions. 
      WRITE(10,*) ma,na
      ! Write the matrix... 
      DO i=1,ma
        WRITE(10,*) (A(i,j), j=1,na)
      ENDDO

   ! If the dimensions are not to be written. 
   ELSEIF (dims_in_file .EQV. .FALSE.) THEN 
     ! Write matrix starting at first line.
      DO i=1,ma
        WRITE(10,*) (A(i,j), j=1,na)
      ENDDO

    ENDIF

    CLOSE(10)

  END SUBROUTINE write_mat


  !-----------------------------------------------
  SUBROUTINE calc_trace(A,                       &
  &                     m,                       &
  &                     trace)
  !-----------------------------------------------
    ! This function calculates the trace of the square matrix A 

    ! Global Vars
    ! -----------
    ! The square matrix A
    real (dp), intent(in)    :: A(:,:)
    ! The dimension of A 
    integer, intent(in)      :: m
    ! The output trace
    real (dp), intent(out)    :: trace
  
    ! Local Vars
    ! ----------
    ! This is just a indexing variable.
    integer                  :: k
  
  
    ! All that is necessary is to loop over the m dimension 
    ! and add all the a_ii diagonal values. 
    ! Initializing tr
    trace = 0.0_dp
    ! Looping over the diagonal of A
    DO k=1,m
      trace = trace + A(k,k)
    END DO 
  
  END SUBROUTINE calc_trace


  !----------------------------------------------
  SUBROUTINE calc_two_norm(v,                   &
  &                        m,                   &
  &                        two_norm) 
  !----------------------------------------------
    ! This function calculates the Euclidean norm of 
    ! vector v. 

    ! Global Vars
    ! -----------
    ! The vector whose norm we are taking.
    real (dp), intent(in)    :: v(:)
    ! The dimension of the vector.
    integer, intent(in)      :: m
    ! The Euclidian norm output.
    real (dp), intent(out)   :: two_norm

    ! Local Vars 
    ! ----------
    ! Indexing variable
    integer                  :: k
    ! The square of the two norm
    real (dp)                :: two_norm_sqrd

    ! Init the two norm 
    two_norm_sqrd = 0.0_dp 
    DO k=1,m
      two_norm_sqrd = two_norm_sqrd + v(k) ** 2
    END DO 

    ! Taking the sqrt to get the 2-norm.
    two_norm = SQRT(two_norm_sqrd) 

  END SUBROUTINE calc_two_norm

  !--------------------------------------------------------
  !-----------------Calculates-Frobenius-Norm--------------
  SUBROUTINE calc_f_norm(ma,                              &
  &                      na,                              &
  &                      A,                               &
  &                      f_norm) 
  !--------------------------------------------------------
    ! This function calculates the frobenius norm of a matrix A
    
    ! Global Vars
    ! -----------
    ! The dimensions of A.
    integer, intent(in)                :: ma, na
    ! The input matrix A.
    real (dp), intent(in)              :: A(ma,na)
    ! The output frobenius norm.
    real (dp), intent(out)             :: f_norm  

    ! Local Vars 
    ! ----------
    ! Indexing 
    integer                            :: i
    ! The trace of AtA
    real (dp)                          :: trace

    ! The start of the routine.
    ! first to calculate the trace of AtA. 
    CALL calc_trace(MATMUL(A, TRANSPOSE(A)), ma, trace)
    ! The frobennius is the sqrt of the trace. 
    f_norm = SQRT(trace)

  END SUBROUTINE calc_f_norm

  !-----------------------------------------------
  SUBROUTINE print_matrix(A,                     &
  &                       m,                     &
  &                       n)                
  !-----------------------------------------------
    ! This function prints the matrix A and its given dimensions

    ! Global Vars
    ! -----------
    ! The matrix to be printed
    real (dp), intent(in)    :: A(:,:)
    ! The number of rows 
    integer, intent(in)      :: m
    ! The  number of columns
    integer, intent(in)      :: n 

    ! Local Vars
    ! ----------
    ! Some Loop Vars
    integer                  :: i, j
    
    print *, 'n_rows:', m
    print *, 'n_cols:', n 
    do i = 1, m
      write(*,*) (A(i,j) , j = 1, n )
    end do


  END SUBROUTINE print_matrix


  !----------------------------------------------
  !--Gaussian Elimination with Partial Pivoting--
  SUBROUTINE gauss_elim_pp(ma,                  &
  &                        mb, nb,              &
  &                        A,                   &
  &                        B,                   &
  &                        singular)
  !----------------------------------------------
    ! This function performs gaussian elimination with 
    ! partial pivoting to calculate the solution to Ax = B.

    ! Global Vars
    ! -----------
    ! ma is the dimension of the square coefficient matrix A.
    integer, intent(in)        :: ma
    ! mb is the number of rows and nb in the number of columns of 
    ! the rhs matrix B.
    integer, intent(in)        :: mb, nb
    ! A is the square coefficient matrix. It will be rewritten as
    ! an upper triangular matrix after gaussian elimination is applied
    real  (dp), intent(inout)  :: A(ma, ma)
    ! B is the rhs matrix. It will be returned as the correspnding and
    ! modified rhs vectors. 
    real (dp), intent(inout)   :: B(mb, nb)
    ! This is the singular or nonsingular flag. 
    ! It returns .TRUE. if singular and
    ! .FALSE. is non-singular. 
    logical, intent(out)       :: singular

    ! Local Vars
    ! ----------
    ! Some indexes
    integer                    :: i,j
    ! Another index for the pivot location. 
    integer                    :: k
    ! The value of a_jj and a_ij
    real (dp)                  :: a_jj, a_ij
    ! Place holder for row to be swapped in A pivot.
    real (dp)                  :: ra_j(ma), ra_k(ma)
    ! The same place holders but for the rows of matrix b
    real (dp)                  :: rb_j(nb), rb_k(nb)

    ! Checking that ma = mb
    IF (ma .NE. mb) THEN 
      STOP('mu != mb error')
    ENDIF
    ! The start of the Gaussian Elimination Algorithm
    ! Loop over the columns
    DO j = 1, ma-1
      ! Finding the k index, index of the row that needs to be swaped
      ! with. Using MAXLOC function to find that index. 
      !  Must add j to index because maxloc finds index at 1.
      k = j + MAXLOC(ABS(A(j:,j)), 1) - 1 
      ! Then checking if k!=j
      IF (k .NE. j) THEN
        ! Saving the rows. 
        ra_j = A(j,:)
        ra_k = A(k,:) 
        rb_j = B(j,:)
        rb_k = B(k,:)
        ! Interchanging the rows.
        A(j,:) = ra_k
        A(k,:) = ra_j
        B(j,:) = rb_k
        B(k,:) = rb_j
      ENDIF 
      ! Check for singular or not. 
      IF (A(j,j) .EQ. 0) THEN
        singular = .TRUE.
      ELSE 
        singular = .FALSE.
      ENDIF 
      ! Loop over the rows under row j. 
      DO i=j+1,ma 
        ! Aij and Ajj
        a_jj = A(j,j)
        a_ij = A(i,j)
        ! Elimination of the leading coefficient of row i. 
        A(i,:) = A(i,:) - (a_ij * A(j,:)) / a_jj
        ! Doing same transformation on the columns of B. 
        B(i,:) = B(i,:) - (a_ij * B(j,:))/ a_jj
      ENDDO
    ENDDO
  END SUBROUTINE gauss_elim_pp 


  !----------------------------------------------
  !-----------Backsubstitution-------------------
  SUBROUTINE gauss_back_substitution(mu,        &
  &                                  mb, nb,    &                
  &                                  U,         &
  &                                  B,         &
  &                                  X) 
  !----------------------------------------------
    ! This function performs backsubstitution using 
    ! the upper triangular matrix U and the rhs matrix
    ! with vector solutions B. 

    ! Global Vars
    ! -----------
    ! This is the dimension of the square matrix u. 
    integer, intent(in)         :: mu
    ! The dimensions of B, where mb=nrows, nb=ncols
    integer, intent(in)         :: mb, nb
    ! The upper triangular matrix U
    real (dp), intent(in)       :: U(mu, mu) 
    ! The rhs matrix B
    real (dp), intent(in)       :: B(mb, nb)
    ! The return of the solution matrix X, should be same 
    ! dimensions as B.
    real (dp), intent(inout)      :: X(mb, nb) 

    ! Local Vars
    ! ----------
    ! Indexing vars, 
    integer                       :: i,k
    ! The sum variable. 
    real (dp)                     :: tot(nb) 
    
    tot = 0.0_dp
    ! Checking that mu = mb
    IF (mu .NE. mb) THEN 
      STOP('mu != mb error')
    ENDIF

    ! Start of backsubstitution
    ! The bottom most row of x is just y/u s.t.
    X(mb,:) = B(mb,:) / U(mb,mb)

    ! Loop over the rows of U. 
    DO i=mu-1, 1, -1 
      ! Check for singularity
      if (U(i,i) .EQ. 0) THEN
        STOP 
      ENDIF
      ! Init the sum
      tot = 0.0_dp
      ! Loop over the columns greater than the
      ! index of the current row.
      DO k=i+1,mu 
        tot = tot + U(i,k)*X(k,:)
      ENDDO 
      ! Solving for X
      X(i,:) = (B(i,:) - tot) / U(i,i)
    ENDDO

  END SUBROUTINE gauss_back_substitution    


  !----------------------------------------------
  !---------LU-Decomposition-Partial-Pivot-------
  SUBROUTINE lu_decomp_pp(ma,                   &
  &                       A,                    &
  &                       singular,             &
  &                       s)
  !----------------------------------------------
    ! This routine solves for the LU decomposition 
    ! of the square matrix A, representing a system
    ! of linear equations. 
    
    ! Global Vars
    ! -----------
    ! This is the dimension of square matrix A. 
    integer, intent(in)      :: ma
    ! This is the square matriz A to be decompoesed
    real (dp), intent(inout) :: A(ma,ma)
    ! This is a logical flag s.t. True=singular, false=nonsingular. 
    logical, intent(inout)   :: singular
    ! This is the permutation vector. 
    integer, intent(inout)   :: s(ma)
    
    ! Local Vars 
    ! ----------
    ! Indexes
    integer                  :: i,j,k
    ! The rows for interchaning in pivot. 
    real (dp)                :: r_j(ma), r_k(ma)
    ! The scalars for recording permutation.
    integer                  :: s_j, s_k

    ! The start of the decomposition algorithm. 
    ! First inititalize the permutation vector. 
    DO j=1,ma 
      s(j) = j
    ENDDO

    ! Loop over the columns of A. 
    DO j=1,ma
      ! Finding the k index, index of the row that needs to be swaped
      ! with. Using MAXLOC function to find that index. 
      !  Must add j to index because maxloc finds index at 1.
      k = j + MAXLOC(ABS(A(j:,j)), 1) - 1 
      ! Pivot if k != j. 
      IF (k .NE. j) THEN 
        ! Saving rows and permutation vals. 
        r_j = A(j,:)
        r_k = A(k,:)
        s_j = s(j)
        s_k = s(k)
        ! Swapping the k and j vals. 
        A(j,:) = r_k
        A(k,:) = r_j
        s(j) = s_k
        s(k) = s_j
      ENDIF 
      ! Checking if singular. 
      IF (A(j,j) .EQ. 0) THEN
        singular = .TRUE. 
      ELSE 
        singular = .FALSE. 
      ENDIF
      ! Loop over the rows of A. 
      DO i = j+1, ma
        !  Here the l_ij is stored back in A. 
        A(i,j) = A(i,j) / A(j,j)
        DO k=j+1,ma 
          A(i,k) = A(i,k) - A(i,j) * A(j,k)
        ENDDO 
      ENDDO
    ENDDO

  END SUBROUTINE lu_decomp_pp


  !--------------------------------------------------------
  !------------------LU-Backsubstitution-------------------
  SUBROUTINE lu_back_substitution(mlu,                    &
  &                               mb, nb,                 &
  &                               LU,                     & 
  &                               B,                      &
  &                               s,                      &
  &                               X)                   
  !--------------------------------------------------------
    ! This function performs lu-backsubstitution on the LU 
    ! matrix to solve for the solution matrix X. The rhs matrix 
    ! B is not necessarily square. 

    ! Global Vars
    ! -----------
    ! This is the dimension of the square LU matrix.
    integer, intent(in)                :: mlu
    ! The dimensions of B.
    integer, intent(in)                :: mb, nb
    ! The LU matrix. 
    real (dp), intent(in)              :: LU(mlu, mlu)
    ! The rhs matrix.
    real (dp), intent(in)              :: B(mb,nb)
    ! The permutation vector. 
    integer, intent(in)                :: s(mlu)
    ! The solution matrix
    real (dp), intent(inout)           :: X(mb,nb)

    ! Local Vars 
    ! ----------
    ! Indexing. 
    integer                            :: i,j,k
    ! The Lower matrix solution vectors. 
    real (dp)                          :: Y(mb,nb)
    !  The sum variable. 
    real (dp)                          :: tot(nb)
    

    Y = 0.0_dp
    ! Start of lu backsubstitution. 
    ! Loop over the rows of Y, B
    DO j=1,mlu 
      ! Init Y with Pb
      Y(j, :) = B(s(j),:)
    ENDDO 
    
    ! The forward substitution of X = L^-1PB. 
    ! Loop over the columns of LU. 
    DO j=1,mlu-1 
      ! Looping over the rows of a. 
      DO i=j+1,mlu
        Y(i,:) = Y(i,:) - Y(j,:) * LU(i,j)
      ENDDO 
    ENDDO

    ! The back subsitution, Ux = Y
    ! Starting at the bottom of the upper triangular.
    DO i=mlu,1,-1
      ! Check for singular. 
      IF (LU(i,i) .EQ. 0) THEN 
        STOP('Singular U')
      ENDIF
      !Init the tot as 0
      tot = 0.0_dp
      ! Loop columns. 
      DO k=i+1,mlu
        tot = tot +  LU(i,k) * X(k,:)
      ENDDO 
      ! Solving for x. 
      X(i,:) = (Y(i,:) - tot) / LU(i,i)
    
    ENDDO

  END SUBROUTINE lu_back_substitution


  !--------------------------------------------------------
  !------------------LU-Backsubstitution-Vector-Sol--------
  SUBROUTINE lu_back_substitution_vectorsol(mlu,          &
  &                                         LU,           & 
  &                                         b,            &
  &                                         s,            &
  &                                         x)                   
  !--------------------------------------------------------
    ! This function performs lu-backsubstitution on the LU 
    ! matrix to solve for the solution vector x. This routine 
    ! is essentially the same as the lu_back_substitution above 
    ! but take a vector as the rhs argument and solves for a 
    ! vector x instead of matrices.

    ! Global Vars
    ! -----------
    ! This is the dimension of the square LU matrix.
    integer, intent(in)                :: mlu
    ! The LU matrix. 
    real (dp), intent(in)              :: LU(mlu, mlu)
    ! The rhs vector.
    real (dp), intent(in)              :: b(mlu)
    ! The permutation vector. 
    integer, intent(in)                :: s(mlu)
    ! The solution vector.
    real (dp), intent(inout)           :: x(mlu)

    ! Local Vars 
    ! ----------
    ! Indexing. 
    integer                            :: i,j,k
    ! The Lower matrix solution vector. 
    real (dp)                          :: y(mlu)
    !  The sum variable. 
    real (dp)                          :: tot
    

    y = 0.0_dp
    ! Start of lu backsubstitution. 
    ! Loop over the rows of Y, B
    DO j=1,mlu 
      ! Init Y with Pb
      y(j) = b(s(j))
    ENDDO 
    
    ! The forward substitution of X = L^-1PB. 
    ! Loop over the columns of LU. 
    DO j=1,mlu-1 
      ! Looping over the rows of a. 
      DO i=j+1,mlu
        y(i) = y(i) - y(j) * LU(i,j)
      ENDDO 
    ENDDO

    ! The back subsitution, Ux = Y
    ! Starting at the bottom of the upper triangular.
    DO i=mlu,1,-1
      ! Check for singular. 
      IF (LU(i,i) .EQ. 0) THEN 
        STOP('Singular U')
      ENDIF
      !Init the tot as 0
      tot = 0.0_dp
      ! Loop columns. 
      DO k=i+1,mlu
        tot = tot +  LU(i,k) * x(k)
      ENDDO 
      ! Solving for x. 
      x(i) = (y(i) - tot) / LU(i,i)
    
    ENDDO

  END SUBROUTINE lu_back_substitution_vectorsol


  !---------------------------------------------------------
  !---------------------Cholesky-Decomposition--------------
  SUBROUTINE cholesky_decomp(ma,                           &
  &                          A,                            &
  &                          singular_or_spd)
  !---------------------------------------------------------
    ! This function performs cholesky decomposition on the 
    ! square matrix A.
    !
    ! Parameters
    ! ----------
    ! ma: Integer
    !     The dimension of the square matrix A. 
    ! A: Float, (ma, ma)
    !     The input matrix argument to be decomposed. 
    !
    ! Returns
    ! -------
    ! A: Float, (ma, ma)
    !     The output lower triangular matrix is stored in the 
    !     the lower triangle of A. 
    ! singular_or_spd: Boolean
    !     This flag is returned as True if the input matrix A
    !     is either singular or not symetric, positive definite.  

    ! Global Vars
    ! -----------
    ! The dimension of A. 
    integer, intent(in)                :: ma
    ! The matrix A. 
    real (dp), intent(inout)           :: A(ma,ma)
    ! The singular or SPD flag. 
    logical, intent(out)               :: singular_or_spd

    ! Local Vars
    ! ----------
    ! Indexing
    integer                            :: i,j,k


    ! The start of the algorithm. 
    ! Loop over the columns.  
    DO j=1,ma
      ! calculate the diagonals larger than a_11. 
      DO k=1, j-1
        A(j,j) = A(j,j) - A(j,k) * A(j,k) 
      ENDDO 
      ! sqrt of the diagonal. 
      A(j,j) = SQRT(A(j,j))

      ! The elements below the diagonal in the current column.
      DO i=j+1, ma
        ! Looping over the sum of the entries in the rows to the right
        ! of the current row. 
        DO k=1,j-1
          ! The solution for the current column, row l, stored back into A. 
          A(i,j) = A(i,j) - A(i,k)*A(j,k)
        ENDDO 
        ! Dividing by the prviously calculated diagonal entry. 
        A(i,j) = A(i,j) / A(j,j) 
      ENDDO
    ENDDO 

  END SUBROUTINE cholesky_decomp


  !--------------------------------------------------------
  !--------------------Cholesky-Backsubstitution-----------
  SUBROUTINE cholesky_back_substitution(ma,               &
  &                                     A,                &
  &                                     b,                &
  &                                     x)
  !--------------------------------------------------------
    ! This function implements the Cholesky Backsubstitution 
    ! scheme as described in the AM213A lecture pdf notes. 
    ! The routine first solves the forward substitution equation
    ! Ly = b for y. Then it solves the back substitution equation
    ! L*x = y for x. Thus giving the solution to the system, x. 
    ! 
    ! The routine assumes that the vector b is actually a matrix, 
    ! b such that each column is a seperate rhs solution to the problem. 
    ! 
    ! Paramters
    ! ---------
    ! ma: Integer 
    !     The dimension of the square matrix A. 
    ! A: Float, (ma, ma)
    !     The cholesky decomposed A matrix. The elements of the lower triangular
    !     matrix L are should be overwritten the lower-triang elements
    !     of the matrix A. 
    ! b: Float, (ma)
    !     The rhs vector b. 
    ! 
    ! Returns
    ! -------
    ! x: Float, (ma)
    !     The solution to the system. Solved via cholesky back substitution. 

    ! Global Vars
    ! -----------
    ! The dimension of A. 
    integer, intent(in)                :: ma
    ! The matrix A. 
    real (dp), intent(in)              :: A(ma,ma)
    ! The matrix B. 
    real (dp), intent(in)              :: b(ma)
    ! The solution matrix X, same size as B. 
    real (dp), intent(inout)           :: x(ma)

    ! Local Vars
    ! ----------
    ! Indexing
    integer                            :: i,j,k
    ! The sum vector, size of the columns of B. 
    real (dp)                          :: tot
    ! The Y matrix for the solution to forward substitution. 
    real (dp)                          :: y(ma)

    ! The beginnning of the routine. 
    ! The forward substitution piece, LY=B. 
    ! Loop over the columns in lower traingle of A. 
    DO i=1,ma
      ! Init the sum for the first row. 
      tot = b(i)
      ! Loop over the columns of L, forward substituting previous sols. 
      DO j=1,i-1
        tot = tot - y(j)*A(i,j)
      ENDDO
      ! Solve for y_i
      y(i) = tot / A(i,i)
    ENDDO

    ! The back substitution piece. 
    ! Looping over columns from the last to first.
    DO i=ma,1, -1
      ! Checking for singular value. 
      IF (A(i,i) .EQ. 0) THEN 
        STOP("Singular Value in Cholesky Backsubstitution")
      ENDIF
      ! Loop over the rows, skip the first i=m index. 
      DO k=i+1,ma
        ! Substituting in the previous solutions.
        y(i) = y(i) - A(k,i)*x(k)
      ENDDO 
      ! Solving for x_i
      x(i) = y(i) / A(i,i) 
    ENDDO

  END SUBROUTINE cholesky_back_substitution


  !--------------------------------------------------------
  !-------------------Householder-QR-Decomposition---------
  SUBROUTINE householder_qr_decomp(ma,                    &
  &                                na,                    &
  &                                A,                     &
  &                                Q)
  !--------------------------------------------------------
    ! This function performs a Householder based QR decomposition 
    ! on the input tall matrix A. 
    ! 
    ! Paramters
    ! ---------
    ! ma: Integer 
    !     The number of rows of A. 
    ! na: Integer
    !     The number of columns of A. 
    ! Id: Float, (ma, ma)
    !     An identity matrix.
    ! A: Float, (ma, na)
    !     The input matrix to be decomposed. Note that ma>na.
    ! Q: Float, (ma,ma)
    !     The Q matrix, input as zeros.
    ! 
    ! Returns
    ! -------
    ! A: Float, (ma, na)
    !     The output of A contains R stored as an upper triangle. 
    ! Q: Float, (ma, na)
    !     The orthogonal matrix of the QR decomposition. 
    
    ! Global Vars
    ! -----------
    ! The dimensions of A. 
    integer, intent(in)                :: ma, na
    ! The matrix A, its upper triangle will become R.  
    real (dp), intent(inout)           :: A(ma,na)
    ! The orthogonal matrix Q 
    real (dp), intent(inout)           :: Q(ma,ma)

    ! Local Vars
    ! ----------
    ! Indexing
    integer                            :: i,j,k
    ! The orthogonal vectors corresponding to each different householder.  
    real (dp)                          :: V(ma,na) 
    ! An idenity matrix with rows and columns equal to the rows of A.  
    real (dp)                          :: Id(ma,ma)
    ! A sum place holder
    real (dp)                          :: tot
    ! The signed norm 
    real (dp)                          :: s

    ! Ensure that V is filled with zeros. 
    V = 0.0_dp
    Id = 0.0_dp
    
    ! Construct the identity matrix for problem. 
    DO i=1,ma
      Id(i,i) = 1.0_dp
    ENDDO
    ! The start of the algorithim. 
    ! Loop over the columns. 
    DO j=1,na
      ! First step is to calculate the signed norm.
      tot = 0 
      ! Loop over the rows below the diagonal of current column and add
      ! to calculate the norm. 
      DO k=j,ma
        tot = tot + A(k,j) ** 2
      ENDDO 
      ! The signed norm
      ! SIGN(A,B) function returns value A with sign of B.
      s = SIGN(SQRT(tot), A(j,j)) 

      ! Now to compute the vector V
      V(j:,j) = A(j:,j)
      ! Adding the signed norm to the diagonal of the vector v. 
      V(j,j) = V(j,j) + s
      ! Compute the norm of the newly calculated v. 
      tot = 0 
      DO k=j,ma
        tot = tot + V(k,j)**2
      ENDDO 
      tot = SQRT(tot)
      ! Computing the normalized vector householder vector. 
      V(:,j) = V(:,j) / tot

      ! Calculate the updated A, R reflected by the householder matrix. 
      ! The reshape comes from the fortran completion of the outer product. 
      ! I couldn't find an intrinsic function for the outer product of two
      ! vectors in fortran bu the following reference helped with the below
      ! formulation: 
      ! http://www.new-npac.org/projects/cdroms/cewes-1998-05/reports/hpfforspmd/html/node30.html
      A = A - 2*(MATMUL(MATMUL(RESHAPE(V(:,j), (/ma,1/)), RESHAPE(V(:,j), (/1, ma/))), A))

    ENDDO

    ! A is now officially R and we need to calculate Q. 
    ! Now to loop over the columns of V and create Q. 
    ! Init Q as the first householder.
    Q = (Id - 2* MATMUL(RESHAPE(V(:,1), (/ma,1/)),RESHAPE(V(:,1), (/1,ma/))))
    DO j=2,na 
      ! Same outer product formulation as above.
      Q = MATMUL(Q,(Id - 2* MATMUL(RESHAPE(V(:,j), (/ma,1/)),RESHAPE(V(:,j), (/1,ma/)))))
    ENDDO
    ! After the algorithim A will be decomposed into Q and R, R will be returned
    ! in A and Q as its own matrix.

  END SUBROUTINE householder_qr_decomp


  !------------------------------------------------------------------
  !----------------------Full-QR-Backsubstitution--------------------
  SUBROUTINE qr_backsubstitution(mr,                                &
  &                              nr,                                &
  &                              Q,                                 &
  &                              R,                                 &
  &                              b,                                 &
  &                              x)
  !------------------------------------------------------------------
    ! This function solves the equation Ax = b, with A being fully decomposed
    ! into Q and R. This gives the true equation to be solved which is
    ! Rx = Qtb
    ! Note: I will be solving the reduced system even though the Q and R are
    ! full, this mean I will only nead nr rows of Q and R. 
    ! 
    ! Parameters
    ! ----------
    ! mr: Integer
    !     The number of rows of R, this is also the numebr of rows of A. 
    ! nr: Integer
    !     The num=ber of  columns of R, also number of cols in A. 
    ! Q: Float, (mr, mr)
    !     The Q orthogonal matrix. 
    ! R: Float,  (mr,nr)
    !     The R matrix, uppeer triangle, this is why we can use back sub to
    !     solve. 
    ! b: Float, (mr)
    !     The rhs of the original Ax=b. Dimension should be num rows of A, R. 
    ! x: Float, (nr)
    !     The input zero matrix. This is what we are solving for. 
    ! 
    ! Returns
    ! -------
    ! x: Float, (nr)
    !     The solution to Rx = Qtb

    ! Global Vars
    ! -----------
    ! dims of R. 
    integer, intent(in)                :: mr, nr
    ! Orthogonal matrix Q. 
    real(dp), intent(in)               :: Q(mr,mr)
    ! Matrix R, upper triangle. 
    real (dp), intent(in)              :: R(mr,nr)
    ! The rhs, b. 
    real (dp), intent(in)              :: b(mr)
    ! What we are solving for
    real (dp), intent(inout)           :: x(nr)

    ! Local Vars
    ! ----------
    ! Indexing: 
    integer                            :: i,j,k
    ! sum var. 
    real (dp)                          :: tot
    ! transpose(Q) = Qt, this it Qtb. 
    real (dp)                          :: Qtb(mr) 

    ! The start of the routine. 
    ! Calculating Qtb.
    Qtb = MATMUL(TRANSPOSE(Q), b)

    ! Now for the back substitution to solve for x. 
    ! Again, note I will only do this for nr rows, ie the reduced 
    ! form the Rx = Qtb. 
    
    ! Loop over the columns, starting at right side, ie 'back'.
    DO i=nr,1,-1
      ! Check for singularity. 
      IF (R(i,i) .EQ. 0.0_dp) THEN
        STOP('Singularity in QR Backsubstitution.')
      ENDIF 
      ! Init the tot. 
      tot = 0.0_dp
      ! Loop over the rows, skip i=nr index. 
      ! Note: only up to row rn, the reduced system.
      DO k=i+1,nr
        tot = tot + R(i,k) * x(k)
      ENDDO
      ! Solve for x. 
      x(i) = (Qtb(i) - tot) / R(i,i) 
    ENDDO

  END SUBROUTINE qr_backsubstitution


  !--------------------------------------------------------
  !------------------Eigenvector-Inverse-Iteration---------
  SUBROUTINE inverse_iteration(ma,                        &
  &                            A,                         &
  &                            e_val,                     & 
  &                            err,                       &
  &                            e_vec)
  !--------------------------------------------------------
    ! This function performs an inverse iteration scheme to 
    ! solve for the the eigen vector of A given an initial guess 
    ! for the corresponding eigenvalue.
    ! 
    ! Parameters
    ! ----------
    ! ma: Integer
    !     The number of rows and columns of square matrix A. 
    ! A: Float, (ma,ma)
    !     The matrix whose eigenvecotr we are solving for.
    ! e_val: Float
    !     The initial eigenvalue whose eigenvector we are solving for.
    !     The closer this is to the true value the faster the method
    !     converges to a solution.
    ! err: Float
    !     The convergence of the method. This is such that ||r||>err
    !     where r is the difference between the iteration of eigenvalues. 
    ! e_vec: Float, (ma)
    !     The eigenvecotor we are solving for. The argument should be given 
    !     as an empty array.
    ! 
    ! Returns
    ! -------
    ! e_vec: Float, (ma)
    !     The solution of the eigen vecotr for the given eigen value within
    !     the desired tolerance.

    ! Global Vars 
    ! -----------
    ! Dimension of A. 
    integer, intent(in)                :: ma 
    ! Matrix A. 
    real (dp), intent(in)              :: A(ma,ma)
    ! The eigen value.
    real (dp), intent(in)              :: e_val
    ! The tolerance for convergence. 
    real (dp), intent(in)              :: err 
    ! The eigenvector. 
    real (dp), intent(inout)           :: e_vec(ma)

    ! Local Vars
    ! ----------
    ! Indexing. 
    integer                            :: i
    ! The vector of the difference between eigenvecotrs of subsequent
    ! iterations.
    real (dp)                          :: r(ma)
    ! The two norm of r. 
    real (dp)                          :: r_norm
    ! The two norm of the eigenvector.
    real (dp)                          :: e_vec_norm
    ! The flag for lu_decomp for singularity. 
    logical                            :: singular
    ! The permutation vector of the LU decomp.  
    integer                            :: s(ma)
    ! The LU decomp filled A matrix.
    real (dp)                          :: LU(ma,ma)
    ! The identity matrix, same dims as A.
    real (dp)                          :: Id(ma,ma)
    ! This is a variable for the storage of the eigenvector before its
    ! normalized. 
    real (dp)                          :: w(ma)

    ! The start of the inverse iteration algorithm.
    
    ! inititalize any local vectors and arrays as zeros. 
    r = 0.0_dp
    s = 0.0_dp
    LU = 0.0_dp 
    Id = 0.0_dp
    w = 0.0_dp

    ! Initialize the norm of r to be greater than the designated err. 
    r_norm = 1.0_dp + err

    ! Inititalize and normalize the eigen vector. 
    e_vec = 1.0_dp
    CALL calc_two_norm(e_vec, ma, e_vec_norm) 
    e_vec = e_vec / e_vec_norm

    ! Initialize the identity matrix. 
    DO i=1,ma
      Id(i,i) = 1.0_dp
    ENDDO

    ! Init LU
    LU = A
    ! This LU decomposition will not change during the iteration step thus I
    ! only need to compute the LU once and perform the backsubstitution each
    ! iteration step.
    ! A becomes A - mu*I
    LU = (LU - e_val* Id)
    ! This returns A as an LU decomposition.
    CALL lu_decomp_pp(ma, LU, singular, s)  


    ! Iterate while the norm of the difference between subsequent eigenvalues is
    ! greater than then the error. 
    DO WHILE (r_norm > err)
      ! This returns the next iteration of the eigen vector. 
      CALL lu_back_substitution_vectorsol(ma, LU, e_vec, s, w)
      ! The new eigen vector is the normalized w. 
      CALL calc_two_norm(w, ma, e_vec_norm)
      w = w / e_vec_norm
      ! Calculate the difference between the old eigenvector (eig_vec) and the
      ! new eigenvector (w).
      r = w - e_vec
      ! Calculate the new norm of r. 
      CALL calc_two_norm(r, ma, r_norm)
      IF (r_norm .EQ. 2.0_dp) THEN
        STOP('STOP')
      ENDIF
      ! Finally replace the old eigenvector with the new one. 
      e_vec = w
    ENDDO 
      
  END SUBROUTINE inverse_iteration 
        

  !--------------------------------------------------------
  !----------------Eigenvalue-QR-Algorithm-----------------
  SUBROUTINE qr_algorithm(shift,                          &
  &                       err,                            &
  &                       ma,                             &
  &                       A,                              &
  &                       N_iter)
  !--------------------------------------------------------
    ! This function computes the eigenvalues of the given matrix 
    ! A  using the QR alogrithm either with or without shifts.
    ! 
    ! Parameters
    ! ----------
    ! shift: Boolean
    !     The flag for using a shift or not. If true then the algorithm uses a
    !     shift, if false it does not. 
    ! ma: Integer
    !     The dimension of the square matrix A. 
    ! A: Float, (ma,ma)
    !     The square matrix for which the eigenvalues are solved for.
    ! err: Float, scalar
    !     The convergence tolerance for the algorithm.
    !
    ! Returns
    ! -------
    ! A: Float, (ma,ma)
    !     The diagonal matrix A, with the eigenvalues in the diagonal.
    ! N_iter: Integer
    !     The number of iterations required for convergence.

    ! Global Vars
    ! -----------
    ! The flag to implement shift or not.
    logical, intent(in)                :: shift
    ! The given err tolerance for convergence.
    real (dp), intent(in)              :: err
    ! Dimension of A. 
    integer, intent(in)                :: ma
    ! matrix A
    real (dp), intent(inout)           :: A(ma,ma)
    ! The number of iterations before convergence. 
    integer, intent(out)               :: N_iter

    ! Local Vars
    ! ----------
    ! Indexing. 
    integer                            :: i
    ! The matrix Q. 
    real (dp)                          :: Q(ma,ma)
    ! The subdiagonal vector. 
    real (dp)                          :: subd_vec(ma-1)
    ! The two norm of the subdiagonal vector. 
    real (dp)                          :: subd_norm
    ! This is the shift value. 
    real (dp)                          :: mu
    ! This is the identity matrix for the qr algorithm with shift. 
    real (dp)                          :: Id(ma,ma)
    ! The A-mu*Id matrix. 
    real (dp)                          :: AmmuId(ma,ma)
    ! This is the congergance criterion value for the qr algorithm w/ shift. 
    real (dp)                          :: convg


    ! The start of the algorithm. 
    ! Init Q as zeros. 
    Q = 0.0_dp
    subd_vec = 0.0_dp
    Id = 0.0_dp
    AmmuId = 0.0_dp
    ! Init the subd_norm such that it is greater than the err since it is the
    ! measure of convergence for the case of the qr algorithm without shifts.
    subd_norm = err + 1.0_dp
    ! Init the convg value of the shifted method to be also greater tham err. 
    convg = err +1.0_dp

    ! If the shift flag is false then do not perform shift. 
    IF (shift .EQV. .FALSE.) THEN 
      ! The counter to count the iterations before convergence.
      N_iter=0
      ! Implement the QR algorithm without shift. 
      DO WHILE (subd_norm > err)
        ! Perform the QR decomposition of A. 
        CALL householder_qr_decomp(ma, ma, A, Q)
        ! Solve for A = RQ
        A = MATMUL(A,Q)
        ! Retreieve the subdiagonal vector of A. 
        DO i=1,ma-1
          subd_vec(i) = A(i+1, i)
        ENDDO 
        ! Compute the two norm. 
        CALL calc_two_norm(subd_vec, ma-1, subd_norm) 
        ! increase the counter. 
        N_iter = 1+N_iter
      ENDDO 
    ENDIF

    ! If the shift flag is true then compute using shift.
    IF (shift .EQV. .TRUE.) THEN 
      ! Construc the identity matrix. 
      DO i=1,ma
        Id(i,i) = 1.0_dp
      ENDDO 
      ! Init the counter to count iterations before convergence.
      N_iter = 0 
      ! Loop over iterations until desired convergence. 
      DO WHILE (convg > err)
        ! The shift is the last eigenvalue on the diagonal of A. 
        mu = A(ma,ma)
        ! Calculate A-mu*I. 
        AmmuId = A - mu*Id
        ! Implement the QR algorithm with shift.
        CALL householder_qr_decomp(ma, ma, AmmuId, Q)
        ! solve for A = RQ + mu*Id, at this point AmmuId is R.
        A = MATMUL(AmmuId,Q) + mu*Id
        ! This is the convergence tolderance metric.
        ! Ian mentions this in the announcement on canvas: 
        ! "Hw4: Convergenve criteria for eigenvalue methods"
        ! The final subdiagonal entry. 
        convg = ABS(A(ma, ma-1))
        ! Increase the ccounter. 
        N_iter = 1 +  N_iter
      ENDDO
    ENDIF


  END SUBROUTINE qr_algorithm

  !--------------------------------------------------------
  !----------------Reduction-To-Hessenburg-Form------------
  SUBROUTINE hessenburg_form(ma,                          &
  &                          A)
  !--------------------------------------------------------
    ! This function computes the reduction of the matrix A 
    ! into its hessenburg form. If A is symmetric then the 
    ! resulting hessenburg form should be tridiagonal.
    ! 
    ! Paramters
    ! ---------
    ! ma: Integer
    !     The dimension of the input square matrix A. 
    ! A: Float, (ma,ma)
    !     The input square matrix whose hessenburg form is solved for.
    ! 
    ! Returns
    ! -------
    ! A: Float, (ma,ma)
    !     The output hessenburg form of the given matrix. Again this will be 
    !     in tridiagonal form if the input A is symmetric. 

    ! Global Vars
    ! -----------
    ! The dimension of A. 
    integer, intent(in)                :: ma
    ! The matrix A, returned in hessenburg form. 
    real (dp), intent(inout)           :: A(ma,ma)

    ! Local Vars
    ! ----------
    ! Indexing. 
    integer                            ::i,j,k
    ! The signed norm.
    real (dp)                          :: s 
    ! The householder vector. 
    real (dp)                          :: v(ma)
    ! A place holder for the sum. 
    real (dp)                          :: tot

    ! The start of the reduction to hessenburg form algorithm.

    ! First loop over the columns of A, up to ma-2, since we are putting in
    ! hessenburg form this means we should end the column before last. 
    DO j=1,ma-2
      ! Init the tot and the houyseholder vector to be zeros.
      tot = 0.0_dp
      v = 0.0_dp 
      ! Compute the Signed norm. .
      ! Loop over the rows below the diagonal of current column and add
      ! to calculate the norm. 
      DO k=j+1,ma
        tot = tot + A(k,j) ** 2
      ENDDO 
      ! The signed norm
      ! SIGN(A,B) function returns value A with sign of B.
      s = SIGN(SQRT(tot), A(j+1,j)) 
      ! Compute the normalized householder vector.
      v(j+1:) = A(j+1:,j)
      ! Adding the signed norm to the v(j+1) value. 
      v(j+1) = v(j+1) + s 
      ! Compute the norm of v. 
      CALL calc_two_norm(v, ma, tot)
      ! Normalize the householder vector. 
      v = v / tot

      ! Make the necessary updates to the matrix A. 
      ! Start on the left. 
      A = A - 2*(MATMUL(MATMUL(RESHAPE(v, (/ma,1/)), RESHAPE(v, (/1, ma/))), A))
      ! Update from the right. 
      A = A - 2*(MATMUL(A, MATMUL(RESHAPE(v, (/ma,1/)), RESHAPE(v, (/1, ma/)))))
    ENDDO

  END SUBROUTINE hessenburg_form 

  
  !--------------------------------------------------------
  !--------------Singular-Value-Compression----------------
  SUBROUTINE svd_compression(ma,                          &
  &                          na,                          &
  &                          k_start,                     &
  &                          k_stop,                      &
  &                          U,                           &
  &                          S,                           & 
  &                          VT,                          &
  &                          A)                           
  !--------------------------------------------------------
    ! This function implements a reduced rank representation to create a 
    ! compressed matrix A, from the SVD of its parent matrix.
    ! 
    ! See AM213A Canvas Announcement "Forming approximate images" by Ian May on
    ! March 7, 2022
    ! 
    ! Parameters
    ! ----------
    ! ma: Integer
    !     This is the number of rows of A. 
    ! na: Integer
    !     The number of columns of A. 
    ! k_start: Intger
    !     This is the starting index of the sum, the loop over the singular 
    !     values and vectors will start at this k.
    ! k_stop: Integer
    !     This is the ending index of the sum, the loop over the singular
    !     values and vector will stop at this k. 
    ! U: Float, (ma,ma)
    !     This is the left side orthogonal matrix of the SVD. The columns of
    !     this matrix will be used in the compression. 
    ! S: Float, (na)
    !     This is the singular value vector, its values are the singular values
    !     of the parent matrix of A. 
    ! VT: Float, (na, na)
    !     This is the right hand side transposed orthogonal matrix of the SVD.
    !     the rows of this matrix will be used to calculate the compression. 
    ! A: Float, (ma, na)
    !     This should be input as either a zero matrix or a previous lower rank
    !     approximation of the compression. 
    ! 
    ! Returns
    ! -------
    ! A: Flaot, (ma, na)
    !     This is the resulting compressesion using the given SVD. 

    ! Global Var
    ! ----------
    ! The dimensions of A. 
    integer, intent(in)                :: ma, na
    ! The start and stop indexes. 
    integer, intent(in)                :: k_start, k_stop
    ! The left side matrix of SVD. 
    real (dp), intent(in)              :: U(ma,ma)
    ! Singular Value matrix. 
    real (dp), intent(in)              :: S(na)
    ! The right side matrix of SVD. 
    real (dp), intent(in)              :: VT(na,na)
    ! The provided matrix A
    real (dp), intent(inout)           :: A(ma,na)

    ! Local Vars
    ! ----------
    ! Indexing. 
    integer                            :: i

    ! Loop over the given start, stop indexes. 
    DO i=k_start,k_stop
      ! This is A + the sinle value times the outer product of the ith vector 
      ! of U and the ith column of VT.
      A = A + S(i) * MATMUL(RESHAPE(U(:,i), (/ma,1/)), RESHAPE(VT(i,:), (/1,na/)))
    ENDDO

  END SUBROUTINE svd_compression

  !--------------------------------------------------------
  !----------Gauss-Jacobi-Iteration-Algorithm--------------
  SUBROUTINE gauss_jacobi(ma,                             &
  &                       A,                              &
  &                       b,                              &
  &                       acc,                            &
  &                       savefile,                       &
  &                       x,                              &
  &                       N_iter)
  !---------------------------------------------------------
    ! This subroutine implements the Gauss-Joacobi Iterativive 
    ! method for the solution of a linear system Ax=b where 
    ! A is square. 
    ! 
    ! Parameters
    ! ----------
    ! ma: Integer
    !     The number of rows and columns of the given matrix A. 
    !     Also the length of b and x vectors. 
    ! A: Flaot, (ma,ma)
    !     The given matrix of Ax=b. 
    ! b: Float, (ma)
    !     The right hand side solution vector. 
    ! acc: Float, Scalar
    !     The tolerance for the accuracy of the solution such that 
    !     ||b-Ax||_2 < acc
    ! x: Float, (ma)
    !     This on input is the initial guess for x. 
    ! 
    ! Results
    ! -------
    ! x: Float, (ma)
    !     This is the returned solution for x for the given convergence
    !     tolerance. 

    ! Global Vars
    ! -----------
    ! Dimension of A. 
    integer, intent(in)                :: ma 
    ! The provided matrix, A. 
    real (dp), intent(in)              :: A(ma,ma)
    ! The b vector. 
    real (dp), intent(in)              :: b(ma)
    ! The desired accuracy for convergence. 
    real (dp), intent(in)              :: acc
    ! The name of the file to save the convergence values at each iteration
    ! into. 
    character(len=*), intent(in)       :: savefile
    ! The x vector solution solved for given convergence. 
    real (dp), intent(inout)           :: x(ma)
    ! The number of iterations for convergence. 
    integer, intent(out)                :: N_iter

    ! Local Vars
    ! ----------
    ! Inexing. 
    integer                            :: i,j,k
    ! The place holder for the previous x vector in the iteration. 
    real (dp)                          :: x_prev(ma)
    ! Error vector, b-Ax. 
    real (dp)                          :: e(ma)
    ! The calculated error of ||b-Ax||_2. 
    real (dp)                          :: convg
    ! A sum place holder. 
    real (dp)                          :: tot

    ! Init the initial convg to be > acc. 
    convg = 1 + acc
    
    ! Init the x_prev to the initial guess for x. 
    x_prev = x

    ! Open the file we are writing the convergence value to. 
    OPEN(10,file=savefile)

    ! The convergence index, k
    k=1
    !  While the convergence tolerance is not met. 
    DO WHILE (convg>acc)
      ! Loop over the column. 
      DO i=1,ma
        ! The dot product of x with column of R, ie A without the diagonal. 
        ! This is done with, x(k), the previous iteration solution for x.
        tot=0.0_dp
        DO j=1,ma
          ! All entries are part of the sum except diagonal. 
          IF (j .NE. i) THEN 
            tot = tot + A(i,j) * x_prev(j)
          ENDIF
        ENDDO

        ! Calculate the next value of x. 
        x(i) = (1/A(i,i)) * (b(i) - tot)
      ENDDO

      ! Calculate the error matrix. 
      e = b - MATMUL(A,x)
      ! Calculate the convergence. 
      CALL calc_two_norm(e, ma, convg)

      ! Write the congence to file.
      WRITE(10,*) convg

      ! Update the x_prev vector. 
      x_prev = x

      ! Increase index. 
      k=k+1

      ! If no convergence: 
      IF (k> 1000) THEN
        PRINT *, '!BAD CONVERGENCE!'
        EXIT
      ENDIF 

    ENDDO

    N_iter = k

    ! Lastly close the file. 
    CLOSE(10)

  END SUBROUTINE gauss_jacobi


  !--------------------------------------------------------
  !-------------Gauss-Seidel-Iteration-Algorithm-----------
  SUBROUTINE gauss_seidel(ma,                             &
  &                       A,                              &
  &                       b,                              &
  &                       acc,                            &
  &                       savefile,                       &
  &                       x,                              &
  &                       N_iter)
  !--------------------------------------------------------
    ! This subroutine implements the Gauss-Seidel iterative 
    ! method for the solution of a linear system Ax=b where 
    ! A is square. 
    ! 
    ! The algortihm follows equation 6.15 in the AM213A Lecture notes. 
    ! 
    ! Parameters
    ! ----------
    ! ma: Integer
    !     The number of rows and columns of the given matrix A. 
    !     Also the length of b and x vectors. 
    ! A: Flaot, (ma,ma)
    !     The given matrix of Ax=b. 
    ! b: Float, (ma)
    !     The right hand side solution vector. 
    ! acc: Float, Scalar
    !     The tolerance for the accuracy of the solution such that 
    !     ||b-Ax||_2 < acc
    ! x: Float, (ma)
    !     This on input is the initial guess for x. 
    ! 
    ! Results
    ! -------
    ! x: Float, (ma)
    !     This is the returned solution for x for the given convergence
    !     tolerance. 

    ! Global Vars
    ! -----------
    ! Dimension of A. 
    integer, intent(in)                :: ma 
    ! The provided matrix, A. 
    real (dp), intent(in)              :: A(ma,ma)
    ! The b vector. 
    real (dp), intent(in)              :: b(ma)
    ! The desired accuracy for convergence. 
    real (dp), intent(in)              :: acc
    ! The name of the file to save the convergence values at each iteration
    ! into. 
    character(len=*), intent(in)       :: savefile
    ! The x vector solution solved for given convergence. 
    real (dp), intent(inout)           :: x(ma)
    ! The number of iterations for convergence. 
    integer, intent(out)                :: N_iter

    ! Local Vars
    ! ----------
    ! Inexing. 
    integer                            :: i,j,k
    ! Error vector, b-Ax. 
    real (dp)                          :: e(ma)
    ! The calculated error of ||b-Ax||_2. 
    real (dp)                          :: convg
    ! A sum place holder. 
    real (dp)                          :: tot1
    ! Another seperate sum place holder. 
    real (dp)                          :: tot2


    ! Init the initial convg to be > acc. 
    convg = 1 + acc
    
    ! Open the file we are writing the convergence value to. 
    OPEN(10,file=savefile)

    ! The convergence index, k
    k=1
    !  While the convergence tolerance is not met. 
    DO WHILE (convg>acc)
    
      ! Loop over the entries of the column vector. 
      DO i=1,ma

        ! This is the left sum in equation 6.15, the sum over the newly 
        ! added entries of x. 
        tot1=0.0_dp
        DO j=1,i-1
          ! All entries are part of the sum except diagonal. 
          IF (j .NE. i) THEN 
            tot1 = tot1 + A(i,j) * x(j)
          ENDIF
        ENDDO

        ! This is the right sum in equation 6.15, the sum over the remaining
        ! entries of x. 
        tot2=0.0_dp
        DO j=i+1,ma
          ! All entries are part of the sum except diagonal. 
          IF (j .NE. i) THEN 
            tot2 = tot2 + A(i,j) * x(j)
          ENDIF
        ENDDO

        ! Update the next value of x. 
        x(i) = (1/A(i,i)) * (b(i) - tot1 - tot2)

      ! End Loop over the entries. 
      ENDDO

      ! Calculate the error matrix. 
      e = b - MATMUL(A,x)
      ! Calculate the convergence. 
      CALL calc_two_norm(e, ma, convg)

      ! Write the congence to file.
      WRITE(10,*) convg

      ! Increase index. 
      k=k+1

      ! If no convergence: 
      IF (k> 1000) THEN
        PRINT *, '!BAD CONVERGENCE!'
        EXIT
      ENDIF 

    ENDDO

    N_iter = k

    ! Lastly close the file. 
    CLOSE(10)

  END SUBROUTINE gauss_seidel


  !--------------------------------------------------------
  !---------------Conjugate-Gradient-Algorithm-------------
  SUBROUTINE conj_grad(ma,                                &
  &                    A,                                 &
  &                    b,                                 &
  &                    acc,                               &
  &                    x,                                 &
  &                    N_iter)
  !--------------------------------------------------------
    ! This routine implements the Conjugate Gradient Algorithm
    ! to solve Ax=b for a real, positive defnite, symmetric matrix A. 
    ! 
    ! Parameters
    ! ----------
    ! ma: Integer
    !     The number of rows and columns of the given matrix A. 
    !     Also the length of b and x vectors. 
    ! A: Flaot, (ma,ma)
    !     The given matrix of Ax=b. 
    ! b: Float, (ma)
    !     The right hand side solution vector. 
    ! acc: Float, Scalar
    !     The tolerance for the accuracy of the solution such that 
    !     ||b-Ax||_2 < acc
    ! x: Float, (ma)
    !     This on input is the initial guess for x. 
    ! 
    ! Results
    ! -------
    ! x: Float, (ma)
    !     This is the returned solution for x for the given convergence
    !     tolerance. 

    ! Global Vars
    ! -----------
    ! Dimension of A. 
    integer, intent(in)                :: ma 
    ! The provided matrix, A. 
    real (dp), intent(in)              :: A(ma,ma)
    ! The b vector. 
    real (dp), intent(in)              :: b(ma)
    ! The desired accuracy for convergence. 
    real (dp), intent(in)              :: acc
    ! The name of the file to save the convergence values at each iteration
    ! into. 
    !character(len=*)                   :: savefile
    ! The x vector solution solved for given convergence. 
    real (dp), intent(inout)           :: x(ma)
    ! The number of iterations for convergence. 
    integer, intent(out)               :: N_iter

    ! Local Vars
    ! ----------
    ! Inexing. 
    integer                            :: i,j,k
    ! The Error vector. 
    real (dp)                          :: r(ma)
    ! The square of the two norm of r. 
    real (dp)                          :: E, E_new
    ! The p value in the conjugate gradient method. 
    real (dp)                          :: p(ma)
    ! Some scalar place holders in CG method. 
    real (dp)                          :: alpha, beta
    ! The y vector. 
    real (dp)                          :: y(ma)


    ! Check if A is symmetric, i.e. that A = TRANSPOSE(A). 
    DO i=1,ma
      DO j=1,ma
        IF (A(i,j) .NE. A(j,i)) THEN 
          STOP("A is NOT symmetric")
        ENDIF
      ENDDO
    ENDDO

    ! Initialize r, p
    r = b - MATMUL(A,x)
    p = r 

    ! Calculate the square of the norm of the r vector. 
    CALL calc_two_norm(r, ma, E) 
    E = E**2

    ! The iteration counter. 
    k = 1
    ! Loop while sqrt(E) > desired accuracy. 
    DO WHILE (SQRT(E) > acc)
      ! Calculate y. 
      y = MATMUL(A, p)
      ! Calculate aplha
      alpha = E / DOT_PRODUCT(p, y)
      ! Calculate and update x
      x = x + alpha * p
      ! Calculate r. 
      r = r - alpha*y 
      
      ! Recalculate E. 
      CALL calc_two_norm(r,ma,E_new)
      E_new = E_new**2

      ! calculate beta. 
      beta = E_new / E 

      ! Calculate and update p. 
      p = r + beta * p 

      ! Update E. 
      E = E_new

      ! increase k. 
      k = k+1

    ENDDO

    N_iter = k

  END SUBROUTINE conj_grad


END MODULE LinAl
