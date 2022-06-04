! This is the driver routine for AM213A Final Part 1 Problem a
! 
PROGRAM Driver_SVD_Compression

  USE LinAl, ONLY: read_mat, write_mat
  USE LinAl, ONLY: svd_compression

  implicit none

  ! Parameters
  ! ----------
  integer, parameter                   :: dp=kind(0.d0)
  ! The name of the original image data file to be read in. 
  character(len=100), parameter        :: imdat_file = "dog_bw_data_lowRes_811x1280.dat"
  ! The number of rows of the original image matrix, A 
  integer, parameter                   :: ma=811
  ! The number of columns of the orignal image matrix, A. 
  integer, parameter                   :: na=1280
  ! The number of k values. 
  integer, parameter                   :: N_k = 8
  ! The list of k-values
  integer, parameter, dimension(N_k)   :: k_list=(/20,40,80,160,320,640,1280,2560/)
  ! Test k to text the algorithm without loopin over ks. 
  integer, parameter                   :: test_k=20

  ! Parameters for the implementation of LAPACK's dgesvd. 
  ! -----------------------------------------------------
  ! See reference for complete dexcription: 
  ! http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html
  ! This is the options for ortho. mat. U, 'A' option returns all columns. 
  character(len=1), parameter          :: JOBU='A'
  ! This is the options for ortho. mat. VT, 'A' returns all rows.
  character(len=1), parameter          :: JOBVT='A'
  ! The dimensions of the inputs to dgesvd
  integer, parameter                   :: M=ma
  integer, parameter                   :: N=na
  ! The leading dimension of A. 
  integer, parameter                   :: LDA=M
  ! The leading dimension of U
  integer, parameter                   :: LDU=M
  ! The leading dimension of VT.
  integer, parameter                   :: LDVT=N

  ! Local Variables
  ! ---------------
  ! The im data matrix A. 
  real (dp), allocatable               :: A(:,:)
  ! The saved version of A. 
  real (dp), allocatable               :: A_s(:,:)
  ! The singular values of A, taken of size min(ma,na) 
  real (dp), allocatable               :: S(:)
  ! The left side orthogonal matrix from the SVD. 
  real (dp), allocatable               :: U(:,:)
  ! The transposed right side orthogonal matrix from the SVD. 
  real (dp), allocatable               :: VT(:,:)
  ! The work matrix of the dgesvd method. 
  real (dp), allocatable               :: WORK(:)
  ! LWORK dimension of the array WORK. 
  ! I will be using the dgesvd optimal LWORK routine Call following 
  ! the demonstration of: 
  ! https://www.intel.com/content/www/us/en/develop/documentation/onemkl-lapack-examples/top/least-squares-and-eigenvalue-problems/singular-value-decomposition/gesvd-function/dgesvd-example/dgesvd-example-fortran.html
  integer                              :: LWORK
  ! The returned info for the dgesvd call.
  integer                              :: INFO 
  ! indexing 
  integer                              :: i
  ! The start and stop of k for the svd compression subroutine. 
  integer                              :: k_start, k_stop
  ! The name of the file that the compressed image will be saved into. 
  character(len=100)                   :: savefile
  ! The character value of k to be formatted into the savefile string
  character(len=6)                     :: k_char
  ! This is a flag for reading in the matrix, set false since the 
  ! first line of the data file is not the dimension of the matrix. 
  logical                              :: dims_in_file


  ! Initializations
  ! ---------------

  ! Array allocation and initialization. 
  ALLOCATE(A(ma,na))
  ALLOCATE(A_s(ma,na))
  ALLOCATE(S(N))
  ALLOCATE(U(LDU,M))
  ALLOCATE(VT(LDVT,N))
  ! Let the work array be allocated to size 1 for the moment, until the 
  ! optimal LWORK query retrieves the optimal LWORK. 
  ALLOCATE(WORK(1))
  
  ! Inititalize the arrays to zero. 
  A=0.0_dp 
  A_s=0.0_dp 
  S=0.0_dp
  U=0.0_dp 
  VT=0.0_dp
  WORK=0.0_dp

  ! Now read in A. 
  dims_in_file = .FALSE.
  CALL read_mat(imdat_file, ma, na, dims_in_file, A)
  ! Save copy of A. 
  A_s = A 


  ! Perform Full SVD
  ! ----------------
  
  ! Perform the workspace query for the optimal size of the 
  ! WORK array, LWORK. 
  ! To implement this query we set LWORK =-1
  LWORK=-1
  CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
  print *, 'DGESVD INFO1:', INFO 

  ! Now that I know the optimal size of WORK, I can deallocate and reallocate it
  ! to that optimal size. 
  ! First let LWORK be the optimal value calculated form previous call. 
  LWORK = INT(WORK(1))
  ! Deallocate previous WORK Array. 
  DEALLOCATE(WORK)
  ! Reallocate and initialize to zero using optimal LWORK. 
  ALLOCATE(WORK(LWORK))
  WORK=0.0_dp
  ! Remake A, since it might have been destroyed in previous call.
  A = A_s
  ! Now to perform the actual SVD. 
  CALL DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
  print *, 'DGESVD INFO2:', INFO 
 

  ! Perform Image Compression
  ! -------------------------

  ! Init A zeros. 
  A = 0.0_dp

  ! Init the k_start to be = 1. 
  k_start = 1
  ! Loop over the different k_values, 
  DO i=1,N_k
    ! The particular stoping point of the loop over the singular value should be 
    ! the given k value. 
    k_stop = k_list(i)

    ! The SVD compression routine. 
    CALL svd_compression(ma, na, k_start, k_stop, U, S, VT, A)

    ! Write the resulting compressed image data to file. 
    ! the name of the file to write the data into. 
    WRITE(k_char,'(I6)'), 100000 + k_stop
    savefile = "Image_appn_" // k_char // ".dat"
    dims_in_file = .TRUE.
    CALL write_mat(savefile, ma, na, dims_in_file, A)

    ! Lastly, increase the start point of the next iteration by one. 
    ! The plus one is beacuse fortran inclused both the designated 
    ! start and stop of DO loop, i.e.
    ! DO i=2,5 ==> i = {2,3,4,5}.
    k_start = k_stop + 1

  ENDDO 

  ! Array allocation and initialization. 
  DEALLOCATE(A)
  DEALLOCATE(A_s)
  DEALLOCATE(S)
  DEALLOCATE(U)
  DEALLOCATE(VT)
  DEALLOCATE(WORK)

END PROGRAM Driver_SVD_Compression
