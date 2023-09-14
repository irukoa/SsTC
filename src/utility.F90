module SsTC_utility

  use, intrinsic :: iso_fortran_env, only: input_unit, &
    output_unit, &
    error_unit

  implicit none

  private

  !Default units for input, output and error.
  integer, public :: stdin = input_unit, &
                     stdout = output_unit, &
                     stderr = error_unit

  integer, parameter, public          :: dp = 8

  integer, public :: timing(8)

  real(kind=dp), parameter, public    :: pi = acos(-1.0_dp)

  complex(kind=dp), parameter, public :: cmplx_0 = cmplx(0.0_dp, 0.0_dp, dp), &
                                         cmplx_i = cmplx(0.0_dp, 1.0_dp, dp)

  !Physical constants.
  real(kind=dp), parameter, public :: e_charge = 1.602176565e-19_dp, &
                                      e_mass = 9.10938291e-31_dp, &
                                      hbar = 1.054571726e-34_dp, &
                                      k_B = 1.3806488e-23_dp, &
                                      bohr_mag = 927.400968e-26_dp, &
                                      eps0 = 8.854187817e-12_dp, &
                                      c_light = 299792458.0_dp, &
                                      hbar_over_e = 6.582119e-16_dp

  !Symmetrization and antisymmetrization.
  integer, dimension(6), parameter, public :: SsTC_alpha_S = (/1, 2, 3, 1, 1, 2/)
  integer, dimension(6), parameter, public :: SsTC_beta_S = (/1, 2, 3, 2, 3, 3/)

  integer, dimension(3), parameter, public :: SsTC_alpha_A = (/2, 3, 1/)
  integer, dimension(3), parameter, public :: SsTC_beta_A = (/3, 1, 2/)

  public :: SsTC_utility_delta
  public :: SsTC_utility_delta_vec
  public :: SsTC_utility_get_degen
  public :: SsTC_utility_diagonalize
  public :: SsTC_utility_schur
  public :: SsTC_utility_SVD
  public :: SsTC_utility_exphs
  public :: SsTC_utility_logu

contains

  !=============================================================!
  function SsTC_utility_delta(x) result(res)
    !========================================================================!
    !Auxiliary routine to approximate the Dirac delta of x.                  !
    !========================================================================!

    implicit none

    real(kind=dp), intent(in) :: x
    real(kind=dp) :: arg
    real(kind=dp) :: res

    arg = min(200.0_dp, x**2)
    res = exp(-arg)/sqrt(pi)

  end function SsTC_utility_delta

  function SsTC_utility_delta_vec(x) result(res)
    !========================================================================!
    !Auxiliary routine to approximate the Dirac delta of an array x.         !
    !For each x_i returns res_i = \delta(x_i).                               !
    !========================================================================!

    implicit none

    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: arg(size(x))
    real(kind=dp) :: res(size(x))

    arg = min(200.0_dp, x**2)
    res = exp(-arg)/sqrt(pi)

  end function SsTC_utility_delta_vec

  !=============================================================!
  function SsTC_utility_get_degen(eig, degen_thr) result(deg)
    !========================================================================!
    !Auxiliary routine to get the degree of degeneracy for a given list eig. !
    !The list is supposed to have it's elements stored in ascending order:   !
    !eig(i + 1) >= eig(i).                                                   !
    !The degeneracy threshold is given by degen_thr.                         !
    !The results are set up such that if deg(i) = N,                         !
    !then eig(i) = eig(i + 1) = ... = eig(i + N - 1)                         !
    !with the equality holding up to degen_thr.                              !
    !If i < j < i + N - 1, then deg(j) = 0.                                  !
    !If the value is nondegenerate, then deg(j) = 1.                         !
    !========================================================================!

    implicit none

    real(kind=dp), intent(in)  :: degen_thr
    real(kind=dp), intent(in)  :: eig(:)

    integer :: deg(size(eig))

    integer                    :: i, j, dim

    deg = 0
    dim = size(eig)

    do i = 1, dim
      do j = i, dim !In ascending order,
        if (abs(eig(j) - eig(i)) .LE. degen_thr) then
          !count number of elements equal to eig(i).
          deg(i) = deg(i) + 1
        endif
      enddo
    enddo

    !In the case of eig(j+2) - eig(j+1) < degen_thr and eig(j+1) - eig(j) < degen_thr,
    !but eig(j+2) - eig(j) > degen_thr (closely packed levels),
    do i = dim - 1, 2, -1
      if ((deg(i) .GT. 1) .AND. (deg(i - 1) .GT. 1)) then
        deg(i - 1) = deg(i) + 1   !increase the degeneracy value of the
        !degenerate level according to the degenerate levels following it.
        deg(i + 1) = 0            !Set the next levels to 0.
      endif
    enddo

    do i = dim - 1, 1, -1
      !At this point, the second index of a closely packed level shall be corrected.
      if ((deg(i) .NE. 1) .AND. (deg(i) .GE. deg(i + 1)) .AND. (deg(i + 1) .GE. 1)) then
        deg(i + 1) = 0
      endif
    enddo

  end function SsTC_utility_get_degen

  !===============WRAPPERS FOR SOME LAPACK ROUTINES===============!

  !===========================================================!
  subroutine SsTC_utility_diagonalize(mat, dim, eig, rot, error)
    !==================================================================!
    !                                                                  !
    !!Given a hermitian dim x dim matrix mat, computes its             !
    !diagonalization mat = rot*D*rot^dagger, where                     !
    !D_nm = delta_nm eig_n for eig_n real.                             !
    !                                                                  !
    !==================================================================!

    implicit none

    integer, intent(in)     :: dim
    complex*16, intent(in)  :: mat(dim, dim)
    real(8), intent(out)    :: eig(dim)      !Eigenvalues.
    complex*16, intent(out) :: rot(dim, dim) !Eigenvectors.
    logical, intent(inout)  :: error

    complex*16, allocatable :: work(:)
    real(8)                 :: rwork(3*dim - 2)
    integer                 :: info, lwork
    external                :: zheev

    !Initialization.
    rot = mat

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call zheev('V', 'U', dim, rot, dim, eig, work, lwork, rwork, info)
    lwork = nint(real(work(1), 8))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call zheev('V', 'U', dim, rot, dim, eig, work, lwork, rwork, info)
    deallocate (work)

    !Check convergence.
    if (info < 0) then
      error = .True.
      write (unit=stderr, fmt='(a, i3, a)') 'Error in utility_diagonalize: THE ', -info, &
        ' ARGUMENT OF ZHEEV HAD AN ILLEGAL VALUE.'
      return
    endif
    if (info > 0) then
      error = .True.
      write (unit=stderr, fmt='(a, i3, a)') 'Error in utility_diagonalize: ', info, &
        ' EIGENVECTORS FAILED TO CONVERGE.'
      return
    endif

  end subroutine SsTC_utility_diagonalize

  !===========================================================!
  subroutine SsTC_utility_schur(mat, dim, T, Z, error, S)
    !==================================================================!
    !                                                                  !
    !!Given a dim x dim matrix mat, computes its Schur decomposition   !
    !mat = Z*S*Z^dagger.                                               !
    !!Note that the Schur decomposition of unitary matrices always     !
    !involves a diagonal Schur form S, S_nm = delta_nm T_n.            !
    !                                                                  !
    !==================================================================!

    implicit none

    integer, intent(in)               :: dim
    complex*16, intent(in)            :: mat(dim, dim)
    complex*16, intent(out)           :: T(dim)      !Eigenvalues (diagonal elements of B).
    complex*16, intent(out)           :: Z(dim, dim) !Eigenvectors.
    logical, intent(inout)            :: error
    complex*16, intent(out), optional :: S(dim, dim) !Schur Form.

    complex*16              :: B(dim, dim)
    complex*16, allocatable :: work(:)
    real(8)                 :: rwork(dim)
    integer                 :: info, lwork, sdim
    logical                 :: bwork(dim), select
    external                :: zgees

    !Initialization.
    B = mat

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call zgees('V', 'N', select, dim, B, dim, sdim, T, Z, dim, work, lwork, rwork, bwork, info)
    lwork = nint(real(work(1), 8))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call zgees('V', 'N', select, dim, B, dim, sdim, T, Z, dim, work, lwork, rwork, bwork, info)
    if (present(S)) S = B
    deallocate (work)

    !Check convergence.
    if (info < 0) then
      error = .True.
      write (unit=stderr, fmt='(a, i3, a)') 'Error in utility_schur: THE ', -info, &
        ' ARGUMENT OF ZGEES HAD AN ILLEGAL VALUE.'
      return
    endif
    if (info > 0) then
      error = .True.
      write (unit=stderr, fmt='(a, i3,a)') 'Error in utility_schur: ', info, &
        ' EIGENVECTORS FAILED TO CONVERGE.'
      return
    endif

  end subroutine SsTC_utility_schur

  subroutine SsTC_utility_SVD(mat, sigma, error, U, V)
    !==================================================================!
    !                                                                  !
    !!Given a m x n matrix mat, computes its Singular value            !
    !decomposition mat = U*sigma*V^dagger.                             !
    !                                                                  !
    !==================================================================!

    implicit none

    complex*16, intent(in)            :: mat(:, :)
    real(8), intent(out)              :: sigma(size(mat(:, 1)), size(mat(1, :)))
    logical, intent(inout)            :: error
    complex*16, intent(out), optional :: U(size(mat(:, 1)), size(mat(:, 1))), &
                                         V(size(mat(1, :)), size(mat(1, :)))

    complex*16              :: B(size(mat(:, 1)), size(mat(1, :)))
    real(8)                 :: sigmaw(min(size(mat(:, 1)), size(mat(1, :))))
    complex*16              :: Uw(size(mat(:, 1)), size(mat(:, 1))), &
                               Vw(size(mat(1, :)), size(mat(1, :)))
    complex*16, allocatable :: work(:)
    real(8)                 :: rwork(5*min(size(mat(1, :)), size(mat(:, 1))))
    external                :: zgesvd

    integer :: m, n
    integer :: info, lwork
    integer :: i

    m = size(mat(:, 1)) !Number of rows of the input matrix.
    n = size(mat(1, :)) !Number of cols of the input matrix.

    !Initialization.
    B = mat
    sigma = 0.d0

    !Query optimal workspace.
    lwork = -1
    allocate (work(1))
    call zgesvd('A', 'A', m, n, B, m, sigmaw, Uw, m, Vw, n, work, lwork, rwork, info)
    lwork = nint(real(work(1), 8))
    deallocate (work)

    !Calculation.
    allocate (work(lwork))
    call zgesvd('A', 'A', m, n, B, m, sigmaw, Uw, m, Vw, n, work, lwork, rwork, info)
    forall (i=1:size(sigmaw)) sigma(i, i) = sigmaw(i)
    if (present(U)) U = Uw
    if (present(V)) V = conjg(transpose(Vw))
    deallocate (work)

    !Check convergence.
    if (info < 0) then
      error = .True.
      write (unit=stderr, fmt='(a, i3, a)') 'Error in utility_svd: THE ', -info, &
        ' ARGUMENT OF CGESVD HAD AN ILLEGAL VALUE.'
      return
    endif
    if (info > 0) then
      error = .True.
      write (unit=stderr, fmt='(a, i3,a)') 'Error in utility_svd: ', info, &
        ' SUPERDIAGONALS FAILED TO CONVERGE TO ZERO.'
      return
    endif

  end subroutine SsTC_utility_SVD

  !========HERMITIAN/ANTIHERMITIAN MATRIX EXPONENTIAL AND UNITARY MATRIX LOGARITHM========!

  function SsTC_utility_exphs(mat, dim, skew, error) result(exphs)
    !==================================================================!
    !                                                                  !
    !Given a Hermitian/Skew-Hermitian dim x dim matrix mat, the routine!
    !computes the dim x dim matrix exphs such that exphs = exp(mat).   !
    !If mat is Skew-Hermitian, then skew = .true. and if Hermitian,    !
    !then, , then skew = .false..                                      !
    !                                                                  !
    !==================================================================!

    implicit none

    integer, intent(in)          :: dim
    complex(kind=dp), intent(in) :: mat(dim, dim)
    logical, intent(in)          :: skew
    logical, intent(inout)       :: error

    complex(kind=dp) :: exphs(dim, dim), &
                        rot(dim, dim)
    real(kind=dp)    :: eig(dim)
    integer          :: i

    exphs = cmplx_0

    if (skew) then
      !Skew-Hermitian matrix.
      exphs = cmplx_i*mat !Now exphs is Hermitian.

      call SsTC_utility_diagonalize(exphs, dim, eig, rot, error)
      if (error) then
        write (unit=stderr, fmt="(a)") "Error in utility_exphs."
        return
      endif
      exphs = cmplx_0

      do i = 1, dim
        exphs(i, i) = exp(-cmplx_i*eig(i))
      enddo

      exphs = matmul(matmul(rot, exphs), conjg(transpose(rot)))

    else
      !Hermitian matrix.

      call SsTC_utility_diagonalize(mat, dim, eig, rot, error)
      if (error) then
        write (unit=stderr, fmt="(a)") "Error in utility_exphs."
        return
      endif

      do i = 1, dim
        exphs(i, i) = exp(eig(i))
      enddo

      exphs = matmul(matmul(rot, exphs), conjg(transpose(rot)))

    endif

  end function SsTC_utility_exphs

  function SsTC_utility_logu(mat, dim, error) result(logu)
    !==================================================================!
    !                                                                  !
    !Given an Unitary dim x dim matrix mat, computes the Skew-Hermitian!
    !dim x dim matrix logu such that logu = log(mat).                  !
    !                                                                  !
    !==================================================================!

    implicit none

    integer, intent(in)          :: dim
    complex(kind=dp), intent(in) :: mat(dim, dim)
    logical, intent(inout)       :: error

    complex(kind=dp) :: logu(dim, dim), &
                        rot(dim, dim), &
                        eig(dim)
    integer          :: i

    logu = 0.d0

    call SsTC_utility_schur(mat, dim, eig, rot, error)
    if (error) then
      write (unit=stderr, fmt="(a)") "Error in utility_logu."
      return
    endif

    do i = 1, dim
      logu(i, i) = log(eig(i))
    enddo

    logu = matmul(matmul(rot, logu), conjg(transpose(rot)))

  end function SsTC_utility_logu

end module SsTC_utility
