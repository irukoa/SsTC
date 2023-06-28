module utility

  implicit none

  integer, parameter, public          :: dp = 8

  real(kind=dp), parameter, public    :: pi = acos(-1.0_dp)

  complex(kind=dp), parameter, public :: cmplx_0 = cmplx(0.0_dp, 0.0_dp, dp), &
                                         cmplx_i = cmplx(0.0_dp, 1.0_dp, dp)

  public :: utility_diagonalize
  public :: utility_schur
  public :: utility_exphs
  public :: utility_logu

  contains

  !===============WRAPPERS FOR SOME LAPACK ROUTINES===============!

  !===========================================================!
  subroutine utility_diagonalize(mat, dim, eig, rot, error, errormsg)
    !==================================================================!
    !                                                                  !
    !!Given a hermitian dim x dim matrix mat, computes its             !
    !diagonalization mat = rot*D*rot^dagger, where                     !
    !D_nm = delta_nm eig_n for eig_n real.                             !
    !                                                                  !
    !==================================================================!

    logical,              intent(inout)            :: error
    character(len=120),   intent(inout)            :: errormsg

    integer,              intent(in)               :: dim
    complex*16,           intent(in)               :: mat(dim,dim)
    real(8),              intent(out)              :: eig(dim)      !Eigenvalues.
    complex*16,           intent(out)              :: rot(dim,dim)  !Eigenvectors.

    complex*16,                        allocatable :: work(:)
    real(8)                                        :: rwork(3*dim-2)
    integer                                        :: info, lwork

    external                                       :: zheev

    !Initialization.
    rot = mat

    !Query optimal workspace.
    lwork= -1
    allocate(work(1))
    call zheev('V','U',dim,rot,dim,eig,work,lwork,rwork,info)
    lwork = nint(real(work(1),8))
    deallocate(work)

    !Calculation.
    allocate(work(lwork))
    call zheev('V','U',dim,rot,dim,eig,work,lwork,rwork,info)
    deallocate(work)

    !Check convergence.
    if (info < 0) then
      error = .True.
      write (errormsg, '(a,i3,a)') 'Error in utility_diagonalize: THE ', -info, &
        ' ARGUMENT OF ZHEEV HAD AN ILLEGAL VALUE'
      return
    endif
    if (info > 0) then
      error = .True.
      write (errormsg, '(i3,a)') 'Error in utility_diagonalize: ', info, &
        ' EIGENVECTORS FAILED TO CONVERGE'
      return
    endif

  end subroutine utility_diagonalize

  !===========================================================!
  subroutine utility_schur(mat, dim, T, Z, error, errormsg, S)
    !==================================================================!
    !                                                                  !
    !!Given a dim x dim matrix mat, computes its Schur decomposition   !
    !mat = Z*S*Z^dagger.                                               !
    !!Note that the Schur decomposition of unitary matrices always     !
    !involves a diagonal Schur form S, S_nm = delta_nm T_n.            !
    !                                                                  !
    !==================================================================!

    logical,              intent(inout)            :: error
    character(len=120),   intent(inout)            :: errormsg

    integer,              intent(in)               :: dim
    complex*16,           intent(in)               :: mat(dim,dim)
    complex*16,           intent(out)              :: T(dim)      !Eigenvalues (diagonal elements of B).
    complex*16,           intent(out)              :: Z(dim,dim)  !Eigenvectors.
    complex*16, optional, intent(out)              :: S(dim,dim)  !Schur Form.

    complex*16                                     :: B(dim,dim)
    complex*16,                        allocatable :: work(:)
    real(8)                                        :: rwork(dim)
    integer                                        :: info, lwork, sdim
    logical                                        :: bwork(dim), select

    external                                       :: zgees

    !Initialization.
    B = mat

    !Query optimal workspace.
    lwork= -1
    allocate(work(1))
    call zgees('V','N',select,dim,B,dim,sdim,T,Z,dim,work,lwork,rwork,bwork,info)
    lwork = nint(real(work(1),8))
    deallocate(work)

    !Calculation.
    allocate(work(lwork))
    call zgees('V','N',select,dim,B,dim,sdim,T,Z,dim,work,lwork,rwork,bwork,info)
    if (present(S)) S = B
    deallocate(work)

    !Check convergence.
    if (info < 0) then
      error = .True.
      write (errormsg, '(a,i3,a)') 'Error in utility_schur: THE ', -info, &
        ' ARGUMENT OF ZGEES HAD AN ILLEGAL VALUE'
      return
    endif
    if (info > 0) then
      error = .True.
      write (errormsg, '(i3,a)') 'Error in utility_schur: ', info, &
        ' EIGENVECTORS FAILED TO CONVERGE'
      return
    endif

  end subroutine utility_schur

  !========HERMITIAN MATRIX EXPONENTIAL AND UNITARY MATRIX LOGARITHM========!

  function utility_exphs(mat, dim, skew, error, errormsg) result(exphs)
    !==================================================================!
    !                                                                  !
    !Given a Hermitian/Skew-Hermitian dim x dim matrix mat, the routine!
    !computes the dim x dim matrix exphs such that exphs = exp(mat).   !
    !If mat is Skew-Hermitian, then skew = .true. and if Hermitian,    !
    !then, , then skew = .false..                                      !
    !                                                                  !
    !==================================================================!

    logical,              intent(inout)            :: error
    character(len=120),   intent(inout)            :: errormsg

    integer,              intent(in)               :: dim
    complex(kind=dp),     intent(in)               :: mat(dim, dim)
    logical,              intent(in)               :: skew

    complex(kind=dp)                               :: exphs(dim, dim), &
                                                      rot(dim, dim)
    real(kind=dp)                                  :: eig(dim)
    integer                                        :: i

    exphs = cmplx_0

    if (skew) then
      !Skew-Hermitian matrix.
      exphs = cmplx_i*mat !Now exphs is Hermitian.

      call utility_diagonalize(exphs, dim, eig, rot, error, errormsg)
      if (error .eqv. .True.) return
      exphs = cmplx_0

      do i = 1, dim
        exphs(i,i) = exp(-cmplx_i*eig(i))
      enddo

      exphs = matmul(matmul(rot, exphs), conjg(transpose(rot)))

    else
      !Hermitian matrix.

      call utility_diagonalize(mat, dim, eig, rot, error, errormsg)
      if (error .eqv. .True.) return

      do i = 1, dim
        exphs(i,i) = exp(eig(i))
      enddo

      exphs = matmul(matmul(rot, exphs),conjg(transpose(rot))) 

    endif

  end function utility_exphs

  function utility_logu(mat, dim, error, errormsg) result(logu)
    !==================================================================!
    !                                                                  !
    !Given an Unitary dim x dim matrix mat, computes the Skew-Hermitian!
    !dim x dim matrix logu such that logu = log(mat).                  !
    !                                                                  !
    !==================================================================! 

    logical,              intent(inout)            :: error
    character(len=120),   intent(inout)            :: errormsg

    integer,              intent(in)               :: dim
    complex(kind=dp),     intent(in)               :: mat(dim, dim)

    complex(kind=dp)                               :: logu(dim,dim), &
                                                      rot(dim,dim), eig(dim)
    integer                                        :: i

    logu = 0.d0

    call utility_schur(mat, dim, eig, rot, error, errormsg)
    if (error .eqv. .True.) return

    do i = 1, dim
      logu(i,i) = log(eig(i))
    enddo

    logu = matmul(matmul(rot, logu),conjg(transpose(rot))) 

  end function utility_logu

end module utility