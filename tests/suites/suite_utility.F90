module suite_utility

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_utility

contains

  subroutine collect_suite_utility(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test Hermitian Matrix Exponential", test_expsh), &
                new_unittest("Test Skew-Hermitian Matrix Exponential", test_expsa), &
                new_unittest("Test Unitary Matrix Logarithm", test_logu), &
                new_unittest("Test SVD Decomposition", test_svd) &
                ]

  end subroutine collect_suite_utility

  subroutine test_expsh(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    complex(kind=dp) :: matrixh(2, 2), expm(2, 2)

    logical :: derror = .false.

    matrixh(1, :) = (/cmplx(2.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 2.0_dp, dp)/)
    matrixh(2, :) = (/cmplx(1.0_dp, -2.0_dp, dp), cmplx(3.0_dp, 0.0_dp, dp)/)

    expm = SsTC_utility_exphs(mat=matrixh, dim=2, skew=.false., error=derror)

    if (abs(real(expm(1, 1), dp) - 47.835805811134101) > 1.0E-5_dp) cond = .false.

    if (.not. cond) then
      call test_failed(error, "Mismatch between the results of the test &
                               &for the exponential of a Hermitian matrix and the reference.")
      return
    end if

  end subroutine test_expsh

  subroutine test_expsa(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    complex(kind=dp) :: matrixa(2, 2), expm(2, 2)

    logical :: derror = .false.

    matrixa(1, :) = (/cmplx(2.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 2.0_dp, dp)/)
    matrixa(2, :) = (/cmplx(1.0_dp, -2.0_dp, dp), cmplx(3.0_dp, 0.0_dp, dp)/)
    matrixa = cmplx(0.0_dp, 1.0_dp, dp)*matrixa

    expm = SsTC_utility_exphs(mat=matrixa, dim=2, skew=.true., error=derror)

    if (abs(real(expm(1, 1), dp) - 0.62669928264646291) > 1.0E-6_dp) cond = .false.

    if (.not. cond) then
      call test_failed(error, "Mismatch between the results of the test &
                               &for the exponential of a skew-Hermitian matrix and the reference.")
      return
    end if

  end subroutine test_expsa

  subroutine test_logu(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    complex(kind=dp) :: matrixa(2, 2), expm(2, 2), logu(2, 2)

    logical :: derror = .false.

    matrixa(1, :) = (/cmplx(2.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 2.0_dp, dp)/)
    matrixa(2, :) = (/cmplx(1.0_dp, -2.0_dp, dp), cmplx(3.0_dp, 0.0_dp, dp)/)
    matrixa = cmplx(0.0_dp, 1.0_dp, dp)*matrixa/10.0_dp

    expm = SsTC_utility_exphs(mat=matrixa, dim=2, skew=.true., error=derror)

    logu = SsTC_utility_logu(mat=expm, dim=2, error=derror)

    if (abs(aimag(logu(1, 1)) - aimag(matrixa(1, 1))) > 1.0E-6_dp) cond = .false.

    if (.not. cond) then
      call test_failed(error, "Mismatch between the results of the test &
                               &for the logarithm of a unitary matrix and the reference.")
      return
    end if

  end subroutine test_logu

  subroutine test_svd(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    complex(kind=dp) :: matrix(4, 5) = cmplx(0.0_dp, 0.0_dp, dp)
    real(kind=dp) :: sigma(4, 5)
    complex(kind=dp) :: u(4, 4), v(5, 5), alt(4, 5)

    logical :: derror = .false.

    matrix(1, 1) = cmplx(1.0_dp, 0.0_dp, dp)
    matrix(1, 5) = cmplx(2.0_dp, 0.0_dp, dp)
    matrix(2, 3) = cmplx(3.0_dp, 0.0_dp, dp)
    matrix(4, 2) = cmplx(2.0_dp, 0.0_dp, dp)

    call SsTC_utility_SVD(mat=matrix, sigma=sigma, error=derror, U=u, V=v)
    alt = matmul(u, matmul(sigma, conjg(transpose(v))))

    if (abs(sigma(2, 2) - sqrt(5.0_dp)) > 1.0E-6_dp) cond = .false.
    if (abs(real(alt(4, 2), dp) - 2.0_dp) > 1.0E-6_dp) cond = .false.

    if (.not. cond) then
      call test_failed(error, "Mismatch between the results of the test &
                               &for the SVD decomposition of a matrix and the reference.")
      return
    end if

  end subroutine test_svd

end module suite_utility
