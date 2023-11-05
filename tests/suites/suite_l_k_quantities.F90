module suite_l_k_quantities

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_l_k_quantities

contains

  subroutine collect_suite_l_k_quantities(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test Wannier Hamiltonian", test_WH) &
                ]

  end subroutine collect_suite_l_k_quantities

  subroutine test_WH(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    complex(kind=dp), allocatable :: H(:, :), rot(:, :)
    real(kind=dp), allocatable :: eig(:)

    logical :: cond = .true., &
               derror = .false.

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    allocate (H(SsTC_toy_model%num_bands, SsTC_toy_model%num_bands), &
              rot(SsTC_toy_model%num_bands, SsTC_toy_model%num_bands), &
              eig(SsTC_toy_model%num_bands))
    rot = cmplx(0.0_dp, 0.0_dp, dp)
    eig = 0.0_dp
    H = SsTC_wannier_hamiltonian(SsTC_toy_model, (/0.0_dp, 0.0_dp, 0.0_dp/))
    call SsTC_utility_diagonalize(H, SsTC_toy_model%num_bands, eig, rot, derror)
    if ((abs(eig(1) + 0.5_dp*sqrt(37.0_dp)) .gt. 1.0E-6_dp) .or. (abs(eig(2) - 0.5_dp*sqrt(37.0_dp)) .gt. 1.0E-6_dp)) cond = .false.
    if (derror) cond = .false.

    deallocate (H)

    if (.not. cond) then
      call test_failed(error, "Mismatch between the eigenvalues of the Hamiltonian and the reference.")
      return
    end if

  end subroutine test_WH

end module suite_l_k_quantities
