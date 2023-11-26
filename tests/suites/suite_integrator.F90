module suite_integrator

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_integrator

contains

  subroutine collect_suite_integrator(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test Integrator: Rectangle Method", test_R), &
                new_unittest("Test Integrator: Extrapolation Method", test_E) &
                ]

  end subroutine collect_suite_integrator

  subroutine test_R(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    real(kind=dp) :: res
    logical :: cond = .true.

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    global_calculator:block

      !g_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_BZ_integral_task) :: test_integral

      call SsTC_BZ_integral_task_constructor(task=test_integral, name="g_rectangle_test", &
                                             g_calculator=g_calc_test, &
                                             method="rectangle", samples=(/33, 1, 1/), &
                                             N_int_ind=1, int_ind_range=(/1/), &
                                             N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                             ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_and_integrate_BZ_integral_task(task=test_integral, &
                                                      system=SsTC_toy_model)

      call SsTC_print_BZ_integral_task(task=test_integral, &
                                       system=SsTC_toy_model)

      res = real(test_integral%result(1, 16), dp) !Retrieve sinc(0.30303030E+001) in the rectangle appr.

      if (abs(res - 5.3173491865169901E-003_dp) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the integration &
                                 &of a global calculator using the rectangle &
                                 &method and the reference.")
        return
      end if

    end block global_calculator

    cond = .true.

    local_calculator:block

      !l_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_BZ_integral_task) :: test_integral

      call SsTC_BZ_integral_task_constructor(task=test_integral, name="l_rectangle_test", &
                                             l_calculator=l_calc_test, &
                                             method="rectangle", samples=(/33, 1, 1/), &
                                             N_int_ind=1, int_ind_range=(/10/), &
                                             N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                             ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_and_integrate_BZ_integral_task(task=test_integral, &
                                                      system=SsTC_toy_model)

      call SsTC_print_BZ_integral_task(task=test_integral, &
                                       system=SsTC_toy_model)

      res = real(test_integral%result(7, 1), dp) !Retrieve sinc(7) in the rectangle appr.

      if (abs(res - 0.11240032628276973) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the integration &
                                 &of a local calculator using the rectangle &
                                 &method and the reference.")
        return
      end if

    end block local_calculator

  end subroutine test_R

  subroutine test_E(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    real(kind=dp) :: res
    logical :: cond = .true.

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    global_calculator:block

      !g_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_BZ_integral_task) :: test_integral

      call SsTC_BZ_integral_task_constructor(task=test_integral, name="extrapolation_test", &
                                             g_calculator=g_calc_test, &
                                             method="extrapolation", samples=(/33, 1, 1/), &
                                             N_int_ind=1, int_ind_range=(/1/), &
                                             N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                             ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_and_integrate_BZ_integral_task(task=test_integral, &
                                                      system=SsTC_toy_model)

      call SsTC_print_BZ_integral_task(task=test_integral, &
                                       system=SsTC_toy_model)

      res = real(test_integral%result(1, 16), dp) !Retrieve sinc(0.30303030E+001) in the extrapolation appr.

      if (abs(res - 3.6649812513784630E-002) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the integration &
                                 &of a global calculator using the extrapolation &
                                 &method and the reference.")
        return
      end if

    end block global_calculator

    cond = .true.

    local_calculator:block

      !l_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_BZ_integral_task) :: test_integral

      call SsTC_BZ_integral_task_constructor(task=test_integral, name="l_extrapolation_test", &
                                             l_calculator=l_calc_test, &
                                             method="extrapolation", samples=(/33, 1, 1/), &
                                             N_int_ind=1, int_ind_range=(/10/), &
                                             N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                             ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_and_integrate_BZ_integral_task(task=test_integral, &
                                                      system=SsTC_toy_model)

      call SsTC_print_BZ_integral_task(task=test_integral, &
                                       system=SsTC_toy_model)

      res = real(test_integral%result(7, 1), dp) !Retrieve sinc(7) in the extrapolation appr.

      if (abs(res - 9.3858749167900907E-002) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the integration &
                                 &of a local calculator using the extrapolation &
                                 &method and the reference.")
        return
      end if

    end block local_calculator

  end subroutine test_E

  function g_calc_test(task, system, k, error) result(u)

    class(SsTC_global_k_data), intent(in) :: task
    type(SsTC_sys), intent(in)            :: system
    real(kind=dp), intent(in)             :: k(3)
    logical, intent(inout)                :: error

    complex(kind=dp) :: u(product(task%integer_indices), &
                          product(task%continuous_indices))

    integer :: i, r, &
               i_arr(size(task%integer_indices)), &
               r_arr(size(task%continuous_indices)), &
               nb

    !The integral of this function is the (unnormalized) sinc(x) function,
    !the role of x is played by task%ext_var_data(1)%data(r_arr(1)).

    u = cmplx(0.0, 0.0, dp)

    nb = system%num_bands

    do i = 1, product(task%integer_indices)
      i_arr = SsTC_integer_memory_element_to_array_element(task, i)
      do r = 1, product(task%continuous_indices)
        r_arr = SsTC_continuous_memory_element_to_array_element(task, r)
        u(i, r) = exp(2.0_dp*cmplx(0.0_dp, 1.0_dp, dp)*k(1)*task%ext_var_data(1)%data(r_arr(1)))
      enddo
    enddo

    error = .false.

  end function g_calc_test

  function l_calc_test(task, system, k, error) result(u)

    class(SsTC_local_k_data), intent(in) :: task
    type(SsTC_sys), intent(in)            :: system
    real(kind=dp), intent(in)             :: k(3)
    logical, intent(inout)                :: error

    complex(kind=dp) :: u(product(task%integer_indices))

    integer :: i, &
               i_arr(size(task%integer_indices)), &
               nb

    !The integral of this function is the (unnormalized) sinc(x) function,
    !the role of x is played by i_arr(1).

    u = cmplx(0.0, 0.0, dp)

    nb = system%num_bands

    do i = 1, product(task%integer_indices)
      i_arr = SsTC_integer_memory_element_to_array_element(task, i)
      u(i) = exp(2.0_dp*cmplx(0.0_dp, 1.0_dp, dp)*k(1)*i_arr(1))
    enddo

    error = .false.

  end function l_calc_test

end module suite_integrator
