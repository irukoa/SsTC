module suite_sampler

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_sampler

contains

  subroutine collect_suite_sampler(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test Sampler: Regular Sampling", test_R), &
                new_unittest("Test Sampler: Predefined Sampling", test_P) &
                ]

  end subroutine collect_suite_sampler

  subroutine test_R(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    real(kind=dp) :: res
    logical :: cond = .true.

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    global_calculator:block

      !g_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_sampling_task) :: test_sampling

      call SsTC_sampling_task_constructor(task=test_sampling, name="gr_sampling_test", &
                                          g_calculator=g_calc_test, &
                                          samples=(/5, 5, 5/), &
                                          N_int_ind=1, int_ind_range=(/1/), &
                                          N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                          ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_sampling_task(task=test_sampling, &
                                     system=SsTC_toy_model)

      call SsTC_print_sampling(task=test_sampling, &
                               system=SsTC_toy_model)

      res = real(test_sampling%BZ_data(1, 16, 2, 1, 1), dp) !Retrieve Re[exp(2i*(-0.25)*0.30303030E+001)]

      if (abs(res - 5.5616100165806738E-002) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                 &of a global calculator and the reference for regular sampling.")
        return
      end if

    end block global_calculator

    cond = .true.

    local_calculator:block

      type(SsTC_sampling_task) :: test_sampling

      call SsTC_sampling_task_constructor(task=test_sampling, name="lr_sampling_test", &
                                          l_calculator=l_calc_test, &
                                          samples=(/5, 5, 5/), &
                                          N_int_ind=1, int_ind_range=(/10/), &
                                          N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                          ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_sampling_task(task=test_sampling, &
                                     system=SsTC_toy_model)

      call SsTC_print_sampling(task=test_sampling, &
                               system=SsTC_toy_model)

      res = real(test_sampling%BZ_data(7, 1, 2, 1, 1), dp) !Retrieve Re[exp(2i*(-0.25)*7)]

      if (abs(res + 0.93645668729079634) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                &of a local calculator and the reference for regular sampling.")
        return
      end if

    end block local_calculator

  end subroutine test_R

  subroutine test_P(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    real(kind=dp) :: vecs(5, 3)
    integer :: i

    real(kind=dp) :: res
    logical :: cond = .true.

    do i = 1, size(vecs(:, 1))
      vecs(i, :) = (/0.25_dp, 0.25_dp, 0.25_dp/)
      vecs(i, 1) = -0.5_dp + real(i - 1, dp)/real(size(vecs(:, 1)) - 1, dp)
    enddo

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    global_calculator:block

      !g_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_sampling_task) :: test_sampling

      call SsTC_sampling_task_constructor(task=test_sampling, name="gp_sampling_test", &
                                          g_calculator=g_calc_test, &
                                          nkpts=size(vecs(:, 1)), kpts=vecs, &
                                          N_int_ind=1, int_ind_range=(/1/), &
                                          N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                          ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_sampling_task(task=test_sampling, &
                                     system=SsTC_toy_model)

      call SsTC_print_sampling(task=test_sampling, &
                               system=SsTC_toy_model)

      res = real(test_sampling%predefined_sampled_data(1, 16, 2), dp) !Retrieve Re[exp(2i*(-0.25)*0.30303030E+001)]

      if (abs(res - 5.5616100165806738E-002) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                 &of a global calculator and the reference for predefined sampling.")
        return
      end if

    end block global_calculator

    cond = .true.

    local_calculator:block

      type(SsTC_sampling_task) :: test_sampling

      call SsTC_sampling_task_constructor(task=test_sampling, name="lp_sampling_test", &
                                          l_calculator=l_calc_test, &
                                          nkpts=size(vecs(:, 1)), kpts=vecs, &
                                          N_int_ind=1, int_ind_range=(/10/), &
                                          N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                          ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_sampling_task(task=test_sampling, &
                                     system=SsTC_toy_model)

      call SsTC_print_sampling(task=test_sampling, &
                               system=SsTC_toy_model)

      res = real(test_sampling%predefined_sampled_data(7, 1, 2), dp) !Retrieve Re[exp(2i*(-0.25)*7)]

      if (abs(res + 0.93645668729079634) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                &of a local calculator and the reference for predefined sampling.")
        return
      end if

    end block local_calculator

  end subroutine test_P

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

end module suite_sampler
