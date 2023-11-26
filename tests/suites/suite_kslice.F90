module suite_kslice

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_kslice

contains

  subroutine collect_suite_kslice(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test Kslice", test_ks) &
                ]

  end subroutine collect_suite_kslice

  subroutine test_ks(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    real(kind=dp) :: vecs(3, 3)

    real(kind=dp) :: res
    logical :: cond = .true.

    vecs(1, :) = (/-0.5_dp, -0.5_dp, 0.0_dp/)
    vecs(2, :) = (/1.0_dp, 0.0_dp, 0.0_dp/)
    vecs(3, :) = (/0.0_dp, 1.0_dp, 0.0_dp/)

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    global_calculator:block

      !g_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_kslice_task) :: test_kslice

      call SsTC_kslice_task_constructor(task=test_kslice, name="g_kslice_test", &
                                        g_calculator=g_calc_test, &
                                        corner=vecs(1, :), vector_a=vecs(2, :), vector_b=vecs(3, :), samples=(/101, 101/), &
                                        N_int_ind=1, int_ind_range=(/1/), &
                                        N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                        ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_kslice_task(task=test_kslice, &
                                   system=SsTC_toy_model)

      call SsTC_print_kslice(task=test_kslice, &
                             system=SsTC_toy_model)

      res = real(test_kslice%kslice_data(1, 16, 26, 5), dp) !Retrieve Re[exp(2i*(-0.25)*0.30303030E+001)]

      if (abs(res - 5.5616100165806738E-002) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                 &of a global calculator and the reference for kslice.")
        return
      end if

    end block global_calculator

    cond = .true.

    local_calculator:block

      type(SsTC_kslice_task) :: test_kslice

      call SsTC_kslice_task_constructor(task=test_kslice, name="l_kslice_test", &
                                        l_calculator=l_calc_test, &
                                        corner=vecs(1, :), vector_a=vecs(2, :), vector_b=vecs(3, :), samples=(/101, 101/), &
                                        N_int_ind=1, int_ind_range=(/10/), &
                                        N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                        ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_sample_kslice_task(task=test_kslice, &
                                   system=SsTC_toy_model)

      call SsTC_print_kslice(task=test_kslice, &
                             system=SsTC_toy_model)

      res = real(test_kslice%kslice_data(7, 1, 26, 5), dp) !Retrieve Re[exp(2i*(-0.25)*7)]

      if (abs(res + 0.93645668729079634) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                &of a local calculator and the reference for kslice.")
        return
      end if

    end block local_calculator

  end subroutine test_ks

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

end module suite_kslice
