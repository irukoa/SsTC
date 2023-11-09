module suite_kpath

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_kpath

contains

  subroutine collect_suite_kpath(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test Kpath", test_k) &
                ]

  end subroutine collect_suite_kpath

  subroutine test_k(error)
    type(error_type), allocatable, intent(out) :: error
    type(SsTC_sys) :: SsTC_toy_model

    real(kind=dp) :: vecs(2, 3)

    real(kind=dp) :: res
    logical :: cond = .true.

    vecs(1, :) = (/-0.5_dp, 0.0_dp, 0.0_dp/)
    vecs(2, :) = (/0.5_dp, 0.0_dp, 0.0_dp/)

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    global_calculator:block

      !g_calc_test implements the integrand of the sinc fucntion.
      type(SsTC_kpath_task) :: test_kpath

      call SsTC_kpath_constructor(task=test_kpath, name="g_kpath_test", &
                                  g_calculator=g_calc_test, &
                                  Nvec=2, vec_coord=vecs, nkpts=(/5/), &
                                  N_int_ind=1, int_ind_range=(/1/), &
                                  N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                  ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_kpath_sampler(task=test_kpath, &
                              system=SsTC_toy_model)

      call SsTC_print_kpath(task=test_kpath, &
                            system=SsTC_toy_model)

      res = real(test_kpath%kpath_data(1, 16, 2), dp) !Retrieve Re[exp(2i*(-0.25)*0.30303030E+001)]

      if (abs(res - 5.5616100165806738E-002) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                 &of a global calculator and the reference for kpath.")
        return
      end if

    end block global_calculator

    cond = .true.

    local_calculator:block

      type(SsTC_kpath_task) :: test_kpath

      call SsTC_kpath_constructor(task=test_kpath, name="l_kpath_test", &
                                  l_calculator=l_calc_test, &
                                  Nvec=2, vec_coord=vecs, nkpts=(/5/), &
                                  N_int_ind=1, int_ind_range=(/10/), &
                                  N_ext_vars=1, ext_vars_start=(/0.0_dp/), &
                                  ext_vars_end=(/20.0_dp/), ext_vars_steps=(/100/))

      call SsTC_kpath_sampler(task=test_kpath, &
                              system=SsTC_toy_model)

      call SsTC_print_kpath(task=test_kpath, &
                            system=SsTC_toy_model)

      res = real(test_kpath%kpath_data(7, 1, 2), dp) !Retrieve Re[exp(2i*(-0.25)*7)]

      if (abs(res + 0.93645668729079634) > 1.0E-6_dp) cond = .false.

      if (.not. cond) then
        call test_failed(error, "Mismatch between the results of the sampling &
                                &of a local calculator and the reference for kpath.")
        return
      end if

    end block local_calculator

  end subroutine test_k

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

end module suite_kpath
