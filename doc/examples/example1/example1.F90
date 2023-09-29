program example01

  USE OMP_LIB
  USE MPI_F08

  use SsTC

  implicit none

  integer, parameter :: dp = 8

  integer :: ierror

  type(SsTC_sys) :: dummy

  call MPI_INIT(ierror)

  call SsTC_init()

  dummy = SsTC_sys_constructor("dummy", "./", efermi = 0.0_dp)

  block

    type(SsTC_BZ_integral_task) :: test_integral

    call SsTC_BZ_integral_task_constructor(task = test_integral, name = "example01-test", &
                                           g_calculator = test_calculator, &
                                           method = "rectangle", samples = (/33, 1, 1/), &
                                           N_int_ind = 2, int_ind_range = (/3, 3/), &
                                           N_ext_vars = 1, ext_vars_start = (/1.0_dp/), &
                                           ext_vars_end = (/10.0_dp/), ext_vars_steps = (/10/))

    call SsTC_sample_and_integrate_BZ_integral_task(task = test_integral, &
                                                    system = dummy)

    call SsTC_print_BZ_integral_task(task = test_integral, &
                                     system = dummy)

  end block

  call MPI_FINALIZE(ierror)

  contains

  function test_calculator(task, system, k, error) result(u)

    class(SsTC_global_k_data), intent(in) :: task
    type(SsTC_sys), intent(in)            :: system
    real(kind=dp), intent(in)             :: k(3)
    logical, intent(inout)                :: error

    complex(kind=dp) :: u(product(task%integer_indices), product(task%continuous_indices))

    integer :: i, r, & !Indices for memory layout.
               i_arr(2), r_arr(1) !Indices for array layout.

    u = cmplx(0.0, 0.0, dp)

    do i = 1, product(task%integer_indices)
      i_arr = SsTC_integer_memory_element_to_array_element(task, i)
      !i_arr(1) now stores the first integer index.
      !i_arr(2) now stores the second integer index.
      do r = 1, product(task%continuous_indices)
        r_arr = SsTC_continuous_memory_element_to_array_element(task, r)
        !r_arr(1) now stores the iteration of the continuous variable.
        u(i, r) = real(i_arr(1) + i_arr(2), dp) * &
                  cmplx(exp(                      &
                              k(1)*               &
                              task%ext_var_data(1)%data(r_arr(1)) & !This references the value r_arr(1) of the continuous variable "1".
                        ))
      enddo
    enddo

  end function test_calculator

end program example01