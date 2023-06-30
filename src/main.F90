program floquet_tight_binding

  USE OMP_LIB

  use, intrinsic :: iso_c_binding 

  use utility
  use system
  use integrator

  include 'fftw3.f03'

  type(sys) :: a

  type(external_vars) :: omega(1)
  type(BZ_integrated_data) :: test, test2

  open(unit=112, action="write", file="exec.out")

  call omp_set_nested(.true.)
  !call OMP_SET_MAX_ACTIVE_LEVELS(2)

  !EXAMPLE OF USAGE.
  omega(1) = external_variable_constructor(start = 0.0_dp,  &
                                           end   = 10.0_dp, &
                                           steps = 11)

  test = task_constructor(name      = "benchmark", &
                          nint      = 2,           &
                          int_range = (/3, 3/),    &
                          ext_vars  = omega,       &
                          method = "extrapolation" )

  a%name = "CHLM"

  call sample_and_integrate_in_BZ_CHLM(task = test,                    &
                                       system = a,                     &
                                       external_variable_data = omega, &
                                       samples = (/65, 9, 9/),         & !Carefull, this corresponds to a 16*product(samples) byte array, if it surpasses 4MB it will yield SIGSEGV. Use "ulimit -s unlimited"
                                       calculator_CHLM = calculator_test_CHLM)

  call print_task_result(task = test, &
                         system = a, &
                         external_variables = omega)

  a%name = "CLHM"

  call sample_and_integrate_in_BZ_CLHM(task = test,                    &
                                       system = a,                     &
                                       external_variable_data = omega, &
                                       samples = (/65, 9, 9/),         & !Carefull, this corresponds to a 16*product(samples)*product(task%continuous_indices) byte array, if it surpasses 4MB it will yield SIGSEGV. Use "ulimit -s unlimited"
                                       calculator_CLHM = calculator_test_CLHM)

  call print_task_result(task = test, &
                         system = a, &
                         external_variables = omega)

  a%name = "C1M3"

  call sample_and_integrate_in_BZ_C1M3(task = test,                    &
                                      system = a,                     &
                                      external_variable_data = omega, &
                                      samples = (/129, 9, 9/),         & !Carefull, this corresponds to a 16*product(samples)*product(task%continuous_indices) byte array, if it surpasses 4MB it will yield SIGSEGV. Use "ulimit -s unlimited"
                                      calculator_C1M3 = calculator_test_C1M3)

  call print_task_result(task = test, &
                        system = a, &
                        external_variables = omega)

  test2 = task_constructor(name      = "benchmarkrec", &
                           nint      = 2,           &
                           int_range = (/3, 3/),    &
                           ext_vars  = omega,       &
                           method = "rectangle" )

  call sample_and_integrate_in_BZ_C1M3(task = test2,                    &
                                       system = a,                     &
                                       external_variable_data = omega, &
                                       samples = (/200000, 9, 9/),         & !Carefull, this corresponds to a 16*product(samples)*product(task%continuous_indices) byte array, if it surpasses 4MB it will yield SIGSEGV. Use "ulimit -s unlimited"
                                       calculator_C1M3 = calculator_test_C1M3)

call print_task_result(task = test2, &
                       system = a, &
                       external_variables = omega)

  close(unit=112)

  contains

  function calculator_test_CLHM(task, system, external_variable_data, i_arr, k) result(u)

    type(BZ_integrated_data), intent(in) :: task
    type(sys),                intent(in) :: system
    type(external_vars),      intent(in) :: external_variable_data(:)
    integer,                  intent(in) :: i_arr(:)
    real(kind=dp),            intent(in) :: k(3)

    complex(kind=dp)                     :: u(product(task%continuous_indices))
    integer                              :: r, r_arr(size(task%continuous_indices))

    u = i_arr(1)*(k(1)**2)*exp(sin(10*k(1)))
    do r = 1, product(task%continuous_indices) !For each continuous index.
      r_arr = continuous_memory_element_to_array_element(task, r) !Pass to array layout.
      u(r) = u(r)*external_variable_data(1)%data(r_arr(1))
    enddo
  end function calculator_test_CLHM

  function calculator_test_C1M3(task, system, external_variable_data, k) result(u)
    type(BZ_integrated_data), intent(in) :: task
    type(sys),                intent(in) :: system
    type(external_vars),      intent(in) :: external_variable_data(:)
    real(kind=dp),            intent(in) :: k(3)

    complex(kind=dp)                     :: u(product(task%integer_indices), product(task%continuous_indices))
    integer                              :: r, r_arr(size(task%continuous_indices)), i, i_arr(size(task%integer_indices))

    u = (k(1)**2)*exp(sin(10*k(1)))
    do i = 1, product(task%integer_indices)
      i_arr = integer_memory_element_to_array_element(task, i)
      do r = 1, product(task%continuous_indices) !For each continuous index.
        r_arr = continuous_memory_element_to_array_element(task, r) !Pass to array layout.
        u(i, r) = u(i, r)*external_variable_data(1)%data(r_arr(1))*i_arr(1)
      enddo
    enddo

  end function calculator_test_C1M3

  
  function calculator_test_CHLM(task, system, external_variable_data, i_arr, r_arr, k) result(u)

    type(BZ_integrated_data), intent(in) :: task
    type(sys),                intent(in) :: system
    type(external_vars),      intent(in) :: external_variable_data(:)
    integer,                  intent(in) :: i_arr(:), r_arr(:)
    real(kind=dp),            intent(in) :: k(3)

    complex(kind=dp)                     :: u

    u = external_variable_data(1)%data(r_arr(1))*i_arr(1)*(k(1)**2)*exp(sin(10*k(1)))

  end function calculator_test_CHLM

  subroutine init_model

    real(kind=dp) :: lattice_const

    real(kind=dp) :: offset, tunnelling

    real(kind=dp), dimension(2) :: dipole

    real(kind=dp), dimension(2, 2) :: brav, recip 

    lattice_const = 1.0_dp !in A.

    offset     = 0.0_dp !In eV.
    tunnelling = 1.0_dp

    dipole     = lattice_const*(/1.0_dp, 0.0_dp/)/sqrt(3.0_dp) !In A.

    brav(1, :) = 0.5_dp*(/sqrt(3.0_dp),  1.0_dp/)!Relative to lattice_const.
    brav(2, :) = 0.5_dp*(/-sqrt(3.0_dp), 1.0_dp/)

    recip(1, :) = 2*pi*(/1.0_dp,  sqrt(3.0_dp)/)/sqrt(3.0_dp)!Relative to 1/lattice_const.
    recip(2, :) = 2*pi*(/-1.0_dp, sqrt(3.0_dp)/)/sqrt(3.0_dp)

    print*, dot_product(brav(1, :), recip(2, :))

  end subroutine init_model

end program floquet_tight_binding