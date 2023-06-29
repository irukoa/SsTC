program floquet_tight_binding

  USE OMP_LIB

  use, intrinsic :: iso_c_binding 

  use utility
  use system
  use integrator

  include 'fftw3.f03'

  type(sys) :: a

  type(external_vars) :: omega(1)
  type(BZ_integrated_data) :: test

  call omp_set_nested(.true.)
  !call OMP_SET_MAX_ACTIVE_LEVELS(2)

  !EXAMPLE OF USAGE.
  omega(1) = external_variable_constructor(start = 0.0_dp,  &
                                           end   = 20.0_dp, &
                                           steps = 10)

  test = task_constructor(name      = "benchmark", &
                          nint      = 2,           &
                          int_range = (/3, 3/),    &
                          ext_vars  = omega)

  a%name = "test"

  call sample_and_integrate_in_BZ(task = test,                    &
                                  system = a,                     &
                                  external_variable_data = omega, &
                                  samples = (/2049, 1, 1/),       &
                                  calculator = calculator_test)

  call print_task_result(task = test, &
                         system = a, &
                         external_variables = omega)

  contains

  function calculator_test(task, system, external_variable_data, i_arr, r_arr, k) result(u)

    type(BZ_integrated_data), intent(in) :: task
    type(sys),                intent(in) :: system
    type(external_vars),      intent(in) :: external_variable_data(:)
    integer,                  intent(in) :: i_arr(:), r_arr(:)
    real(kind=dp),            intent(in) :: k(3)

    complex(kind=dp)                     :: u

    u = r_arr(1)*i_arr(1)*(k(1)**2)*exp(sin(10*k(1)))

  end function calculator_test

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