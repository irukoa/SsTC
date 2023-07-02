program floquet_tight_binding

  USE OMP_LIB

  use, intrinsic :: iso_c_binding 

  use utility
  use data_structures
  use integrator

  use calculator_test

  include 'fftw3.f03'

  type(sys) :: a

  type(BZ_integral_task) :: test, test2

  open(unit=112, action="write", file="exec.out")

  call omp_set_nested(.true.)
  !call OMP_SET_MAX_ACTIVE_LEVELS(2)

  a = sys_constructor("GaAs", "./")

  !EXAMPLE OF USAGE.
  test = task_constructor(name           = "ext_ben", &
                          calculator     = calculator_test_C1M3, &
                          N_int_ind      = 2, &
                          int_ind_range  = (/3, 3/), &
                          N_ext_vars     = 1, &
                          ext_vars_start = (/0.0_dp/), &
                          ext_vars_end   = (/10.0_dp/), &
                          ext_vars_steps = (/11/), &
                          method         = "extrapolation", & !Required memory: 16*product(samples)*product(int_ind_range)*product(ext_vars_steps)
                          samples        = (/65, 65, 65/))

  a%name = "C1M3"

  call sample_and_integrate_in_BZ(task = test, &
                                  system = a)                    

  call print_task_result(task = test, &
                         system = a)

  test2 = task_constructor(name           = "rec_ben", &
                           calculator     = calculator_test_C1M3, &
                           N_int_ind      = 2, &
                           int_ind_range  = (/3, 3/), &
                           N_ext_vars     = 1, &
                           ext_vars_start = (/0.0_dp/), &
                           ext_vars_end   = (/10.0_dp/), &
                           ext_vars_steps = (/11/), &
                           method         = "rectangle", & !Required memory: 16*product(int_ind_range)*product(ext_vars_steps)
                           samples        = (/200000, 9, 9/), &
                           part_int_comp  = (/2, 1/))

  call sample_and_integrate_in_BZ(task = test2, &
                                  system = a)

call print_task_result(task = test2, &
                       system = a)

  close(unit=112)

  contains

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