program floquet_tight_binding

  USE OMP_LIB

  use, intrinsic :: iso_c_binding 

  use utility
  use data_structures
  use integrator

  use kpath
  use kslice
  use calculator_test

  use local_k_quantities!TEST

  include 'fftw3.f03'

  type(sys) :: a

  type(BZ_integral_task) :: test, test2

  type(k_path) :: path

  real(kind=dp) :: kvecs(4, 3)

  kvecs(1, :) = (/0.0_dp, 0.0_dp, 0.0_dp/)
  kvecs(2, :) = (/0.5_dp, 0.0_dp, 0.0_dp/)
  kvecs(3, :) = (/1.0_dp/3, 1.0_dp/3, 0.0_dp/)
  kvecs(4, :) = (/0.0_dp, 0.0_dp, 0.0_dp/)

  !kvecs(2, :) = (/0.00000,      0.50000,      0.00000/)!L
  !kvecs(3, :) = (/0.00000,      0.62500,      0.37500/)!K
  !kvecs(4, :) = (/0.00000,      0.50000,      0.50000/)!X

  open(unit=112, action="write", file="exec.out")

  !call OMP_SET_NUM_THREADS(4)
  call OMP_SET_MAX_ACTIVE_LEVELS(2)

  !a = sys_constructor("GaAs", "./")
  a = sys_constructor("HM", "./")

  path = kpath_constructor(name = "bands", &
                           l_calculator = bands, &
                           Nvec = 4, &
                           vec_coord = kvecs, &
                           nkpts = (/100, 100, 100/), &
                           N_int_ind = 1, &
                           int_ind_range = (/a%num_bands/), &
                           N_ext_vars = 1)

  call kpath_sampler(path, a)
  call print_kpath(path, a)

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
                           samples        = (/400000, 1, 1/), &
                           part_int_comp  = (/2, 1/))

  call sample_and_integrate_in_BZ(task = test2, &
                                  system = a)

call print_task_result(task = test2, &
                       system = a)

  close(unit=112)

end program floquet_tight_binding