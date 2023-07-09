program floquet_tight_binding

  USE OMP_LIB

  use, intrinsic :: iso_c_binding 

  use utility
  use data_structures
  use integrator

  use kpath
  use kslice
  use calculator_test
  use calculators_general
  use calculators_floquet

  use local_k_quantities!TEST

  include 'fftw3.f03'

  type(sys) :: a

  type(BZ_integral_task) :: test, test2, test3

  type(k_path_task) :: path

  real(kind=dp) :: kvecs(4, 3)
  integer :: i

  kvecs(1, :) = (/0.0_dp, 0.0_dp, 0.0_dp/)
  kvecs(2, :) = (/0.5_dp, 0.0_dp, 0.0_dp/)
  kvecs(3, :) = (/1.0_dp/3, 1.0_dp/3, 0.0_dp/)
  kvecs(4, :) = (/0.0_dp, 0.0_dp, 0.0_dp/)

  !kvecs(2, :) = (/0.00000,      0.50000,      0.00000/)!L
  !kvecs(3, :) = (/0.00000,      0.62500,      0.37500/)!K
  !kvecs(4, :) = (/0.00000,      0.50000,      0.50000/)!X

  open(unit=112, action="write", file="exec.out")
  open(unit=113, action="write", file="exec.err")

  !call OMP_SET_NUM_THREADS(1) !SERIAL.
  call OMP_SET_MAX_ACTIVE_LEVELS(2)

  a = sys_constructor("GaAs", "./")
  !a = sys_constructor("HM", "./")

  path = bands_kpath_task_constructor(system = a, &
                                      Nvec = 4, &
                                      vec_coord = kvecs, &
                                      nkpts = (/100, 100, 100/))

  call kpath_sampler(path, a)
  call print_kpath(path, a)

  !EXAMPLE OF USAGE.
  test = task_constructor(name           = "ext_ben", &
                          g_calculator   = calculator_test_C1M3, &
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
                           g_calculator   = calculator_test_C1M3, &
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

  path = quasienergy_kpath_task_constructor(system = a, &
                                            Nvec = 2, &
                                            vec_coord = kvecs(3:4, :), &
                                            nkpts = (/1/), &
                                            Nharm = 1, &
                                            axstart = (/1.0E4_dp/), axend = (/1.0E5_dp/), axsteps = (/1/), &
                                            pxstart = (/0.0_dp/), pxend = (/0.0_dp/), pxsteps = (/1/), &
                                            aystart = (/0.0_dp/), ayend = (/0.0_dp/), aysteps = (/1/), &
                                            pystart = (/0.0_dp/), pyend = (/0.0_dp/), pysteps = (/1/), &
                                            azstart = (/0.0_dp/), azend = (/0.0_dp/), azsteps = (/1/), &
                                            pzstart = (/0.0_dp/), pzend = (/0.0_dp/), pzsteps = (/1/), &
                                            omegastart = 0.5_dp, omegaend = 5.0_dp, omegasteps = 100, &
                                            t0start = 0.0_dp, t0end = 0.0_dp, t0steps = 1)

call kpath_sampler(path, a)
call print_kpath(path, a)

  close(unit=112)
  close(unit=113)

end program floquet_tight_binding