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

  type(local_k_data), allocatable :: H(:), POS(:)!TEST

  integer :: i_mem
  integer, allocatable :: i_arr(:)

  open(unit=112, action="write", file="exec.out")

  !call OMP_SET_NUM_THREADS(4)
  call OMP_SET_MAX_ACTIVE_LEVELS(2)

  a = sys_constructor("HM", "./")

  call get_hamiltonian(system = a, k = (/0.0_dp, 0.0_dp, 0.0_dp/), H = H, Nder_i = 1, only_i = .false.)!TEST
  call get_position(system = a, k = (/0.0_dp, 0.0_dp, 0.0_dp/), A = POS, Nder_i = 1, only_i = .false.)!TEST
  print*, H(1)%k_data!Hamiltonian 0th derivative.
  print*, ""
  allocate(i_arr(product(H(2)%integer_indices)))
  do i_mem = 1, product(H(2)%integer_indices)
    write(*, fmt="(3I2, 2F8.5)"), integer_memory_element_to_array_element(H(2), i_mem), &
    real(H(2)%k_data(i_mem), dp), aimag(H(2)%k_data(i_mem))!Hamiltonian 1st Derivative.
  enddo
  deallocate(i_arr)
  print*, ""
  allocate(i_arr(product(POS(2)%integer_indices)))
  do i_mem = 1, product(POS(2)%integer_indices)
    write(*, fmt="(4I2, 2F8.5)"), integer_memory_element_to_array_element(POS(2), i_mem), &
    real(POS(2)%k_data(i_mem), dp), aimag(POS(2)%k_data(i_mem))!Pos. op. 1st Derivative.
  enddo

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