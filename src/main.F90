program floquet_tight_binding

  USE OMP_LIB
  use, intrinsic :: iso_c_binding 
  use utility
  use calculator
  use integrator
  include 'fftw3.f03'

  type(BZ_integrated_data) :: test_result(2)

  !Declare that task "test" aims to calculate a quantitiy which has 2 integer indices of size 3 and a contonous index of size 1.
  test_result(1)%task = "test"
  allocate(test_result(1)%integer_indices(2))
  test_result(1)%integer_indices(1) = 3
  test_result(1)%integer_indices(2) = 3
  allocate(test_result(1)%continuous_indices(1))
  test_result(1)%continuous_indices(1) = 1
  allocate(test_result(1)%res(product(test_result(1)%integer_indices),&
  product(test_result(1)%continuous_indices)))

  !Declare that task "test" aims to calculate a quantitiy which has 2 integer indices of size 3 and a contonous index of size 1,
  !but we only weant to calculate the 2nd integer component.
  test_result(2)%task = "test"
  allocate(test_result(2)%integer_indices(2))
  test_result(2)%integer_indices(1) = 3
  test_result(2)%integer_indices(2) = 3
  allocate(test_result(2)%continuous_indices(1))
  test_result(2)%continuous_indices(1) = 1
  allocate(test_result(2)%res(product(test_result(2)%integer_indices),&
  product(test_result(2)%continuous_indices)))
  test_result(2)%particular_integer_component = integer_array_element_to_memory_element(test_result(2), (/3, 2/))

  call omp_set_nested(.true.)
  !call OMP_SET_MAX_ACTIVE_LEVELS(2)

  call sample_and_integrate_in_BZ(test_result(1), (/2049, 1, 1/))

  print*, "Task:", test_result(1)%task
  print*, test_result(1)%res

  call sample_and_integrate_in_BZ(test_result(2), (/2049, 1, 1/))

  print*, "Task:", test_result(2)%task
  print*, test_result(2)%res

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