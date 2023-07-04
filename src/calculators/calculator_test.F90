module calculator_test

  use utility
  use data_structures
  use local_k_quantities

  implicit none

  public :: calculator_test_C1M3
  public :: bands

  contains

  function calculator_test_C1M3(task, system, k) result(u)
    type(global_k_data), intent(in) :: task
    type(sys),              intent(in) :: system
    real(kind=dp),          intent(in) :: k(3)

    complex(kind=dp)                   :: u(product(task%integer_indices), product(task%continuous_indices))
    integer                            :: r, r_arr(size(task%continuous_indices)), i, i_arr(size(task%integer_indices))

    real(kind=dp)                      :: part

    u = 0.0_dp
    part = ((k(1)+0.5_dp)**2)*exp(sin(10*(k(1)+0.5_dp)))
    do i = 1, product(task%integer_indices)
      if ((task%particular_integer_component.ne.0).and.(i.ne.task%particular_integer_component)) cycle
      i_arr = integer_memory_element_to_array_element(task, i)
      do r = 1, product(task%continuous_indices) !For each continuous index.
        r_arr = continuous_memory_element_to_array_element(task, r) !Pass to array layout.
        u(i, r) = part*task%ext_var_data(1)%data(r_arr(1))*i_arr(1)
      enddo
    enddo

  end function calculator_test_C1M3

  function bands(k_data, system, k) result(u)
    class(local_k_data), intent(in) :: k_data
    type(sys),          intent(in)  :: system
    real(kind=dp),      intent(in)  :: k(3)

    complex(kind=dp)                :: u(product(k_data%integer_indices))

    type(local_k_data), allocatable :: H(:)
    complex(kind=dp) :: hamiltonian(k_data%integer_indices(1), k_data%integer_indices(1)),&
                        rot(k_data%integer_indices(1), k_data%integer_indices(1))
    real(kind=dp) :: eig(k_data%integer_indices(1))
    integer, allocatable :: i_arr(:)
    integer :: i_mem

    logical :: error
    character(len=120) :: errormsg

    call get_hamiltonian(system, k, H)

    do i_mem = 1, product(H(1)%integer_indices)
      i_arr = integer_memory_element_to_array_element(H(1), i_mem)
      hamiltonian(i_arr(1), i_arr(2)) = H(1)%k_data(i_mem)
    enddo
    deallocate(H, i_arr)
    call utility_diagonalize(hamiltonian, k_data%integer_indices(1), eig, rot, error, errormsg)
    u = eig
  end function bands

end module calculator_test