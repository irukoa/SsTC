module calculator_test

  use utility
  use data_structures
  use local_k_quantities

  implicit none

  public :: calculator_test_C1M3

  contains

  function calculator_test_C1M3(task, system, k) result(u)
    type(BZ_integral_task), intent(in) :: task
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

end module calculator_test