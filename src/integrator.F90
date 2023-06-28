module integrator

  USE OMP_LIB
  use utility
  use extrapolation_integration
  use calculator

  public :: sample_and_integrate_in_BZ
  public :: integer_array_element_to_memory_element
  public :: integer_memory_element_to_array_element
  public :: continuous_array_element_to_memory_element
  public :: continuous_memory_element_to_array_element

  contains

  subroutine sample_and_integrate_in_BZ(task, samples)

    type(BZ_integrated_data), intent(inout) :: task
    integer,                  intent(in)    :: samples(:)

    complex(kind=dp), allocatable :: data_k(:, :, :), sdata_k(:)

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r
    integer :: i_arr(size(task%integer_indices)), r_arr(size(task%continuous_indices))

    do i = 1, product(task%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(task, i) !Pass to array layout.

      if ((task%particular_integer_component.ne.0).and.(i/=task%particular_integer_component)) cycle !Check if only some component is wanted.

      do r = 1, product(task%continuous_indices) !For each continuous index.

        r_arr = continuous_memory_element_to_array_element(task, r) !Pass to array layout.

        allocate(data_k(samples(1), samples(2), samples(3)), sdata_k(product(samples)))

        !!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ik1, ik2, ik3)!For each k-point.
          do ik1=1, samples(1)
            k(1) = real(ik1-1,dp)/real(samples(1)-1,dp)
            do ik2 = 1, samples(2)
              k(2) = real(ik2-1,dp)/real(samples(2)-1,dp)
              do ik3 = 1, samples(3)
                k(3) = real(ik3-1,dp)/real(samples(3)-1,dp)
      
                !TODO: Check why it fails on parallel. Make a progress bar?
      
                call calculator_dict(task, i_arr, r_arr, k, data_k(ik1, ik2, ik3))
              
              enddo
            enddo
          enddo
        !!$OMP END PARALLEL DO

        call shrink_array(data_k, sdata_k, info) !Pass data array to memory layout.
        call integral_extrapolation(sdata_k, samples, (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/), task%res(i, r), info) !Integrate, if possible extrapolation method.
        deallocate(data_k, sdata_k)

      enddo

    enddo
  end subroutine sample_and_integrate_in_BZ

  function integer_array_element_to_memory_element(task, i_arr) result (i_mem)
    !Get integer indices from array layout to memory layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: i_arr(size(task%integer_indices))
    integer :: i_mem

    integer :: counter

    i_mem = 1
    do counter = 0, size(task%integer_indices) - 1
      i_mem = (i_mem-1)*task%integer_indices(counter+1) + i_arr(counter+1)
    enddo

  end function integer_array_element_to_memory_element

  function integer_memory_element_to_array_element(task, i_mem) result (i_arr)
    !Get integer indices from memory layout to array layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: i_mem
    integer :: i_arr(size(task%integer_indices))

    integer :: counter, reduction

    reduction = i_mem
    do counter = size(task%integer_indices), 1, -1
      if (counter == 1) then
        i_arr(counter) = reduction
      else
        i_arr(counter) = modulo(reduction, task%integer_indices(counter))
        if (i_arr(counter)==0) i_arr(counter) = task%integer_indices(counter)
        reduction = int((reduction-i_arr(counter))/task%integer_indices(counter)) + 1
      endif
    enddo

  end function integer_memory_element_to_array_element

  function continuous_array_element_to_memory_element(task, r_arr) result (r_mem)
    !Get continuous indices from array layout to memory layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: r_arr(size(task%continuous_indices))
    integer :: r_mem

    integer :: counter

    r_mem = 1
    do counter = 0, size(task%continuous_indices) - 1
      r_mem = (r_mem-1)*task%continuous_indices(counter+1) + r_arr(counter+1)
    enddo

  end function continuous_array_element_to_memory_element

  function continuous_memory_element_to_array_element(task, r_mem) result (r_arr)
    !Get continuous indices from memory layout to array layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: r_mem
    integer :: r_arr(size(task%continuous_indices))

    integer :: counter, reduction

    reduction = r_mem
    do counter = size(task%continuous_indices), 1, -1
      if (counter == 1) then
        r_arr(counter) = reduction
      else
        r_arr(counter) = modulo(reduction, task%continuous_indices(counter))
        if (r_arr(counter)==0) r_arr(counter) = task%continuous_indices(counter)
        reduction = int((reduction-r_arr(counter))/task%continuous_indices(counter)) + 1
      endif
    enddo

  end function continuous_memory_element_to_array_element

end module integrator