module integrator

  USE OMP_LIB
  use utility
  use extrapolation_integration
  use calculator

  public :: sample_and_integrate_in_BZ

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

end module integrator