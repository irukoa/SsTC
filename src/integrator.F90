module integrator

  USE OMP_LIB

  use utility
  use system
  use extrapolation_integration

  public :: sample_and_integrate_in_BZ

  contains

  !TODO: Check why it fails on parallel. Make a progress bar? Check bounds of the BZ integration.
  subroutine sample_and_integrate_in_BZ(task, system, external_variable_data, samples, calculator)

    type(BZ_integrated_data), intent(inout) :: task
    type(sys),                intent(in)    :: system
    type(external_vars),      intent(in)    :: external_variable_data(:)
    integer,                  intent(in)    :: samples(3)

    interface
      function calculator(task, system, external_variable_data, i_arr, r_arr, k) result(u)
        import :: BZ_integrated_data, sys, external_vars, dp
        type(BZ_integrated_data), intent(in) :: task
        type(sys),                intent(in) :: system
        type(external_vars),      intent(in) :: external_variable_data(:)
        integer,                  intent(in) :: i_arr(:), r_arr(:)
        real(kind=dp),            intent(in) :: k(3)

        complex(kind=dp)                     :: u
      end function calculator
    end interface

    complex(kind=dp) :: data_k(samples(1), samples(2), samples(3)), sdata_k(product(samples))

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r
    integer :: i_arr(size(task%integer_indices)), r_arr(size(task%continuous_indices))

    do i = 1, product(task%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(task, i) !Pass to array layout.

      if ((task%particular_integer_component.ne.0).and.(i/=task%particular_integer_component)) cycle !Check if only some component is wanted.

      do r = 1, product(task%continuous_indices) !For each continuous index.

        r_arr = continuous_memory_element_to_array_element(task, r) !Pass to array layout.

        !!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ik1, ik2, ik3)!For each k-point.
          do ik1=1, samples(1)
            k(1) = real(ik1-1,dp)/real(samples(1)-1,dp)
            do ik2 = 1, samples(2)
              k(2) = real(ik2-1,dp)/real(samples(2)-1,dp)
              do ik3 = 1, samples(3)
                k(3) = real(ik3-1,dp)/real(samples(3)-1,dp)
                
                data_k(ik1, ik2, ik3) = calculator(task, system, external_variable_data, i_arr, r_arr, k)
                
              enddo
            enddo
          enddo
        !!$OMP END PARALLEL DO

        !Pass data array to memory layout.
        call shrink_array(data_k, sdata_k, info)
        !Integrate, if possible extrapolation method.
        call integral_extrapolation(sdata_k, samples, (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/), task%res(i, r), info)

      enddo

    enddo
  end subroutine sample_and_integrate_in_BZ

end module integrator