module integrator

  USE OMP_LIB

  use utility
  use system
  use extrapolation_integration

  public :: sample_and_integrate_in_BZ_CHLM
  public :: sample_and_integrate_in_BZ_CLHM

  contains

  !TODO: Check why it fails on parallel. Make a progress bar? Check bounds of the BZ integration.
  !This subroutine is computationally heavy (CH) and light memorywise (LM).
  subroutine sample_and_integrate_in_BZ_CHLM(task, system, external_variable_data, samples, calculator_CHLM)

    USE OMP_LIB

    type(BZ_integrated_data), intent(inout) :: task
    type(sys),                intent(in)    :: system
    type(external_vars),      intent(in)    :: external_variable_data(:)
    integer,                  intent(in)    :: samples(3)

    interface
      function calculator_CHLM(task, system, external_variable_data, i_arr, r_arr, k) result(u)
        import :: BZ_integrated_data, sys, external_vars, dp
        type(BZ_integrated_data), intent(in) :: task
        type(sys),                intent(in) :: system
        type(external_vars),      intent(in) :: external_variable_data(:)
        integer,                  intent(in) :: i_arr(:), r_arr(:)
        real(kind=dp),            intent(in) :: k(3)

        complex(kind=dp)                     :: u
      end function calculator_CHLM
    end interface

    complex(kind=dp), allocatable :: data_k(:, :, :), sdata_k(:)

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r
    integer :: i_arr(size(task%integer_indices)), r_arr(size(task%continuous_indices))

    integer :: TID 

    write(unit=112, fmt=*) "Starting BZ sampling and integration subroutine CHML."
    write(unit=112, fmt=*) "Integrating task: "//trim(task%task)//" in the BZ for the system "//trim(system%name)//"."
    write(unit=112, fmt=*) "The required memory for the integration process is approximately,"
    write(unit=112, fmt=*) 16.0_dp*real(product(samples),dp)/1024.0_dp**2, "MB."
    write(unit=112, fmt=*) "Some computers limit the maximum memory an array can allocate."
    write(unit=112, fmt=*) "If this is your case and SIGSEGV triggers try using the next command before executing tb.x:"
    write(unit=112, fmt=*) "ulimit -s unlimited"

    allocate(data_k(samples(1), samples(2), samples(3)), sdata_k(product(samples)))

    do i = 1, product(task%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(task, i) !Pass to array layout.

      if ((task%particular_integer_component.ne.0).and.(i/=task%particular_integer_component)) cycle !Check if only some component is wanted.

      do r = 1, product(task%continuous_indices) !For each continuous index.

        r_arr = continuous_memory_element_to_array_element(task, r) !Pass to array layout.

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)
        TID = OMP_GET_THREAD_NUM()
        IF ((TID .EQ. 0).AND.(r .eq. 1)) write(unit=112, fmt=*) "Component ", i_arr, ". Running on ", OMP_GET_NUM_THREADS(), " threads."
        !$OMP DO
          do ik1=1, samples(1)
            k(1) = real(ik1-1,dp)/real(samples(1)-1,dp)
            do ik2 = 1, samples(2)
              k(2) = real(ik2-1,dp)/real(samples(2)-1,dp)
              do ik3 = 1, samples(3)
                k(3) = real(ik3-1,dp)/real(samples(3)-1,dp)
                
                data_k(ik1, ik2, ik3) = calculator_CHLM(task, system, external_variable_data, i_arr, r_arr, k)
                
              enddo
            enddo
          enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !Pass data array to memory layout.
        call shrink_array(data_k, sdata_k, info)
        !Integrate, if possible extrapolation method.
        call integral_extrapolation(sdata_k, samples, (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/), task%res(i, r), info)

      enddo

    enddo

    deallocate(data_k, sdata_k)

  end subroutine sample_and_integrate_in_BZ_CHLM

  !This subroutine is computationally light (CL) and heavy memorywise (HM).
  subroutine sample_and_integrate_in_BZ_CLHM(task, system, external_variable_data, samples, calculator_CLHM)

    type(BZ_integrated_data), intent(inout) :: task
    type(sys),                intent(in)    :: system
    type(external_vars),      intent(in)    :: external_variable_data(:)
    integer,                  intent(in)    :: samples(3)

    interface
      function calculator_CLHM(task, system, external_variable_data, i_arr, k) result(u)
        import :: BZ_integrated_data, sys, external_vars, dp
        type(BZ_integrated_data), intent(in) :: task
        type(sys),                intent(in) :: system
        type(external_vars),      intent(in) :: external_variable_data(:)
        integer,                  intent(in) :: i_arr(:)
        real(kind=dp),            intent(in) :: k(3)

        complex(kind=dp)                     :: u(product(task%continuous_indices))
      end function calculator_CLHM
    end interface

    complex(kind=dp) :: data_k(samples(1), samples(2), samples(3), product(task%continuous_indices)), sdata_k(product(samples), product(task%continuous_indices))

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r
    integer :: i_arr(size(task%integer_indices))

    do i = 1, product(task%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(task, i) !Pass to array layout.

      if ((task%particular_integer_component.ne.0).and.(i/=task%particular_integer_component)) cycle !Check if only some component is wanted.

      !!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ik1, ik2, ik3)!For each k-point.
        do ik1=1, samples(1)
          k(1) = real(ik1-1,dp)/real(samples(1)-1,dp)
          do ik2 = 1, samples(2)
            k(2) = real(ik2-1,dp)/real(samples(2)-1,dp)
            do ik3 = 1, samples(3)
              k(3) = real(ik3-1,dp)/real(samples(3)-1,dp)
              
              data_k(ik1, ik2, ik3, :) = calculator_CLHM(task, system, external_variable_data, i_arr, k)
              
            enddo
          enddo
        enddo
      !!$OMP END PARALLEL DO

      do r = 1, product(task%continuous_indices) !For each continuous index.
      !Pass data array to memory layout.
      call shrink_array(data_k(:, :, :, r), sdata_k(:, r), info)
      !Integrate, if possible extrapolation method.
      call integral_extrapolation(sdata_k(:, r), samples, (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/), task%res(i, r), info)
      enddo

    enddo
  end subroutine sample_and_integrate_in_BZ_CLHM

  !Sub to integrate calculators which return a complex array with integer and ontinuous indices.
  subroutine sample_and_integrate_in_BZ_C1M3(task, system, external_variable_data, samples, calculator_C1M3)

    type(BZ_integrated_data), intent(inout) :: task
    type(sys),                intent(in)    :: system
    type(external_vars),      intent(in)    :: external_variable_data(:)
    integer,                  intent(in)    :: samples(3)

    interface
      function calculator_C1M3(task, system, external_variable_data, k) result(u)
        import :: BZ_integrated_data, sys, external_vars, dp
        type(BZ_integrated_data), intent(in) :: task
        type(sys),                intent(in) :: system
        type(external_vars),      intent(in) :: external_variable_data(:)
        real(kind=dp),            intent(in) :: k(3)

        complex(kind=dp)                     :: u(product(task%integer_indices), product(task%continuous_indices))
      end function calculator_C1M3
    end interface

    complex(kind=dp), allocatable :: data_k(: , :, :, :, :), &
                                     sdata_k(:, :, :), temp_res(:, :)

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r

    integer :: TID

    if (task%method == "extrapolation") then 

      write(unit=112, fmt="(A)") "Starting BZ sampling and integration subroutine C1M3 with extrapolation method."
      write(unit=112, fmt="(A)") "Integrating task: "//trim(task%task)//" in the BZ for the system "//trim(system%name)//"."
      write(unit=112, fmt="(A)") "The required memory for the integration process is approximately,"
      write(unit=112, fmt="(F15.3, A)") 16.0_dp*real(product(samples)*product(task%integer_indices)*product(task%continuous_indices),dp)/1024.0_dp**2, "MB."
      write(unit=112, fmt="(A)") "Some computers limit the maximum memory an array can allocate."
      write(unit=112, fmt="(A)") "If this is your case and SIGSEGV triggers try using the next command before executing tb.x:"
      write(unit=112, fmt="(A)") "ulimit -s unlimited"

      allocate(data_k(samples(1), samples(2), samples(3), product(task%integer_indices), product(task%continuous_indices)), &
              sdata_k(product(samples), product(task%integer_indices), product(task%continuous_indices)))

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)

        TID = OMP_GET_THREAD_NUM()
        IF (TID .EQ. 0) THEN
          write(unit=112, fmt="(A, I6, A)") "Running on ", OMP_GET_NUM_THREADS(), " threads."
        ENDIF

        !$OMP DO
        do ik1=1, samples(1)

          k(1) = real(ik1-1,dp)/real(samples(1)-1,dp)
          do ik2 = 1, samples(2)
            k(2) = real(ik2-1,dp)/real(samples(2)-1,dp)
            do ik3 = 1, samples(3)
              k(3) = real(ik3-1,dp)/real(samples(3)-1,dp)
              
              data_k(ik1, ik2, ik3, :, :) = calculator_C1M3(task, system, external_variable_data, k)
              
            enddo
          enddo

        enddo
      !$OMP END DO
      !$OMP END PARALLEL

      write(unit=112, fmt="(A)") "Sampling done. Starting integration with extrapolation method."

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
        !Pass data array to memory layout.
        call shrink_array(data_k(:, :, :, i, r), sdata_k(:, i, r), info)
        !Integrate, if possible extrapolation method.
        call integral_extrapolation(sdata_k(:, i, r), samples, (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/), task%res(i, r), info)
        enddo
      enddo

      if (info==1) then
        write(unit=112, fmt="(A)") "Integral done. Extrapolation sucessfull."
      else
        write(unit=112, fmt="(A)") "Integral done. Extrapolation fail and return rectangle approximation."
      endif

      deallocate(data_k, sdata_k)

    elseif (task%method == "rectangle") then

      write(unit=112, fmt="(A)") "Starting BZ sampling and integration subroutine C1M3 with rectangle method."
      write(unit=112, fmt="(A)") "Integrating task: "//trim(task%task)//" in the BZ for the system "//trim(system%name)//"."
      write(unit=112, fmt="(A)") "The required memory for the integration process is approximately,"
      write(unit=112, fmt="(F15.3, A)") 16.0_dp*real(product(task%integer_indices)*product(task%continuous_indices),dp)/1024.0_dp**2, "MB."

      allocate(temp_res(product(task%integer_indices), product(task%continuous_indices)))

      !$OMP PARALLEL

      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        write(unit=112, fmt="(A, I6, A)") "Running on ", OMP_GET_NUM_THREADS(), " threads."
      ENDIF

      !$OMP DO REDUCTION(+:temp_res)
      do ik1=1, samples(1)

        k(1) = real(ik1-1,dp)/real(samples(1)-1,dp)
        do ik2 = 1, samples(2)
          k(2) = real(ik2-1,dp)/real(samples(2)-1,dp)
          do ik3 = 1, samples(3)
            k(3) = real(ik3-1,dp)/real(samples(3)-1,dp)
            
            temp_res = temp_res + calculator_C1M3(task, system, external_variable_data, k)
            
          enddo
        enddo

      enddo
    !$OMP END DO
    !$OMP END PARALLEL

    task%res = temp_res/product(samples)

    deallocate(temp_res)

    write(unit=112, fmt="(A)") "Sampling done. Starting integration with rectangle method."
    endif

  end subroutine sample_and_integrate_in_BZ_C1M3

end module integrator