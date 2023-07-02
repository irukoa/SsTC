module integrator

  USE OMP_LIB

  use data_structures
  use extrapolation_integration

  integer, parameter, private :: dp = 8

  public :: sample_and_integrate_in_BZ

  contains

  !Sub to integrate calculators which return a complex array with integer and continuous indices.
  subroutine sample_and_integrate_in_BZ(task, system)

    type(BZ_integral_task), intent(inout) :: task
    type(sys),                intent(in)    :: system

    complex(kind=dp), allocatable :: data_k(: , :, :, :, :), &
                                     sdata_k(:, :, :), temp_res(:, :)

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r

    integer :: TID

    if (task%method == "extrapolation") then 

      write(unit=112, fmt="(A)") "Starting BZ sampling and integration subroutine C1M3 with extrapolation method."
      write(unit=112, fmt="(A)") "Integrating task: "//trim(task%name)//" in the BZ for the system "//trim(system%name)//"."
      write(unit=112, fmt="(A)") "The required memory for the integration process is approximately,"
      write(unit=112, fmt="(F15.3, A)") 16.0_dp*real(product(task%samples)*product(task%integer_indices)*product(task%continuous_indices),dp)/1024.0_dp**2, "MB."
      write(unit=112, fmt="(A)") "Some computers limit the maximum memory an array can allocate."
      write(unit=112, fmt="(A)") "If this is your case and SIGSEGV triggers try using the next command before executing tb.x:"
      write(unit=112, fmt="(A)") "ulimit -s unlimited"

      allocate(data_k(task%samples(1), task%samples(2), task%samples(3), product(task%integer_indices), product(task%continuous_indices)), &
              sdata_k(product(task%samples), product(task%integer_indices), product(task%continuous_indices)))

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)

        TID = OMP_GET_THREAD_NUM()
        IF (TID .EQ. 0) THEN
          write(unit=112, fmt="(A, I6, A)") "Running on ", OMP_GET_NUM_THREADS(), " threads."
        ENDIF

        !$OMP DO
        do ik1=1, task%samples(1)

          k(1) = real(ik1-1,dp)/real(task%samples(1)-1,dp)
          do ik2 = 1, task%samples(2)
            k(2) = real(ik2-1,dp)/real(task%samples(2)-1,dp)
            do ik3 = 1, task%samples(3)
              k(3) = real(ik3-1,dp)/real(task%samples(3)-1,dp)
              
              data_k(ik1, ik2, ik3, :, :) = task%calculator(task, system, k)
              
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
        call integral_extrapolation(sdata_k(:, i, r), task%samples, (/0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp/), task%result(i, r), info)
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
      write(unit=112, fmt="(A)") "Integrating task: "//trim(task%name)//" in the BZ for the system "//trim(system%name)//"."
      write(unit=112, fmt="(A)") "The required memory for the integration process is approximately,"
      write(unit=112, fmt="(F15.3, A)") 16.0_dp*real(product(task%integer_indices)*product(task%continuous_indices),dp)/1024.0_dp**2, "MB."

      allocate(temp_res(product(task%integer_indices), product(task%continuous_indices)))

      !$OMP PARALLEL

      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        write(unit=112, fmt="(A, I6, A)") "Running on ", OMP_GET_NUM_THREADS(), " threads."
      ENDIF

      !$OMP DO REDUCTION(+:temp_res)
      do ik1=1, task%samples(1)

        k(1) = real(ik1-1,dp)/real(task%samples(1)-1,dp)
        do ik2 = 1, task%samples(2)
          k(2) = real(ik2-1,dp)/real(task%samples(2)-1,dp)
          do ik3 = 1, task%samples(3)
            k(3) = real(ik3-1,dp)/real(task%samples(3)-1,dp)
            
            temp_res = temp_res + task%calculator(task, system, k)
            
          enddo
        enddo

      enddo
    !$OMP END DO
    !$OMP END PARALLEL

    task%result = temp_res/product(task%samples)

    deallocate(temp_res)

    write(unit=112, fmt="(A)") "Sampling done. Starting integration with rectangle method."
    endif

  end subroutine sample_and_integrate_in_BZ

end module integrator