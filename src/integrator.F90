module integrator

  USE OMP_LIB

  use utility
  use data_structures
  use extrapolation_integration

  implicit none

  type, extends(global_k_data) :: BZ_integral_task
    !Integration method.
    character(len=120)                             :: method
    !Integration samples.
    integer                                        :: samples(3)
    !Result of the integration, contains the integer 
    !index and the continuous index in memory layout, respectively.
    complex(kind=dp), allocatable                  :: result(:, :)
  end type BZ_integral_task

  public :: task_constructor
  public :: sample_and_integrate_in_BZ
  public :: print_task_result

  contains

  function task_constructor(name, &
                            l_calculator, g_calculator, &
                            method, samples, &
                            N_int_ind, int_ind_range, &
                            N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                            part_int_comp) result(task)

    character(len=*) :: name

    procedure (local_calculator),  optional :: l_calculator
    procedure (global_calculator), optional :: g_calculator

    character(len=*), optional, intent(in) :: method
    integer,          optional, intent(in) :: samples(3)

    integer, optional, intent(in) :: N_int_ind
    integer, optional, intent(in) :: int_ind_range(:)


    integer,       optional, intent(in) :: N_ext_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(:), ext_vars_end(:)
    integer,       optional, intent(in) :: ext_vars_steps(:)

    integer, optional,  intent(in) :: part_int_comp(:)

    type(BZ_integral_task) :: task

    integer :: i

    !Set name.
    task%name = name

    !Set integer index data.
    if (((N_int_ind).ge.1).and.(present(int_ind_range))) then
      allocate(task%integer_indices(N_int_ind))
      do i = 1, N_int_ind
        task%integer_indices(i) = int_ind_range(i)
      enddo
    else
      allocate(task%integer_indices(1))
      task%integer_indices(1) = 1
    endif

    !Set external variable data.
    if (((N_ext_vars).ge.1).and.(present(ext_vars_start)).and.(present(ext_vars_end)).and.(present(ext_vars_steps))) then
      allocate(task%continuous_indices(N_ext_vars), task%ext_var_data(N_ext_vars))
      do i = 1, N_ext_vars
        task%continuous_indices(i) = ext_vars_steps(i)
        allocate(task%ext_var_data(i)%data(ext_vars_steps(i)))
        task%ext_var_data(i) = external_variable_constructor(ext_vars_start(i), ext_vars_end(i), ext_vars_steps(i))
      enddo
    else
      allocate(task%continuous_indices(1), task%ext_var_data(1))
      task%continuous_indices(1) = 1.0_dp
      allocate(task%ext_var_data(1)%data(1))
      task%ext_var_data(1)%data = (/1.0_dp/)
    endif

    if (present(g_calculator).and.(.not.present(l_calculator))) then
      !Set calculator pointer (function alias).
      task%global_calculator => g_calculator
      nullify(task%local_calculator)
    elseif (present(l_calculator).and.(.not.present(g_calculator))) then
      !Set calculator pointer (function alias).
      task%local_calculator => l_calculator
      nullify(task%global_calculator)
    endif

    !Set result.
    allocate(task%result(product(task%integer_indices), product(task%continuous_indices)))
    task%result = cmplx_0

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = integer_array_element_to_memory_element(task, part_int_comp)

    !Set integration method.
    if (present(method)) then
      if (method == "extrapolation") then
        task%method = "extrapolation"
        write(unit=112, fmt="(A)") "Warning: To employ the extrapolation method all the elements of the array 'samples' must be expressible as either 1 or 2^n + 1 for n = 0, 1, ..."
      elseif (method == "rectangle") then
        task%method = "rectangle"
      else
        print*, "Integration method not recognized. Setting up rectangle method"
        task%method = "rectangle"
      endif
    else
      task%method = "rectangle"
    endif

    !Set number of integration samples (discretization of BZ).
    if(present(samples)) then
      task%samples = samples
    else
      task%samples = (/10, 10, 10/)
    endif

  end function task_constructor

  !Sub to integrate calculators which return a complex array with integer and continuous indices.
  !The interface for the generic calculator function is given in data_structures module.
  subroutine sample_and_integrate_in_BZ(task, system)

    type(BZ_integral_task), intent(inout) :: task
    type(sys),              intent(in)    :: system

    complex(kind=dp), allocatable :: data_k(: , :, :, :, :), &
                                     sdata_k(:, :, :), &
                                     temp_res(:, :)

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, info, i, r

    integer :: TID

    if (task%method == "extrapolation") then !Extraplation case. Requires large RAM.

      write(unit=112, fmt = "(A)") "Starting BZ sampling and integration subroutine C1M3 with extrapolation method."
      write(unit=112, fmt = "(A)") "Integrating task: "//trim(task%name)//" in the BZ for the system "//trim(system%name)//"."
      write(unit=112, fmt = "(A)") "The required memory for the integration process is approximately,"
      write(unit=112, fmt = "(F15.3, A)") 16.0_dp*real(product(task%samples)*product(task%integer_indices)*product(task%continuous_indices),dp)/1024.0_dp**2, "MB."
      write(unit=112, fmt = "(A)") "Some computers limit the maximum memory an array can allocate."
      write(unit=112, fmt = "(A)") "If this is your case and SIGSEGV triggers try using the next command before executing tb.x:"
      write(unit=112, fmt = "(A)") "ulimit -s unlimited"

      allocate(data_k(task%samples(1), task%samples(2), task%samples(3), product(task%integer_indices), product(task%continuous_indices)), &
              sdata_k(product(task%samples), product(task%integer_indices), product(task%continuous_indices)))

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k) !!TODO: Check if the comporobation for task%samples() == 1 slows down the code too much.

        TID = OMP_GET_THREAD_NUM()
        IF (TID .EQ. 0) THEN
          write(unit=112, fmt = "(A, I6, A)") "Running on ", OMP_GET_NUM_THREADS(), " threads."
        ENDIF

        !$OMP DO
        do ik1=1, task%samples(1)
          if (task%samples(1) == 1) then
            k(1) = 0.0_dp
          else
            k(1) = -0.5_dp + real(ik1 - 1,dp)/real(task%samples(1) - 1,dp)
          endif
          do ik2 = 1, task%samples(2)
            if (task%samples(2) == 1) then
              k(2) = 0.0_dp
            else
              k(2) = -0.5_dp + real(ik2 - 1,dp)/real(task%samples(2) - 1,dp)
            endif
            do ik3 = 1, task%samples(3)
              if (task%samples(3) == 1) then
                k(3) = 0.0_dp
              else
                k(3) = -0.5_dp + real(ik3 - 1,dp)/real(task%samples(3) - 1,dp)
              endif

              !Gather data.
              if (associated(task%local_calculator)) then
                data_k(ik1, ik2, ik3, :, 1) = task%local_calculator(task, system, k)
              elseif (associated(task%global_calculator)) then 
                data_k(ik1, ik2, ik3, :, :) = task%global_calculator(task, system, k)
              endif
              
            enddo
          enddo
        enddo
      !$OMP END DO
      !$OMP END PARALLEL

      write(unit=112, fmt = "(A)") "Sampling done. Starting integration with extrapolation method."

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
        !Pass data array to memory layout.
        call shrink_array(data_k(:, :, :, i, r), sdata_k(:, i, r), info)
        !Integrate, if possible extrapolation method.
        call integral_extrapolation(sdata_k(:, i, r), task%samples, (/-0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp/), task%result(i, r), info)
        enddo
      enddo

      if (info == 1) then
        write(unit=112, fmt="(A)") "Integral done. Extrapolation sucessfull."
      else
        write(unit=112, fmt="(A)") "Integral done. Extrapolation fail and return rectangle approximation."
      endif

      deallocate(data_k, sdata_k)

    elseif (task%method == "rectangle") then !Rectangle method approximation case.

      write(unit=112, fmt = "(A)") "Starting BZ sampling and integration subroutine C1M3 with rectangle method."
      write(unit=112, fmt = "(A)") "Integrating task: "//trim(task%name)//" in the BZ for the system "//trim(system%name)//"."
      write(unit=112, fmt = "(A)") "The required memory for the integration process is approximately,"
      write(unit=112, fmt = "(F15.3, A)") 16.0_dp*real(product(task%integer_indices)*product(task%continuous_indices),dp)/1024.0_dp**2, "MB."

      allocate(temp_res(product(task%integer_indices), product(task%continuous_indices)))
      temp_res = cmplx(0.0_dp, 0.0_dp, dp)

      !$OMP PARALLEL DEFAULT (SHARED) PRIVATE (k) REDUCTION (+: temp_res) !!TODO: Check if the comporobation for task%samples() == 1 slows down the code too much.

      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        write(unit=112, fmt="(A, I6, A)") "Running on ", OMP_GET_NUM_THREADS(), " threads."
      ENDIF

      !$OMP DO
      do ik1 = 1, task%samples(1)
        if (task%samples(1) == 1) then
          k(1) = 0.0_dp
        else
          k(1) = -0.5_dp + real(ik1 - 1,dp)/real(task%samples(1) - 1,dp)
        endif
        do ik2 = 1, task%samples(2)
          if (task%samples(2) == 1) then
            k(2) = 0.0_dp
          else
            k(2) = -0.5_dp + real(ik2 - 1,dp)/real(task%samples(2) - 1,dp)
          endif
          do ik3 = 1, task%samples(3)
            if (task%samples(3) == 1) then
              k(3) = 0.0_dp
            else
              k(3) = -0.5_dp + real(ik3 - 1,dp)/real(task%samples(3) - 1,dp)
            endif

            !Gather data.
            if (associated(task%local_calculator)) then
              temp_res(:, 1) = temp_res(:, 1) + task%local_calculator(task, system, k)
            elseif (associated(task%global_calculator)) then 
              temp_res = temp_res + task%global_calculator(task, system, k)
            endif
            
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

  subroutine print_task_result(task, system)
    !Subroutine to format and output files related to the result of the task "task".
    type(BZ_integral_task), intent(in) :: task
    type(sys),              intent(in) :: system
    
    character(len=400) :: filename
    integer :: i_arr(size(task%integer_indices)), &
    r_arr(size(task%continuous_indices))
    integer :: i_mem, r_mem, count
    
    do i_mem = 1, product(task%integer_indices) !For each integer index.
    
    i_arr = integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.
    
    filename = trim(system%name)//'-'//trim(task%name)//'_'
    do count = 1, size(task%integer_indices)
    filename = trim(filename)//achar(48 + i_arr(count))
    enddo
    filename = trim(filename)//'.dat'
    
    open(unit=111, action="write", file=filename)
    
    do r_mem = 1, product(task%continuous_indices) !For each continuous index.
    
    r_arr = continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.
    
    write(unit=111, fmt=*) (task%ext_var_data(count)%data(r_arr(count)), count = 1, size(task%continuous_indices)),&
                            real(task%result(i_mem, r_mem), dp), aimag(task%result(i_mem, r_mem))!TODO: SET FORMAT.
    
    enddo
    
    close(unit=111)
    
    enddo
    
    end subroutine print_task_result

end module integrator