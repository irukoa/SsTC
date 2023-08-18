#include 'cond_comp.h'
module SsTC_integrator

  USE OMP_LIB

  use SsTC_utility
  use SsTC_data_structures
  use extrapolation_integration

  implicit none

  private

  type, extends(SsTC_global_k_data) :: SsTC_BZ_integral_task
    !Integration method.
    character(len=120)                             :: method
    !Integration samples.
    integer                                        :: samples(3)
    !Result of the integration, contains the integer
    !index and the continuous index in memory layout, respectively.
    complex(kind=dp), allocatable                  :: result(:, :)
  end type SsTC_BZ_integral_task

  public :: SsTC_BZ_integral_task

  public :: SsTC_BZ_integral_task_constructor
  public :: SsTC_sample_and_integrate_BZ_integral_task
  public :: SsTC_print_BZ_integral_task

contains

  subroutine SsTC_BZ_integral_task_constructor(task, name, &
                                               l_calculator, g_calculator, &
                                               method, samples, &
                                               N_int_ind, int_ind_range, &
                                               N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                                               part_int_comp)

    character(len=*) :: name

    procedure(SsTC_local_calculator), optional  :: l_calculator
    procedure(SsTC_global_calculator), optional :: g_calculator

    character(len=*), optional, intent(in) :: method
    integer, optional, intent(in)          :: samples(3)

    integer, optional, intent(in) :: N_int_ind
    integer, optional, intent(in) :: int_ind_range(:)

    integer, optional, intent(in)       :: N_ext_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(:), ext_vars_end(:)
    integer, optional, intent(in)       :: ext_vars_steps(:)

    integer, optional, intent(in) :: part_int_comp(:)

    class(SsTC_BZ_integral_task), intent(out) :: task

    integer :: i

    !Set name.
    task%name = name

    write (unit=stdout, fmt="(a)") "          Creating BZ integral task "//trim(task%name)//"."

    !Set integer index data.
    if (((N_int_ind) .ge. 1) .and. (present(int_ind_range))) then
      allocate (task%integer_indices(N_int_ind))
      do i = 1, N_int_ind
        task%integer_indices(i) = int_ind_range(i)
      enddo
    else
      allocate (task%integer_indices(1))
      task%integer_indices(1) = 1
    endif

    !Set external variable data.
    if (((N_ext_vars) .ge. 1) .and. (present(ext_vars_start)) .and. (present(ext_vars_end)) .and. (present(ext_vars_steps))) then
      allocate (task%continuous_indices(N_ext_vars), task%ext_var_data(N_ext_vars))
      do i = 1, N_ext_vars
        task%continuous_indices(i) = ext_vars_steps(i)
        allocate (task%ext_var_data(i)%data(ext_vars_steps(i)))
        task%ext_var_data(i) = SsTC_external_variable_constructor(ext_vars_start(i), ext_vars_end(i), ext_vars_steps(i))
      enddo
    else
      allocate (task%continuous_indices(1), task%ext_var_data(1))
      task%continuous_indices(1) = 1.0_dp
      allocate (task%ext_var_data(1)%data(1))
      task%ext_var_data(1)%data = (/1.0_dp/)
    endif

    if (present(g_calculator) .and. (.not. present(l_calculator))) then
      !Set calculator pointer (function alias).
      task%global_calculator => g_calculator
      nullify (task%local_calculator)
    elseif (present(l_calculator) .and. (.not. present(g_calculator))) then
      !Set calculator pointer (function alias).
      task%local_calculator => l_calculator
      nullify (task%global_calculator)
    endif

    !Set result.
    allocate (task%result(product(task%integer_indices), product(task%continuous_indices)))
    task%result = cmplx_0

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = &
      SsTC_integer_array_element_to_memory_element(task, part_int_comp)

    !Set integration method.
    if (present(method)) then
      if (method == "extrapolation") then
        task%method = "extrapolation"
        write (unit=stdout, fmt="(a)") "          Warning: To employ the extrapolation method &
          & all the elements of the array 'samples' must be expressible as either 1 &
          & or 2^n + 1 for n = 0, 1, ..."
      elseif (method == "rectangle") then
        task%method = "rectangle"
      else
        write (unit=stdout, fmt="(a)") "          Integration method not recognized. Setting up rectangle method"
        task%method = "rectangle"
      endif
    else
      task%method = "rectangle"
    endif

    !Set number of integration samples (discretization of BZ).
    if (present(samples)) then
      task%samples = samples
    else
      task%samples = (/10, 10, 10/)
    endif

    write (unit=stdout, fmt="(a)") "          Done."
    write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_BZ_integral_task_constructor

  subroutine SsTC_sample_and_integrate_BZ_integral_task(task, system)

    class(SsTC_BZ_integral_task), intent(inout) :: task
    type(SsTC_sys), intent(in)                  :: system

    complex(kind=dp), allocatable :: data_k(:, :, :, :, :), &
                                     sdata_k(:, :, :), &
                                     temp_res(:, :)

    real(kind=dp) :: k(3)

    integer :: ik1, ik2, ik3, &
               info, &
               i, r

    integer :: TID, report_step
    integer :: progress = 0
    real(kind=dp) :: start_time, end_time

    logical :: error = .false.

    report_step = nint(real(product(task%samples)/100, dp)) + 1

    start_time = omp_get_wtime() !Start timer.

    if (task%method == "extrapolation") then !Extraplation case. Requires large RAM.

      write (unit=stdout, fmt="(a)") "          Starting BZ sampling and integration subroutine with extrapolation method."
      write (unit=stdout, fmt="(a)") "          Integrating task: "//trim(task%name)//" in the BZ for the system "&
        &//trim(system%name)//"."
      write (unit=stdout, fmt="(a)") "          The required memory for the integration process is approximately,"
      write (unit=stdout, fmt="(a, f15.3, a)") "          ", 16.0_dp*real(product(task%samples)*product(task%integer_indices)*&
      & product(task%continuous_indices), dp)/1024.0_dp**2, "MB."
      write (unit=stdout, fmt="(a)") "          Some computers limit the maximum memory an array can allocate."
      write (unit=stdout, fmt="(a)") "          If this is your case and SIGSEGV triggers&
      & try using the next command before execution:"
      write (unit=stdout, fmt="(a)") "          ulimit -s unlimited"

      allocate (data_k(task%samples(1), task%samples(2), task%samples(3), &
                       product(task%integer_indices), product(task%continuous_indices)), &
                sdata_k(product(task%samples), &
                        product(task%integer_indices), product(task%continuous_indices)))

      !_OMPTGT_(PARALLEL DEFAULT(SHARED) PRIVATE(k))

      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        write (unit=stdout, fmt="(a, i5, a)") "         Running on ", OMP_GET_NUM_THREADS(), " threads."
      ENDIF

      !_OMPTGT_(DO COLLAPSE(3))
      do ik1 = 1, task%samples(1)
        do ik2 = 1, task%samples(2)
          do ik3 = 1, task%samples(3)

            if (task%samples(1) == 1) then
              k(1) = 0.0_dp
            else
              k(1) = -0.5_dp + real(ik1 - 1, dp)/real(task%samples(1) - 1, dp)
            endif
            if (task%samples(2) == 1) then
              k(2) = 0.0_dp
            else
              k(2) = -0.5_dp + real(ik2 - 1, dp)/real(task%samples(2) - 1, dp)
            endif
            if (task%samples(3) == 1) then
              k(3) = 0.0_dp
            else
              k(3) = -0.5_dp + real(ik3 - 1, dp)/real(task%samples(3) - 1, dp)
            endif

            !Gather data.
            if (associated(task%local_calculator)) then
              data_k(ik1, ik2, ik3, :, 1) = task%local_calculator(task, system, k, error)
            elseif (associated(task%global_calculator)) then
              data_k(ik1, ik2, ik3, :, :) = task%global_calculator(task, system, k, error)
            endif

            if (error) then
              write (unit=stderr, fmt="(a, 3e18.8e3)") "Error when sampling k-point", k
              write (unit=stderr, fmt="(a)") "Stopping..."
              stop
            endif

            !_OMPTGT_(ATOMIC UPDATE)
            progress = progress + 1

            if (modulo(progress, report_step) == report_step/2) then !Update progress.
              write (unit=stdout, fmt="(a, i12, a, i12, a)") &
              & "          Progress: ", progress, "/", product(task%samples),&
              & " kpts sampled."
            endif

          enddo
        enddo
      enddo
      !_OMPTGT_(END DO)
      !_OMPTGT_(END PARALLEL)
      progress = 0

      write (unit=stdout, fmt="(a)") "          Sampling done. Starting integration with extrapolation method."

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
          !Pass data array to memory layout.
          call shrink_array(data_k(:, :, :, i, r), sdata_k(:, i, r), info)
          !Integrate, if possible extrapolation method.
          call integral_extrapolation(sdata_k(:, i, r), task%samples, &
                                      (/-0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp/), task%result(i, r), info)
        enddo
      enddo

      if (info == 1) then
        write (unit=stdout, fmt="(a)") "          Integral done. Extrapolation sucessfull."
      else
        write (unit=stdout, fmt="(a)") "          Integral done. Extrapolation failed. Returning rectangle approximation."
      endif
      write (unit=stdout, fmt="(a)") ""

      deallocate (data_k, sdata_k)

    elseif (task%method == "rectangle") then !Rectangle method approximation case.

      write (unit=stdout, fmt="(a)") "          Starting BZ sampling and integration subroutine with rectangle method."
      write (unit=stdout, fmt="(a)") "          Integrating task: "//trim(task%name)//" in the BZ for the system "&
      &//trim(system%name)//"."
      write (unit=stdout, fmt="(a)") "          The required memory for the integration process is approximately,"
      write (unit=stdout, fmt="(a, f15.3, a)") "          ", 16.0_dp*real(product(task%integer_indices)*&
      & product(task%continuous_indices), dp)/1024.0_dp**2, "MB."

      allocate (temp_res(product(task%integer_indices), product(task%continuous_indices)))
      temp_res = cmplx(0.0_dp, 0.0_dp, dp)

      !_OMPTGT_(PARALLEL DEFAULT (SHARED) PRIVATE (k) REDUCTION (+: temp_res))

      TID = OMP_GET_THREAD_NUM()
      IF (TID .EQ. 0) THEN
        write (unit=stdout, fmt="(a, I5, a)") "         Running on ", OMP_GET_NUM_THREADS(), " threads."
      ENDIF

      !_OMPTGT_(DO COLLAPSE(3))
      do ik1 = 1, task%samples(1)
        do ik2 = 1, task%samples(2)
          do ik3 = 1, task%samples(3)

            if (task%samples(1) == 1) then
              k(1) = 0.0_dp
            else
              k(1) = -0.5_dp + real(ik1 - 1, dp)/real(task%samples(1) - 1, dp)
            endif
            if (task%samples(2) == 1) then
              k(2) = 0.0_dp
            else
              k(2) = -0.5_dp + real(ik2 - 1, dp)/real(task%samples(2) - 1, dp)
            endif
            if (task%samples(3) == 1) then
              k(3) = 0.0_dp
            else
              k(3) = -0.5_dp + real(ik3 - 1, dp)/real(task%samples(3) - 1, dp)
            endif

            !Gather data.
            if (associated(task%local_calculator)) then
              temp_res(:, 1) = temp_res(:, 1) + task%local_calculator(task, system, k, error)
            elseif (associated(task%global_calculator)) then
              temp_res = temp_res + task%global_calculator(task, system, k, error)
            endif

            if (error) then
              write (unit=stderr, fmt="(a, 3e18.8e3)") "Error when sampling k-point", k
              write (unit=stderr, fmt="(a)") "Stopping..."
              stop
            endif

            !_OMPTGT_(ATOMIC UPDATE)
            progress = progress + 1

            if (modulo(progress, report_step) == report_step/2) then !Update progress.
              write (unit=stdout, fmt="(a, i12, a, i12, a)") &
              &"          Progress: ", progress, "/", product(task%samples), &
              &" kpts sampled."
            endif

          enddo
        enddo
      enddo
      !_OMPTGT_(END DO)
      !_OMPTGT_(END PARALLEL)
      progress = 0

      end_time = omp_get_wtime() !End timer.

      write (unit=stdout, fmt="(a)") "          Integral done."
      write (unit=stdout, fmt="(a, f15.3, a)") "          Total execution time: ", end_time - start_time, " s."

      task%result = temp_res/product(task%samples)

      deallocate (temp_res)

      write (unit=stdout, fmt="(a)") ""
    endif

  end subroutine SsTC_sample_and_integrate_BZ_integral_task

  subroutine SsTC_print_BZ_integral_task(task, system)
    !Subroutine to format and output files related to the result of the task "task".
    class(SsTC_BZ_integral_task), intent(in) :: task
    type(SsTC_sys), intent(in)               :: system

    character(len=400) :: filename, fmtf
    integer            :: i_arr(size(task%integer_indices)), &
                          r_arr(size(task%continuous_indices))
    integer            :: i_mem, r_mem, count
    integer            :: printunit

    write (unit=stdout, fmt="(a)") "          Printing BZ integral task: "&
    &//trim(task%name)//" for the system "//trim(system%name)//"."

    if (associated(task%local_calculator)) then

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count = 1, size(task%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count))
        enddo
        filename = trim(filename)//'.dat'

        open (newunit=printunit, action="write", file=filename)

        write (unit=printunit, fmt="(2e18.8e3)") real(task%result(i_mem, 1), dp), aimag(task%result(i_mem, 1))

        close (unit=printunit)

      enddo

    elseif (associated(task%global_calculator)) then

      write (fmtf, "(i10)") size(task%continuous_indices) + 2
      fmtf = '('//trim(adjustl(trim(fmtf)))//'e18.8e3)'

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count = 1, size(task%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count))
        enddo
        filename = trim(filename)//'.dat'

        open (newunit=printunit, action="write", file=filename)

        do r_mem = 1, product(task%continuous_indices) !For each continuous index.

          r_arr = SsTC_continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.

          write (unit=printunit, fmt=fmtf) (task%ext_var_data(count)%data(r_arr(count)), count=1, size(task%continuous_indices)), &
            real(task%result(i_mem, r_mem), dp), aimag(task%result(i_mem, r_mem))

        enddo

        close (unit=printunit)

      enddo

    endif

    write (unit=stdout, fmt="(a)") "          Printing done."
    write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_print_BZ_integral_task

end module SsTC_integrator
