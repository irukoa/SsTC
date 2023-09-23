#include "cond_comp.h"
module SsTC_integrator

  USE OMP_LIB
  USE MPI_F08

  use SsTC_utility
  use SsTC_mpi_comms
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

    implicit none

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

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Creating BZ integral task "//trim(task%name)//"."

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
        if (rank == 0) write (unit=stdout, fmt="(a)") "          Warning: To employ the extrapolation method &
          & all the elements of the array 'samples' must be expressible as either 1 &
          & or 2^n + 1 for n = 0, 1, ..."
      elseif (method == "rectangle") then
        task%method = "rectangle"
      else
        if (rank == 0) write (unit=stdout, fmt="(a)") "          Integration method not recognized. Setting up rectangle method"
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

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Done."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_BZ_integral_task_constructor

  subroutine SsTC_sample_and_integrate_BZ_integral_task(task, system)

    implicit none

    class(SsTC_BZ_integral_task), intent(inout) :: task
    type(SsTC_sys), intent(in)                  :: system

    complex(kind=dp), allocatable :: data_k(:, :, :), &
                                     simd_tmp(:, :), &
                                     local_data_k(:, :, :), &
                                     temp_res(:, :), res(:, :)

    real(kind=dp) :: k(3)

    integer :: ik, k_ind(3), ik1, ik2, ik3, &
               info, &
               i, r

    real(kind=dp) :: start_time, end_time

    logical :: error = .false.

    type(SsTC_local_k_data) :: sampling_info

    integer, allocatable :: counts(:), displs(:)

    start_time = MPI_WTIME()

    allocate (sampling_info%integer_indices(3))
    sampling_info%integer_indices = task%samples

    allocate (counts(0:nProcs - 1), displs(0:nProcs - 1))
    call get_MPI_task_partition(product(task%samples), nProcs, counts, displs)

    if (task%method == "extrapolation") then !Extrapolation case.

      if (rank == 0) write (unit=stdout, fmt="(a)") &
        "          Starting BZ sampling and integration subroutine with extrapolation method."
      if (rank == 0) write (unit=stdout, fmt="(a)") "          Integrating task: "//trim(task%name)//" in the BZ for the system " &
        //trim(system%name)//"."

      allocate (local_data_k(displs(rank) + 1:displs(rank) + counts(rank), &
                             product(task%integer_indices), product(task%continuous_indices)), &
                data_k(product(task%samples), &
                       product(task%integer_indices), product(task%continuous_indices)), &
                simd_tmp(product(task%integer_indices), product(task%continuous_indices)))

      local_data_k = cmplx(0.0_dp, 0.0_dp, dp)
      data_k = cmplx(0.0_dp, 0.0_dp, dp)
      simd_tmp = cmplx(0.0_dp, 0.0_dp, dp)

      !_OMPOFFLOADTGT_(TARGET TEAMS)
      !_OMPTGT_(PARALLEL DO REDUCTION (.or.: error) &)
      !_OMPTGT_(SHARED(task, system, displs, counts, rank, sampling_info, local_data_k) &)
      !_OMPTGT_(PRIVATE(ik, k_ind, ik1, ik2, ik3, k, simd_tmp))
      !_OMPTGT_(SIMD)
      do ik = displs(rank) + 1, displs(rank) + counts(rank)

        k_ind = SsTC_integer_memory_element_to_array_element(sampling_info, ik)
        ik1 = k_ind(1)
        ik2 = k_ind(2)
        ik3 = k_ind(3)

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
          simd_tmp(:, 1) = task%local_calculator(task, system, k, error)
          local_data_k(ik, :, 1) = simd_tmp(:, 1)
        elseif (associated(task%global_calculator)) then
          simd_tmp = task%global_calculator(task, system, k, error)
          local_data_k(ik, :, :) = simd_tmp
        endif

      enddo
      !_OMPOFFLOADTGT_(END TARGET)

      if (error) write (unit=stderr, fmt="(a, i5, a)") "Rank: ", rank, ". ERROR in subroutine &
      & SsTC_sample_and_integrate_BZ_integral_task. Check error log."

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
          call MPI_ALLGATHERV(local_data_k(:, i, r), size(local_data_k(:, i, r)), MPI_COMPLEX16, data_k(:, i, r), &
                              counts, &
                              displs, &
                              MPI_COMPLEX16, MPI_COMM_WORLD, ierror)
        enddo
      enddo

      if (rank == 0) write (unit=stdout, fmt="(a)") "          Sampling done. Starting integration with extrapolation method."

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
          !Integrate, if possible extrapolation method.
          call integral_extrapolation(data_k(:, i, r), task%samples, &
                                      (/-0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp, -0.5_dp, 0.5_dp/), task%result(i, r), info)
        enddo
      enddo

      if (info == 1) then
        if (rank == 0) write (unit=stdout, fmt="(a)") &
        &"          Integral done. Extrapolation sucessfull."
      else
        if (rank == 0) write (unit=stdout, fmt="(a)") &
        &"          Integral done. Extrapolation failed. Returning rectangle approximation."
      endif
      if (rank == 0) write (unit=stdout, fmt="(a)") ""

      deallocate (local_data_k, data_k, simd_tmp)

    elseif (task%method == "rectangle") then !Rectangle method approximation case.

      if (rank == 0) write (unit=stdout, fmt="(a)") &
        "          Starting BZ sampling and integration subroutine with rectangle method."
      if (rank == 0) write (unit=stdout, fmt="(a)") "          Integrating task: "//trim(task%name)//" in the BZ for the system "&
      &//trim(system%name)//"."

      allocate (temp_res(product(task%integer_indices), product(task%continuous_indices)), &
                res(product(task%integer_indices), product(task%continuous_indices)), &
                simd_tmp(product(task%integer_indices), product(task%continuous_indices)))

      !_OMPOFFLOADTGT_(TARGET TEAMS)
      !_OMPTGT_(PARALLEL DO REDUCTION (+: temp_res) REDUCTION (.or.: error) &)
      !_OMPTGT_(SHARED(task, system, displs, counts, rank, sampling_info) &)
      !_OMPTGT_(PRIVATE(ik, k_ind, ik1, ik2, ik3, k, simd_tmp))
      !_OMPTGT_(SIMD)
      do ik = displs(rank) + 1, displs(rank) + counts(rank)

        k_ind = SsTC_integer_memory_element_to_array_element(sampling_info, ik)
        ik1 = k_ind(1)
        ik2 = k_ind(2)
        ik3 = k_ind(3)

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
          simd_tmp(:, 1) = task%local_calculator(task, system, k, error)
          temp_res(:, 1) = temp_res(:, 1) + simd_tmp(:, 1)
        elseif (associated(task%global_calculator)) then
          simd_tmp = task%global_calculator(task, system, k, error)
          temp_res = temp_res + simd_tmp
        endif

      enddo
      !_OMPOFFLOADTGT_(END TARGET)

      if (error) write (unit=stderr, fmt="(a, i5, a)") "Rank: ", rank, ". ERROR in subroutine &
      & SsTC_sample_and_integrate_BZ_integral_task. Check error log."

      call MPI_ALLREDUCE(temp_res, res, size(temp_res), MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD, ierror)
      task%result = res/product(task%samples)

      if (rank == 0) write (unit=stdout, fmt="(a)") "          Integral done."

      deallocate (temp_res, res, simd_tmp)
      deallocate (sampling_info%integer_indices, counts, displs)

      if (rank == 0) write (unit=stdout, fmt="(a)") ""
    endif

    end_time = MPI_WTIME() !End timer.
    if (rank == 0) write (unit=stdout, fmt="(a, f15.3, a)") "          Total execution time: ", end_time - start_time, " s."

  end subroutine SsTC_sample_and_integrate_BZ_integral_task

  subroutine SsTC_print_BZ_integral_task(task, system)

    implicit none

    !Subroutine to format and output files related to the result of the task "task".
    class(SsTC_BZ_integral_task), intent(in) :: task
    type(SsTC_sys), intent(in)               :: system

    character(len=400) :: filename, fmtf, num_label
    integer            :: i_arr(size(task%integer_indices)), &
                          r_arr(size(task%continuous_indices))
    integer            :: i_mem, r_mem, count
    integer            :: printunit

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Printing BZ integral task: "&
    &//trim(task%name)//" for the system "//trim(system%name)//"."

    if (associated(task%local_calculator)) then

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count = 1, size(task%integer_indices)
          write (num_label, fmt="(i0)") i_arr(count)
          filename = trim(filename)//trim(num_label)
        enddo
        filename = trim(filename)//'.dat'

        if (rank == 0) open (newunit=printunit, action="write", file=filename)

        if (rank == 0) write (unit=printunit, fmt="(2e18.8e3)") &
        &real(task%result(i_mem, 1), dp), aimag(task%result(i_mem, 1))

        if (rank == 0) close (unit=printunit)

      enddo

    elseif (associated(task%global_calculator)) then

      write (fmtf, "(i10)") size(task%continuous_indices) + 2
      fmtf = '('//trim(adjustl(trim(fmtf)))//'e18.8e3)'

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count = 1, size(task%integer_indices)
          write (num_label, fmt="(i0)") i_arr(count)
          filename = trim(filename)//trim(num_label)
        enddo
        filename = trim(filename)//'.dat'

        if (rank == 0) open (newunit=printunit, action="write", file=filename)

        do r_mem = 1, product(task%continuous_indices) !For each continuous index.

          r_arr = SsTC_continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.

          if (rank == 0) write (unit=printunit, fmt=fmtf) &
            (task%ext_var_data(count)%data(r_arr(count)), count=1, size(task%continuous_indices)), &
            real(task%result(i_mem, r_mem), dp), aimag(task%result(i_mem, r_mem))

        enddo

        if (rank == 0) close (unit=printunit)

      enddo

    endif

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Printing done."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_print_BZ_integral_task

end module SsTC_integrator
