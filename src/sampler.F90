#include "cond_comp.h"
module SsTC_sampler

  USE OMP_LIB
  USE MPI_F08

  use SsTC_utility
  use SsTC_mpi_comms
  use SsTC_data_structures
  use extrapolation_integration

  implicit none

  private

  type, extends(SsTC_global_k_data) :: SsTC_sampling_task
    integer                       :: samples(3) = (/10, 10, 10/)
    !Integer index, continuous index and kpt index 1, 2 and 3 respectively.
    complex(kind=dp), allocatable :: BZ_data(:, :, :, :, :)
    logical                       :: predefined_set_of_kpts = .false.
    real(kind=dp), allocatable    :: kpts(:, :)
    !Integer index, continuous index and kpt index respectively.
    complex(kind=dp), allocatable :: predefined_sampled_data(:, :, :)
  end type SsTC_sampling_task

  public :: SsTC_sampling_task

  public :: SsTC_sampling_task_constructor
  public :: SsTC_sample_sampling_task
  public :: SsTC_print_sampling

contains

  subroutine SsTC_sampling_task_constructor(task, name, &
                                            l_calculator, g_calculator, &
                                            samples, nkpts, kpts, &
                                            N_int_ind, int_ind_range, &
                                            N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                                            part_int_comp)

    implicit none

    character(len=*) :: name

    procedure(SsTC_local_calculator), optional  :: l_calculator
    procedure(SsTC_global_calculator), optional :: g_calculator

    integer, optional, intent(in) :: samples(3)
    integer, optional, intent(in) :: nkpts
    real(kind=dp), optional, intent(in) :: kpts(:, :)

    integer, intent(in)           :: N_int_ind
    integer, optional, intent(in) :: int_ind_range(N_int_ind)

    integer, intent(in)                 :: N_ext_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(N_ext_vars), ext_vars_end(N_ext_vars)
    integer, optional, intent(in)       :: ext_vars_steps(N_ext_vars)

    integer, optional, intent(in) :: part_int_comp(N_int_ind)

    class(SsTC_sampling_task), intent(out) :: task

    integer :: i

    !Set name.
    task%name = name

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Creating BZ sampling task "//trim(task%name)//"."

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
    else
      nullify (task%local_calculator)
      nullify (task%global_calculator)
    endif

    !Set number of samples (discretization of BZ).
    if (present(samples)) task%samples = samples

    !Set list of kpoiints and results for the case of a predefined set of k-points.
    if ((present(nkpts)) .and. (present(kpts))) then
      if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") &
        "          This task will sample only over the predefined set of provided k points."
      task%predefined_set_of_kpts = .true.
      allocate (task%kpts(nkpts, 3))
      task%kpts = kpts
      allocate (task%predefined_sampled_data(product(task%integer_indices), product(task%continuous_indices), &
                                             size(task%kpts(:, 1))))
      task%predefined_sampled_data = cmplx_0
    endif

    !Set result.
    if (.not. task%predefined_set_of_kpts) then
      allocate (task%BZ_data(product(task%integer_indices), product(task%continuous_indices), &
                             task%samples(1), task%samples(2), task%samples(3)))
      task%BZ_data = cmplx_0
    endif

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = &
      SsTC_integer_array_element_to_memory_element(task, part_int_comp)

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Done."
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_sampling_task_constructor

  subroutine SsTC_sample_sampling_task(task, system)

    implicit none

    class(SsTC_sampling_task), intent(inout) :: task
    type(SsTC_sys), intent(in)               :: system

    complex(kind=dp), allocatable :: data_k(:, :, :), &
                                     local_data_k(:, :, :)

    real(kind=dp) :: k(3)

    integer :: ik, k_ind(3), &
               ik1, ik2, ik3, &
               i, r, info

    real(kind=dp) :: start_time, end_time

    logical :: error = .false.

    type(SsTC_local_k_data) :: sampling_info

    integer, allocatable :: counts(:), displs(:)

    start_time = MPI_WTIME() !Start timer.

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Starting BZ sampling subroutine."

    if (.not. task%predefined_set_of_kpts) then !Sampling on whole BZ specified by task%samples.

      allocate (sampling_info%integer_indices(3))
      sampling_info%integer_indices = task%samples

      allocate (counts(0:nProcs - 1), displs(0:nProcs - 1))
      call get_MPI_task_partition(product(task%samples), nProcs, counts, displs)

      if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a, a, a, a, a)") "          Sampling task: "//trim(task%name)// &
        " in the BZ for the system "//trim(system%name)//"."

      allocate (local_data_k(displs(rank) + 1:displs(rank) + counts(rank), &
                             product(task%integer_indices), product(task%continuous_indices)), &
                data_k(product(task%samples), &
                       product(task%integer_indices), product(task%continuous_indices)))

      local_data_k = cmplx(0.0_dp, 0.0_dp, dp)
      data_k = cmplx(0.0_dp, 0.0_dp, dp)

      !_OMPOFFLOADTGT_(TARGET TEAMS)
      !_OMPTGT_(PARALLEL DO REDUCTION (.or.: error) &)
      !_OMPTGT_(SHARED(task, system, displs, counts, rank, sampling_info, local_data_k) &)
      !_OMPTGT_(PRIVATE(ik, k_ind, ik1, ik2, ik3, k))
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
          local_data_k(ik, :, 1) = task%local_calculator(task, system, k, error)
        elseif (associated(task%global_calculator)) then
          local_data_k(ik, :, :) = task%global_calculator(task, system, k, error)
        endif

      enddo
      !_OMPOFFLOADTGT_(END TARGET)

      if (error) write (unit=stderr, fmt="(a, i5, a)") "Rank: ", rank, ". ERROR in subroutine &
      & SsTC_sample_sampling_task. Check error log."

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
          call MPI_ALLGATHERV(local_data_k(:, i, r), size(local_data_k(:, i, r)), MPI_COMPLEX16, data_k(:, i, r), &
                              counts, &
                              displs, &
                              MPI_COMPLEX16, MPI_COMM_WORLD, ierror)
        enddo
      enddo

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
          call expand_array(data_k(:, i, r), task%BZ_data(i, r, :, :, :), info)
        enddo
      enddo

      deallocate (local_data_k, data_k, counts, displs)

    else !Sampling along a predefined set of kpoints specified by task%kpts.

      allocate (counts(0:nProcs - 1), displs(0:nProcs - 1))
      call get_MPI_task_partition(size(task%kpts(:, 1)), nProcs, counts, displs)

      if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a, a, a, a, a)") "          Sampling task: "//trim(task%name)// &
        " in a predefined set of kpoints for the system "//trim(system%name)//"."

      allocate (local_data_k(displs(rank) + 1:displs(rank) + counts(rank), &
                             product(task%integer_indices), product(task%continuous_indices)), &
                data_k(size(task%kpts(:, 1)), &
                       product(task%integer_indices), product(task%continuous_indices)))

      local_data_k = cmplx(0.0_dp, 0.0_dp, dp)

      !_OMPOFFLOADTGT_(TARGET TEAMS)
      !_OMPTGT_(PARALLEL DO REDUCTION (.or.: error) &)
      !_OMPTGT_(SHARED(task, system, displs, counts, rank, local_data_k) &)
      !_OMPTGT_(PRIVATE(ik, k))
      do ik = displs(rank) + 1, displs(rank) + counts(rank)

        k = task%kpts(ik, :)

        !Gather data.
        if (associated(task%local_calculator)) then
          local_data_k(ik, :, 1) = task%local_calculator(task, system, k, error)
        elseif (associated(task%global_calculator)) then
          local_data_k(ik, :, :) = task%global_calculator(task, system, k, error)
        endif

      enddo
      !_OMPOFFLOADTGT_(END TARGET)

      do i = 1, product(task%integer_indices) !For each integer index.
        do r = 1, product(task%continuous_indices) !For each continuous index.
          call MPI_ALLGATHERV(local_data_k(:, i, r), size(local_data_k(:, i, r)), MPI_COMPLEX16, data_k(:, i, r), &
                              counts, &
                              displs, &
                              MPI_COMPLEX16, MPI_COMM_WORLD, ierror)
          task%predefined_sampled_data(i, r, :) = data_k(:, i, r)
        enddo
      enddo

      deallocate (local_data_k, data_k, counts, displs)

    endif

    end_time = MPI_WTIME() !End timer.

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Sampling done."
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a, f15.3, a)") &
      "          Total execution time: ", end_time - start_time, " s."
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_sample_sampling_task

  subroutine SsTC_print_sampling(task, system)

    implicit none

    !Subroutine to format and output files related to the result of the task "task".
    class(SsTC_sampling_task), intent(in) :: task
    type(SsTC_sys), intent(in)            :: system

    real(kind=dp) :: k(3)

    character(len=400) :: filename, fmtf, num_label
    integer            :: i_arr(size(task%integer_indices)), &
                          r_arr(size(task%continuous_indices))
    integer            :: i_mem, r_mem, &
                          count, &
                          ik, ik1, ik2, ik3
    integer            :: printunit

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") &
      "          Printing sampling task: "//trim(task%name)//" for the system "//trim(system%name)//"."

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

        if (.not. task%predefined_set_of_kpts) then !Print the sampling on whole BZ specified by task%samples.

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

                if (rank == 0) write (unit=printunit, fmt="(5e18.8e3)") k, real(task%BZ_data(i_mem, 1, ik1, ik2, ik3), dp), &
                  aimag(task%BZ_data(i_mem, 1, ik1, ik2, ik3))

              enddo
            enddo
          enddo

        else !Print the sampling along a predefined set of kpoints specified by task%kpts.

          do ik = 1, size(task%kpts(:, 1))

            k = task%kpts(ik, :)

            if (rank == 0) write (unit=printunit, fmt="(i10, 5e18.8e3)") ik, k, &
              real(task%predefined_sampled_data(i_mem, 1, ik), dp), &
              aimag(task%predefined_sampled_data(i_mem, 1, ik))

          enddo

        endif

        if (rank == 0) close (unit=printunit)

      enddo

    elseif (associated(task%global_calculator)) then

      write (fmtf, "(I10)") size(task%continuous_indices) + 5
      if (task%predefined_set_of_kpts) then
        fmtf = '(i10, '//trim(adjustl(trim(fmtf)))//'E18.8E3)'
      else
        fmtf = '('//trim(adjustl(trim(fmtf)))//'E18.8E3)'
      endif

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

          if (.not. task%predefined_set_of_kpts) then !Print the sampling on whole BZ specified by task%samples.

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

                  if (rank == 0) write (unit=printunit, fmt=fmtf) k, &
                    (task%ext_var_data(count)%data(r_arr(count)), count=1, size(task%continuous_indices)), &
                    real(task%BZ_data(i_mem, r_mem, ik1, ik2, ik3), dp), aimag(task%BZ_data(i_mem, r_mem, ik1, ik2, ik3))

                enddo
              enddo
            enddo

          else !Print the sampling along a predefined set of kpoints specified by task%kpts.

            do ik = 1, size(task%kpts(:, 1))

              k = task%kpts(ik, :)

              if (rank == 0) write (unit=printunit, fmt=fmtf) ik, k, &
                (task%ext_var_data(count)%data(r_arr(count)), count=1, size(task%continuous_indices)), &
                real(task%predefined_sampled_data(i_mem, r_mem, ik), dp), aimag(task%predefined_sampled_data(i_mem, r_mem, ik))

            enddo

          endif

        enddo

        if (rank == 0) close (unit=printunit)

      enddo

    endif

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Printing done."
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_print_sampling

end module SsTC_sampler
