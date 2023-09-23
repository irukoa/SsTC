#include "cond_comp.h"
module SsTC_kslice

  USE OMP_LIB
  USE MPI_F08

  use SsTC_utility
  use SsTC_mpi_comms
  use SsTC_data_structures
  use extrapolation_integration

  implicit none

  private

  type, extends(SsTC_global_k_data) :: SsTC_kslice_task
    !Default: sample k_z = 0 slice in a 100x100 mesh.
    real(kind=dp) :: corner(3) = (/-0.5_dp, -0.5_dp, 0.0_dp/), &
                     vector(2, 3) = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/), (/2, 3/))
    integer       :: samples(2) = (/100, 100/)
    !Integer index, continuous index and kpt index 1 and 2 respectively.
    complex(kind=dp), allocatable :: kslice_data(:, :, :, :)
  end type SsTC_kslice_task

  public :: SsTC_kslice_task

  public :: SsTC_kslice_task_constructor
  public :: SsTC_sample_kslice_task
  public :: SsTC_print_kslice

contains

  subroutine SsTC_kslice_task_constructor(task, name, &
                                          l_calculator, g_calculator, &
                                          corner, vector_a, vector_b, samples, &
                                          N_int_ind, int_ind_range, &
                                          N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                                          part_int_comp)

    implicit none

    character(len=*) :: name

    procedure(SsTC_local_calculator), optional  :: l_calculator
    procedure(SsTC_global_calculator), optional :: g_calculator

    real(kind=dp), optional, intent(in) :: corner(3), vector_a(3), vector_b(3)
    integer, optional, intent(in)       :: samples(2)

    integer, optional, intent(in) :: N_int_ind
    integer, optional, intent(in) :: int_ind_range(:)

    integer, optional, intent(in)       :: N_ext_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(:), ext_vars_end(:)
    integer, optional, intent(in)       :: ext_vars_steps(:)

    integer, optional, intent(in) :: part_int_comp(:)

    class(SsTC_kslice_task), intent(out) :: task

    logical :: error = .False.

    integer          :: i
    complex(kind=dp) :: lindep(2, 3)
    real(kind=dp)    :: sigma(2, 3)
    logical          :: dep = .False.

    !Set name.
    task%name = name

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Creating kpath task "//trim(task%name)//"."

    !Set kslice info.
    if (present(corner)) task%corner = corner

    if (present(vector_a) .and. present(vector_b)) then
      !Check linear depenbdence:
      do i = 1, 3
        lindep(1, i) = cmplx(vector_a(i), 0.0_dp, dp)
        lindep(2, i) = cmplx(vector_b(i), 0.0_dp, dp)
      enddo

      !get SVD decomp,
      call SsTC_utility_SVD(lindep, sigma, error)
      if (error) then
        write (unit=stderr, fmt="(a, i5, a)") "          Rank ", rank, &
          ": Error in function kslice_task_constructor when checking linear dependence of vectors."
        write (unit=stderr, fmt="(a, i5, a)") "          Rank ", rank, ": Error in function kslice_task_constructor&
        & when checking linear dependence of vectors."
        if (rank == 0) write (unit=stdout, fmt="(a)") "Supposing linearly dependent input vectors."
        dep = .True.
        goto 2
      endif

      !check if rank (number of nonzero singular values) is at least 2.
      do i = 1, 2
        if (abs(sigma(i, i)) .LE. 1.0E-6_dp) dep = .True.
      enddo

2     continue
      if (dep) then
        if (rank == 0) write (unit=stdout, fmt="(a)") "          Selected vectors are not linearly independent."
        if (rank == 0) write (unit=stdout, fmt="(a)") "          Keeping default vectors (1, 0, 0), (0, 1, 0)."
      else
        task%vector(1, :) = vector_a
        task%vector(2, :) = vector_b
      endif
    endif

    !Set number of samples.
    if (present(samples)) then
      task%samples = samples
    endif
    if ((task%samples(1) .LE. 1) .or. (task%samples(2) .LE. 1)) task%samples = (/2, 2/)

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

    !Set kdata.
    allocate (task%kslice_data(product(task%integer_indices), product(task%continuous_indices), task%samples(1), task%samples(2)))
    task%kslice_data = cmplx_0

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = &
      SsTC_integer_array_element_to_memory_element(task, part_int_comp)

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Done."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_kslice_task_constructor

  subroutine SsTC_sample_kslice_task(task, system)

    implicit none

    class(SsTC_kslice_task), intent(inout) :: task
    type(SsTC_sys), intent(in)             :: system

    complex(kind=dp), allocatable :: data_k(:, :, :), &
                                     simd_tmp(:, :), &
                                     local_data_k(:, :, :)

    real(kind=dp) :: k(3)

    integer :: ik, k_ind(2), &
               ik1, ik2, &
               i, r, info

    real(kind=dp) :: start_time, end_time

    logical :: error = .false.

    type(SsTC_local_k_data) :: sampling_info

    integer, allocatable :: counts(:), displs(:)

    start_time = MPI_WTIME() !Start timer.

    allocate (sampling_info%integer_indices(2))
    sampling_info%integer_indices = task%samples

    allocate (counts(0:nProcs - 1), displs(0:nProcs - 1))
    call get_MPI_task_partition(product(task%samples), nProcs, counts, displs)

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Starting BZ sampling subroutine kslice."
    if (rank == 0) write (unit=stdout, fmt="(a)") "          Sampling task: "//trim(task%name)//&
    &" in the BZ for the system "//trim(system%name)//"."

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
    !_OMPTGT_(PRIVATE(ik, k_ind, ik1, ik2, k, simd_tmp))
    !_OMPTGT_(SIMD)
    do ik = displs(rank) + 1, displs(rank) + counts(rank)

      k_ind = SsTC_integer_memory_element_to_array_element(sampling_info, ik)
      ik1 = k_ind(1)
      ik2 = k_ind(2)

      k = task%corner + &
          task%vector(1, :)*real(ik1 - 1, dp)/real(task%samples(1) - 1, dp) + &
          task%vector(2, :)*real(ik2 - 1, dp)/real(task%samples(2) - 1, dp)

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
    & SsTC_sample_kslice_task. Check error log."

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
        call expand_array(data_k(:, i, r), task%kslice_data(i, r, :, :), info)
      enddo
    enddo

    end_time = MPI_WTIME() !End timer.

    deallocate (local_data_k, data_k)

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Sampling done."
    if (rank == 0) write (unit=stdout, fmt="(a, f15.3, a)") "          Total execution time: ", end_time - start_time, " s."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_sample_kslice_task

  subroutine SsTC_print_kslice(task, system)

    implicit none

    !Subroutine to format and output files related to the result of the task "task".

    class(SsTC_kslice_task), intent(in) :: task
    type(SsTC_sys), intent(in)          :: system

    real(kind=dp) :: k(3)

    character(len=400) :: filename, fmtf, num_label
    integer            :: i_arr(size(task%integer_indices)), &
                          r_arr(size(task%continuous_indices))
    integer            :: i_mem, r_mem, &
                          count, &
                          ik1, ik2
    integer            :: printunit

    if (rank == 0) write (unit=stdout, fmt="(a, a, a, a, a)") &
      "          Printing kslice task: "//trim(task%name)// &
      " for the system "//trim(system%name)//"."

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

        do ik1 = 1, task%samples(1)
          do ik2 = 1, task%samples(2)

            k = task%corner + &
                task%vector(1, :)*real(ik1 - 1, dp)/real(task%samples(1) - 1, dp) + &
                task%vector(2, :)*real(ik2 - 1, dp)/real(task%samples(2) - 1, dp)

            if (rank == 0) write (unit=printunit, fmt="(5e18.8e3)") k, real(task%kslice_data(i_mem, 1, ik1, ik2), dp), &
              aimag(task%kslice_data(i_mem, 1, ik1, ik2))

          enddo
          if (rank == 0) write (unit=printunit, fmt=*) ""
        enddo

        if (rank == 0) close (unit=printunit)

      enddo
    elseif (associated(task%global_calculator)) then

      write (fmtf, "(I10)") size(task%continuous_indices) + 5
      fmtf = '('//trim(adjustl(trim(fmtf)))//'E18.8E3)'

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

          do ik1 = 1, task%samples(1)
            do ik2 = 1, task%samples(2)

              k = task%corner + &
                  task%vector(1, :)*real(ik1 - 1, dp)/real(task%samples(1) - 1, dp) + &
                  task%vector(2, :)*real(ik2 - 1, dp)/real(task%samples(2) - 1, dp)

              if (rank == 0) write (unit=printunit, fmt=fmtf) k, &
                (task%ext_var_data(count)%data(r_arr(count)), count=1, size(task%continuous_indices)), &
                real(task%kslice_data(i_mem, r_mem, ik1, ik2), dp), aimag(task%kslice_data(i_mem, r_mem, ik1, ik2))

            enddo
            if (rank == 0) write (unit=printunit, fmt=*) ""
          enddo

        enddo

        if (rank == 0) close (unit=printunit)

      enddo
    endif

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Printing done."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_print_kslice

end module SsTC_kslice
