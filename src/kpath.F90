#include "cond_comp.h"
module SsTC_kpath

  USE OMP_LIB
  USE MPI_F08

  use SsTC_utility
  use SsTC_mpi_comms
  use SsTC_data_structures

  implicit none

  private

  type, extends(SsTC_global_k_data) :: SsTC_kpath_task
    !1st index is the id of the vector in the path.
    !2nd index corresponds to the component of the vector
    !in the path in coordinates relative to the reciprocal lattice.
    real(kind=dp), allocatable :: vectors(:, :)
    !Its size is the number of vectors in the BZ - 1.
    !number_of_pts(i) contains the number of k-points between vector i and vector i+1.
    integer, allocatable :: number_of_pts(:)
    !Array to store data with integer index,
    !continuous index and kpt index respectively.
    complex(kind=dp), allocatable :: kpath_data(:, :, :)
  end type SsTC_kpath_task

  public :: SsTC_kpath_task

  public :: SsTC_kpath_constructor
  public :: SsTC_kpath_sampler
  public :: SsTC_print_kpath

contains

  subroutine SsTC_kpath_constructor(task, &
                                    name, &
                                    l_calculator, g_calculator, &
                                    Nvec, vec_coord, nkpts, &
                                    N_int_ind, int_ind_range, &
                                    N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                                    part_int_comp)

    character(len=*) :: name

    procedure(SsTC_local_calculator), optional  :: l_calculator
    procedure(SsTC_global_calculator), optional :: g_calculator

    integer, intent(in)       :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in)       :: nkpts(Nvec - 1)

    integer, optional, intent(in) :: N_int_ind
    integer, optional, intent(in) :: int_ind_range(:)

    integer, optional, intent(in)       :: N_ext_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(:), ext_vars_end(:)
    integer, optional, intent(in)       :: ext_vars_steps(:)

    integer, optional, intent(in) :: part_int_comp(:)

    class(SsTC_kpath_task), intent(out) :: task

    integer :: i

    !Set name.
    task%name = name

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Creating kpath task "//trim(task%name)//"."

    !Set vector info
    allocate (task%vectors(Nvec, 3), task%number_of_pts(Nvec - 1))
    task%vectors = vec_coord
    task%number_of_pts = nkpts

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
    allocate (task%kpath_data(product(task%integer_indices), product(task%continuous_indices), &
                              sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2)))
    task%kpath_data = cmplx_0

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = &
      SsTC_integer_array_element_to_memory_element(task, part_int_comp)

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Done."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_kpath_constructor

  subroutine SsTC_kpath_sampler(task, system)

    class(SsTC_kpath_task), intent(inout) :: task
    type(SsTC_sys), intent(in)            :: system

    complex(kind=dp), allocatable :: data_k(:, :, :), &
                                     local_data_k(:, :, :)

    integer       :: ivec, isampling
    real(kind=dp) :: k(3)

    integer :: ik, &
               i, r

    integer :: TID
    real(kind=dp) :: start_time, end_time

    logical :: error = .false.

    integer, allocatable :: counts(:), displs(:)

    integer :: icount, count

    start_time = MPI_WTIME() !Start timer.

    allocate (counts(0:nProcs - 1), displs(0:nProcs - 1))
    call get_MPI_task_partition(sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2), nProcs, counts, displs)

    if (rank == 0) write (unit=stdout, fmt="(a)") &
      "          Sampling kpath task: "//trim(task%name)// &
      " for the system "//trim(system%name)//"."

    allocate (local_data_k(displs(rank) + 1:displs(rank) + counts(rank), &
                           product(task%integer_indices), product(task%continuous_indices)), &
              data_k(sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2), &
                     product(task%integer_indices), product(task%continuous_indices)))

    !_OMPTGT_(PARALLEL DEFAULT(SHARED) PRIVATE(TID, ik, count, icount, ivec, isampling, k))

    TID = OMP_GET_THREAD_NUM()
    IF ((TID .EQ. 0) .and. (rank .EQ. 0)) THEN
      write (unit=stdout, fmt="(a, i5, a, i5, a)") &
      &"         Running on ", nProcs, " process(es) each using ", OMP_GET_NUM_THREADS(), " threads."
    ENDIF

    !_OMPTGT_(DO)
    do ik = displs(rank) + 1, displs(rank) + counts(rank)

      !Retrieve ivec, the label corresponding to the ivec-th vector and isampling,
      !the discretization index corresponding to the vector starting at ivec-th vector
      !and ending in the ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
      if (ik == (sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2))) then
        ivec = size(task%vectors(:, 1)) - 1
        isampling = task%number_of_pts(ivec)
        goto 100
      endif
      count = 0
      do icount = 1, size(task%vectors(:, 1)) - 1
        count = count + task%number_of_pts(icount) - 1 !<--With this mapping we will get isampling int the range [1, task%number_of_pts(icount)).
        !However, since we are considering a task of size, sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2) (note the -2 instead of -1),
        !we will consider one last vector, thus the last ik corresponds to the ivec = size(task%vectors(:, 1))th vector.
        if (count .GE. ik) then
          ivec = icount
          isampling = ik - (count - (task%number_of_pts(icount) - 1))
          exit
        endif
      enddo
100   continue

      !Define a local vector from ivec-th vector to ivec+1-th vector discretized in task%number_of_pts(ivec) steps.

      if (task%number_of_pts(ivec) == 1) then
        k = task%vectors(ivec, :)
      else
        k = task%vectors(ivec, :) + &
            (task%vectors(ivec + 1, :) - task%vectors(ivec, :))*real(isampling - 1, dp) &
            /real(task%number_of_pts(ivec) - 1, dp)
      endif

      !Gather data.
      if (associated(task%local_calculator)) then
        local_data_k(ik, :, 1) = task%local_calculator(task, system, k, error)
      elseif (associated(task%global_calculator)) then
        local_data_k(ik, :, :) = task%global_calculator(task, system, k, error)
      endif

      if (error) then
        write (unit=stderr, fmt="(a, i5, a, 3e18.8e3)") "Rank ", rank, ": Error when sampling k-point", k
        write (unit=stderr, fmt="(a)") "Stopping..."
        stop
      endif

    enddo
    !_OMPTGT_(END DO)
    !_OMPTGT_(END PARALLEL)

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
        task%kpath_data(i, r, :) = data_k(:, i, r)
      enddo
    enddo

    deallocate (data_k, local_data_k)

    end_time = MPI_WTIME() !End timer.

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Sampling done."
    if (rank == 0) write (unit=stdout, fmt="(a, f15.3, a)") "          Total execution time: ", end_time - start_time, " s."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_kpath_sampler

  subroutine SsTC_print_kpath(task, system)
    !Subroutine to format and output files related to the result of the task "task".
    class(SsTC_kpath_task), intent(in) :: task
    type(SsTC_sys), intent(in)         :: system

    character(len=400) :: filename, fmtf
    integer            :: i_arr(size(task%integer_indices)), &
                          r_arr(size(task%continuous_indices))
    integer            :: i_mem, r_mem, &
                          count_int, &
                          ivec, isampling, &
                          count, icount, ik
    integer            :: printunit

    real(kind=dp) :: k(3)

    if (rank == 0) write (unit=stdout, fmt="(a)") &
      "          Printing kpath task: "//trim(task%name)// &
      " for the system "//trim(system%name)//"."

    if (associated(task%local_calculator)) then

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count_int = 1, size(task%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count_int))
        enddo
        filename = trim(filename)//'.dat'

        if (rank == 0) open (newunit=printunit, action="write", file=filename)

        do ik = 1, sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2)

          !Retrieve ivec, the label corresponding to the ivec-th vector and isampling,
          !the discretization index corresponding to the vector starting at ivec-th vector
          !and ending in the ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
          if (ik == (sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2))) then
            ivec = size(task%vectors(:, 1)) - 1
            isampling = task%number_of_pts(ivec)
            goto 200
          endif
          count = 0
          do icount = 1, size(task%vectors(:, 1)) - 1
            count = count + task%number_of_pts(icount) - 1 !<--With this mapping we will get isampling int the range [1, task%number_of_pts(icount)).
            !However, since we are considering a task of size, sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2) (note the -2 instead of -1),
            !we will consider one last vector, thus the last ik corresponds to the ivec = size(task%vectors(:, 1))th vector.
            if (count .GE. ik) then
              ivec = icount
              isampling = ik - (count - (task%number_of_pts(icount) - 1))
              exit
            endif
          enddo
200       continue

          !Define a local vector from ivec-th vector to ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
          if (task%number_of_pts(ivec) == 1) then
            k = task%vectors(ivec, :)
          else
            k = task%vectors(ivec, :) + &
                (task%vectors(ivec + 1, :) - task%vectors(ivec, :))*real(isampling - 1, dp) &
                /real(task%number_of_pts(ivec) - 1, dp)
          endif

          if (rank == 0) write (unit=printunit, fmt="(6E18.8E3)") real(ik, dp), k, real(task%kpath_data(i_mem, 1, ik), dp), &
            aimag(task%kpath_data(i_mem, 1, ik))

        enddo

        if (rank == 0) close (unit=printunit)

      enddo

    elseif (associated(task%global_calculator)) then

      write (fmtf, "(I10)") size(task%continuous_indices) + 6
      fmtf = '('//trim(adjustl(trim(fmtf)))//'E18.8E3)'

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count_int = 1, size(task%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count_int))
        enddo
        filename = trim(filename)//'.dat'

        if (rank == 0) open (newunit=printunit, action="write", file=filename)

        do r_mem = 1, product(task%continuous_indices) !For each continuous index.

          r_arr = SsTC_continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.

          do ik = 1, sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2)

            !Retrieve ivec, the label corresponding to the ivec-th vector and isampling,
            !the discretization index corresponding to the vector starting at ivec-th vector
            !and ending in the ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
            if (ik == (sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2))) then
              ivec = size(task%vectors(:, 1)) - 1
              isampling = task%number_of_pts(ivec)
              goto 300
            endif
            count = 0
            do icount = 1, size(task%vectors(:, 1)) - 1
              count = count + task%number_of_pts(icount) - 1 !<--With this mapping we will get isampling int the range [1, task%number_of_pts(icount)).
              !However, since we are considering a task of size, sum(task%number_of_pts) - (size(task%vectors(:, 1)) - 2) (note the -2 instead of -1),
              !we will consider one last vector, thus the last ik corresponds to the ivec = size(task%vectors(:, 1))th vector.
              if (count .GE. ik) then
                ivec = icount
                isampling = ik - (count - (task%number_of_pts(icount) - 1))
                exit
              endif
            enddo
300         continue

            !Define a local vector from ivec-th vector to ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
            if (task%number_of_pts(ivec) == 1) then
              k = task%vectors(ivec, :)
            else
              k = task%vectors(ivec, :) + &
                  (task%vectors(ivec + 1, :) - task%vectors(ivec, :))*real(isampling - 1, dp)/ &
                  real(task%number_of_pts(ivec) - 1, dp)
            endif

            if (rank == 0) write (unit=printunit, fmt=fmtf) real(ik, dp), k, &
              (task%ext_var_data(count_int)%data(r_arr(count_int)), count_int=1, size(task%continuous_indices)), &
              real(task%kpath_data(i_mem, r_mem, ik), dp), aimag(task%kpath_data(i_mem, r_mem, ik))

          enddo

        enddo

        if (rank == 0) close (unit=printunit)

      enddo

    endif

    if (rank == 0) write (unit=stdout, fmt="(a)") "          Printing done."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_print_kpath

end module SsTC_kpath
