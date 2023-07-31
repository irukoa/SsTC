module kpath

  USE OMP_LIB

  use utility
  use data_structures

  implicit none

  type, extends(global_k_data) :: k_path_task
    !1st index is the id of the vector in the path.
    !2nd index corresponds to the component of the vector
    !in the path in coordinates relative to the reciprocal lattice.
    real(kind=dp), allocatable :: vectors(:, :)
    !Its size is the number of vectors in the BZ - 1.
    !vector(i) contains the number of k-points between vector i and vector i+1.
    integer, allocatable :: number_of_pts(:)
    !Array to store data with integer index,
    !continuous index and kpt index respectively.
    complex(kind=dp), allocatable :: kpath_data(:, :, :)
  end type k_path_task

  public :: kpath_constructor
  public :: kpath_sampler
  public :: print_kpath

contains

  subroutine kpath_constructor(task, &
                               name, &
                               l_calculator, g_calculator, &
                               Nvec, vec_coord, nkpts, &
                               N_int_ind, int_ind_range, &
                               N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                               part_int_comp)
    character(len=*) :: name

    procedure(local_calculator), optional :: l_calculator
    procedure(global_calculator), optional :: g_calculator

    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)

    integer, optional, intent(in) :: N_int_ind
    integer, optional, intent(in) :: int_ind_range(:)

    integer, optional, intent(in) :: N_ext_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(:), ext_vars_end(:)
    integer, optional, intent(in) :: ext_vars_steps(:)

    integer, optional, intent(in) :: part_int_comp(:)

    class(k_path_task), intent(out) :: task

    integer :: i

    !Set name.
    task%name = name

    write (unit=stdout, fmt="(a)") "Creating kpath task "//trim(task%name)//"."

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
        task%ext_var_data(i) = external_variable_constructor(ext_vars_start(i), ext_vars_end(i), ext_vars_steps(i))
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
    allocate (task%kpath_data(product(task%integer_indices), product(task%continuous_indices), sum(task%number_of_pts)))
    task%kpath_data = cmplx_0

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = integer_array_element_to_memory_element(task, part_int_comp)

    write (unit=stdout, fmt="(a)") "Done."
    write (unit=stdout, fmt="(a)") ""

  end subroutine kpath_constructor

  subroutine kpath_sampler(task, system)

    class(k_path_task), intent(inout) :: task
    type(sys), intent(in)    :: system

    integer       :: ivec, isampling
    real(kind=dp) :: k(3)

    complex(kind=dp), allocatable :: temp_res(:, :, :)

    logical :: error
    error = .false.

    allocate (temp_res(product(task%integer_indices), product(task%continuous_indices), sum(task%number_of_pts)))

    write (unit=stdout, fmt="(a)") "Sampling kpath task: "//trim(task%name)//" for the system "//trim(system%name)//"."

    !Sampling.
    do ivec = 1, size(task%vectors(:, 1)) - 1 !For each considered vector except the last one.
!$OMP       PARALLEL SHARED(temp_res) PRIVATE(isampling, k)
!$OMP       DO
      do isampling = 1, task%number_of_pts(ivec)
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
          temp_res(:, 1, sum(task%number_of_pts(1:ivec - 1)) + isampling) = task%local_calculator(task, system, k, error)
        elseif (associated(task%global_calculator)) then
          temp_res(:, :, sum(task%number_of_pts(1:ivec - 1)) + isampling) = task%global_calculator(task, system, k, error)
        endif

        if (error) then
          write (unit=stderr, fmt="(a, 3e18.8e3)") "Error when sampling k-point", k
          write (unit=stderr, fmt="(a)") "Stopping..."
          stop
        endif

      enddo
!$OMP       END DO
!$OMP       END PARALLEL
    enddo

    write (unit=stdout, fmt="(a)") "Sampling done."
    write (unit=stdout, fmt="(a)") ""

    task%kpath_data = temp_res

    deallocate (temp_res)

  end subroutine kpath_sampler

  subroutine print_kpath(task, system)
    !Subroutine to format and output files related to the result of the task "task".
    class(k_path_task), intent(in) :: task
    type(sys), intent(in) :: system

    character(len=400) :: filename, fmtf
    integer :: i_arr(size(task%integer_indices)), r_arr(size(task%continuous_indices))
    integer :: i_mem, r_mem, count, countk, ivec, isampling
    integer :: printunit
    real(kind=dp) :: k(3)

    write (unit=stdout, fmt="(a)") "Printing kpath task: "//trim(task%name)//" for the system "//trim(system%name)//"."

    if (associated(task%local_calculator)) then

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count = 1, size(task%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count))
        enddo
        filename = trim(filename)//'.dat'

        open (newunit=printunit, action="write", file=filename)

        countk = 0
        do ivec = 1, size(task%vectors(:, 1)) - 1 !For each considered vector except the last one.
          do isampling = 1, task%number_of_pts(ivec)
            countk = countk + 1
            !Define a local vector from ivec-th vector to ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
            if (task%number_of_pts(ivec) == 1) then
              k = task%vectors(ivec, :)
            else
              k = task%vectors(ivec, :) + &
                  (task%vectors(ivec + 1, :) - task%vectors(ivec, :))*real(isampling - 1, dp) &
                  /real(task%number_of_pts(ivec) - 1, dp)
            endif
            write (unit=printunit, fmt="(6E18.8E3)") real(countk, dp), k, real(task%kpath_data(i_mem, 1, countk), dp), &
              aimag(task%kpath_data(i_mem, 1, countk))
          enddo
        enddo

        close (unit=printunit)

      enddo
    elseif (associated(task%global_calculator)) then

      write (fmtf, "(I10)") size(task%continuous_indices) + 6
      fmtf = '('//trim(adjustl(trim(fmtf)))//'E18.8E3)'

      do i_mem = 1, product(task%integer_indices) !For each integer index.

        i_arr = integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(task%name)//'_'
        do count = 1, size(task%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count))
        enddo
        filename = trim(filename)//'.dat'

        open (newunit=printunit, action="write", file=filename)

        do r_mem = 1, product(task%continuous_indices) !For each continuous index.

          r_arr = continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.

          countk = 0
          do ivec = 1, size(task%vectors(:, 1)) - 1 !For each considered vector except the last one.
            do isampling = 1, task%number_of_pts(ivec)
              countk = countk + 1
              !Define a local vector from ivec-th vector to ivec+1-th vector discretized in task%number_of_pts(ivec) steps.
              if (task%number_of_pts(ivec) == 1) then
                k = task%vectors(ivec, :)
              else
                k = task%vectors(ivec, :) + &
                    (task%vectors(ivec + 1, :) - task%vectors(ivec, :))*real(isampling - 1, dp)/ &
                    real(task%number_of_pts(ivec) - 1, dp)
              endif

              write (unit=printunit, fmt=fmtf) real(countk, dp), k, &
                (task%ext_var_data(count)%data(r_arr(count)), count=1, size(task%continuous_indices)), &
                real(task%kpath_data(i_mem, r_mem, countk), dp), aimag(task%kpath_data(i_mem, r_mem, countk))
            enddo
          enddo

        enddo

        close (unit=printunit)

      enddo
    endif

    write (unit=stdout, fmt="(a)") "Printing done."
    write (unit=stdout, fmt="(a)") ""

  end subroutine print_kpath

end module kpath
