module kpath

  USE OMP_LIB

  use utility
  use data_structures

  implicit none

  type, extends(global_k_data) :: k_path
    real(kind=dp), allocatable :: vectors(:, :) !1st index is the index of the vector in the path. 2nd index corresponds to the component of the vector in the path, so its size is vectors(size(number_of_vecotrs), 3).
    integer,       allocatable :: number_of_pts(:) !Its size is the number of vectors in the BZ, vector(i) contains the number of k-points between vector i and vector i+1.
    complex(kind=dp), allocatable :: kpath_data(:, :, :)!Integer index, continuous index and kpt index respectively.
  end type k_path

  contains

  function kpath_constructor(name, &
                             l_calculator, g_calculator, &
                             Nvec, vec_coord, nkpts, &
                             N_int_ind, int_ind_range, &
                             N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                             part_int_comp) result(path)

    character(len=*) :: name

    procedure (local_calculator) :: l_calculator
    procedure (global_calculator), optional :: g_calculator

    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)

    integer, intent(in) :: N_int_ind
    integer, intent(in) :: int_ind_range(N_int_ind)

    integer, intent(in)                     :: N_ext_vars
    real(kind=dp),    optional,  intent(in) :: ext_vars_start(N_ext_vars), ext_vars_end(N_ext_vars)
    integer,          optional,  intent(in) :: ext_vars_steps(N_ext_vars)

    integer,          optional,  intent(in) :: part_int_comp(N_int_ind)

    type(k_path) :: path

    integer :: i

    !Set name.
    path%name = name

    !Set vector info
    allocate(path%vectors(Nvec, 3), path%number_of_pts(Nvec-1))
    path%vectors = vec_coord
    path%number_of_pts = nkpts

    !Set integer index data.
    allocate(path%integer_indices(N_int_ind))
    do i = 1, N_int_ind
      path%integer_indices(i) = int_ind_range(i)
    enddo

    !Set external variable data.
    if (((N_ext_vars).ge.1).and.(present(ext_vars_start)).and.(present(ext_vars_end)).and.(present(ext_vars_steps))) then
      allocate(path%continuous_indices(N_ext_vars), path%ext_var_data(N_ext_vars))
      do i = 1, N_ext_vars
        path%continuous_indices(i) = ext_vars_steps(i)
        allocate(path%ext_var_data(i)%data(ext_vars_steps(i)))
        path%ext_var_data(i) = external_variable_constructor(ext_vars_start(i), ext_vars_end(i), ext_vars_steps(i))
      enddo
    else
      allocate(path%continuous_indices(1), path%ext_var_data(1))
      path%continuous_indices(1) = 1.0_dp
      allocate(path%ext_var_data(1)%data(1))
      path%ext_var_data(1)%data = (/1.0_dp/)
    endif

    if (present(g_calculator)) then
      !Set calculator pointer (function alias).
      path%global_calculator => g_calculator
      nullify(path%local_calculator)
    else
      !Set calculator pointer (function alias).
      path%local_calculator => l_calculator
      nullify(path%global_calculator)
    endif

    !Set kdata.
    allocate(path%kpath_data(product(path%integer_indices), product(path%continuous_indices), sum(path%number_of_pts)))
    path%kpath_data = cmplx_0

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) path%particular_integer_component = integer_array_element_to_memory_element(path, part_int_comp)

  end function kpath_constructor

  subroutine kpath_sampler(path, system)

    type(k_path), intent(inout) :: path
    type(sys), intent(in)    :: system

    integer :: ivec, isampling
    real(kind=dp) :: k(3)

    complex(kind=dp), allocatable :: temp_res(:, :, :)

    allocate(temp_res(product(path%integer_indices), product(path%continuous_indices), sum(path%number_of_pts)))

    !$OMP PARALLEL SHARED(temp_res) PRIVATE(ivec, isampling, k)
    !$OMP DO
    do ivec = 1, size(path%vectors(:, 1)) - 1 !For each considered vector except the last one.
      do isampling = 1, path%number_of_pts(ivec)
        !Define a local vector from ivec-th vector to ivec+1-th vector discretized in path%number_of_pts(ivec) steps.
        k = path%vectors(ivec, :) + (path%vectors(ivec + 1, :) - path%vectors(ivec, :))*real(isampling - 1, dp)/real(path%number_of_pts(ivec) - 1, dp)
        !Gather data.
        if (associated(path%local_calculator)) then
          temp_res(:, 1, sum(path%number_of_pts(1:ivec-1)) + isampling) = path%local_calculator(path, system, k)
        else
          temp_res(:, :, sum(path%number_of_pts(1:ivec-1)) + isampling) = path%global_calculator(path, system, k)
        endif
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    path%kpath_data = temp_res
    deallocate(temp_res)

  end subroutine kpath_sampler

  subroutine print_kpath(path, system)
    !Subroutine to format and output files related to the result of the path "path".
    type(k_path), intent(in) :: path
    type(sys), intent(in) :: system

    character(len=400) :: filename
    integer :: i_arr(size(path%integer_indices)), r_arr(size(path%continuous_indices))
    integer :: i_mem, r_mem, count, countk, ivec, isampling
    real(kind=dp) :: k(3)

    if (associated(path%local_calculator)) then
      do i_mem = 1, product(path%integer_indices) !For each integer index.

        i_arr = integer_memory_element_to_array_element(path, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(path%name)//'_'
        do count = 1, size(path%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count))
        enddo
        filename = trim(filename)//'.dat'

        open(unit=111, action="write", file=filename)

        countk = 0
        do ivec = 1, size(path%vectors(:, 1)) - 1 !For each considered vector except the last one.
          do isampling = 1, path%number_of_pts(ivec)
            countk = countk + 1
            !Define a local vector from ivec-th vector to ivec+1-th vector discretized in path%number_of_pts(ivec) steps.
            k = path%vectors(ivec, :) + (path%vectors(ivec + 1, :) - path%vectors(ivec, :))*real(isampling - 1, dp)/real(path%number_of_pts(ivec) - 1, dp)
            write(unit=111, fmt="(6E18.8E3)") real(countk, dp), k, real(path%kpath_data(i_mem, 1, countk), dp), aimag(path%kpath_data(i_mem, 1, countk))
          enddo
        enddo

        close(unit=111)

      enddo
    else!global calculator is associated.
      do i_mem = 1, product(path%integer_indices) !For each integer index.

        i_arr = integer_memory_element_to_array_element(path, i_mem) !Pass to array layout.

        filename = trim(system%name)//'-'//trim(path%name)//'_'
        do count = 1, size(path%integer_indices)
          filename = trim(filename)//achar(48 + i_arr(count))
        enddo
        filename = trim(filename)//'.dat'

        open(unit=111, action="write", file=filename)

        do r_mem = 1, product(path%continuous_indices) !For each continuous index.

          r_arr = continuous_memory_element_to_array_element(path, r_mem) !Pass to array layout.

          countk = 0
          do ivec = 1, size(path%vectors(:, 1)) - 1 !For each considered vector except the last one.
            do isampling = 1, path%number_of_pts(ivec)
              countk = countk + 1
              !Define a local vector from ivec-th vector to ivec+1-th vector discretized in path%number_of_pts(ivec) steps.
              k = path%vectors(ivec, :) + (path%vectors(ivec + 1, :) - path%vectors(ivec, :))*real(isampling - 1, dp)/real(path%number_of_pts(ivec) - 1, dp)

              write(unit=111, fmt=*) real(countk, dp), k, (path%ext_var_data(count)%data(r_arr(count)), count = 1, size(path%continuous_indices)), &
              real(path%kpath_data(i_mem, r_mem, countk), dp), aimag(path%kpath_data(i_mem, r_mem, countk))
            enddo
          enddo
  
        enddo

        close(unit=111)

      enddo
    endif
  end subroutine print_kpath

end module kpath