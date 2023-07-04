module kpath

  USE OMP_LIB

  use utility
  use data_structures

  implicit none

  type, extends(local_k_data) :: k_path !TODO: MAKE CONSTRUCTOR.
    real(kind=dp), allocatable :: vectors(:, :) !1st index is the index of the vector in the path. 2nd index corresponds to the component of the vector in the path, so its size is vectors(size(number_of_vecotrs), 3).
    integer,       allocatable :: number_of_pts(:) !Its size is the number of vectors in the BZ, vector(i) contains the number of k-points between vector i and vector i+1.
    complex(kind=dp), allocatable :: kpath_data(:, :)
  end type k_path

  contains

  function kpath_constructor(name, calculator, Nvec, vec_coord, nkpts, N_int_ind, int_ind_range) result(path)
    character(len=*) :: name
    procedure (local_calculator) :: calculator
    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)
    integer, intent(in) :: N_int_ind
    integer, intent(in) :: int_ind_range(N_int_ind)

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

    !Set calculator pointer (function alias).
    path%local_calculator => calculator

    !Set kdata.
    allocate(path%kpath_data(product(path%integer_indices), sum(path%number_of_pts)))
    path%kpath_data = cmplx_0

  end function kpath_constructor

  subroutine kpath_sampler(path, system)

    type(k_path), intent(inout) :: path
    type(sys), intent(in)    :: system

    integer :: ivec, isampling
    real(kind=dp) :: k(3)

    complex(kind=dp), allocatable :: temp_res(:, :)

    allocate(temp_res(product(path%integer_indices), sum(path%number_of_pts)))

    !$OMP PARALLEL SHARED(temp_res) PRIVATE(ivec, isampling, k)
    !$OMP DO
    do ivec = 1, size(path%vectors(:, 1)) - 1 !For each considered vector except the last one.
      do isampling = 1, path%number_of_pts(ivec)
        !Define a local vector from ivec-th vector to ivec+1-th vector discretized in path%number_of_pts(ivec) steps.
        k = path%vectors(ivec, :) + (path%vectors(ivec + 1, :) - path%vectors(ivec, :))*real(isampling - 1, dp)/real(path%number_of_pts(ivec) - 1, dp)
        !Gather data.
        temp_res(:, sum(path%number_of_pts(1:ivec-1)) + isampling) = path%local_calculator(path, system, k)
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
    integer :: i_arr(size(path%integer_indices))
    integer :: i_mem, count, ivec, isampling
    real(kind=dp) :: k(3)

    do i_mem = 1, product(path%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(path, i_mem) !Pass to array layout.

      filename = trim(system%name)//'-'//trim(path%name)//'_'
      do count = 1, size(path%integer_indices)
        filename = trim(filename)//achar(48 + i_arr(count))
      enddo
      filename = trim(filename)//'.dat'

      open(unit=111, action="write", file=filename)

      count = 0
      do ivec = 1, size(path%vectors(:, 1)) - 1 !For each considered vector except the last one.
        do isampling = 1, path%number_of_pts(ivec)
          count = count + 1
          !Define a local vector from ivec-th vector to ivec+1-th vector discretized in path%number_of_pts(ivec) steps.
          k = path%vectors(ivec, :) + (path%vectors(ivec + 1, :) - path%vectors(ivec, :))*real(isampling - 1, dp)/real(path%number_of_pts(ivec) - 1, dp)
          write(unit=111, fmt="(6E18.8E3)") real(count, dp), k, real(path%kpath_data(i_mem, count), dp), aimag(path%kpath_data(i_mem, count))
        enddo
      enddo

      close(unit=111)

    enddo

  end subroutine print_kpath

end module kpath