module SsTC_data_structures

  use SsTC_utility
  use SsTC_mpi_comms

  implicit none

  private

  type SsTC_sys
    character(len=120)            :: name
    integer                       :: num_bands                                !Number of bands.
    real(kind=dp)                 :: direct_lattice_basis(3, 3)               !Direct lattice basis vectors in A.
    !1st index is vector label, 2nd index is vector component.
    real(kind=dp)                 :: metric_tensor(3, 3)                      !Metric tensor of the direct lattice basis.
    real(kind=dp)                 :: cell_volume                              !Cell volume.
    integer                       :: num_R_points                             !Number of R points given in the *_tb.dat file.
    integer, allocatable          :: R_point(:, :)                            !Memory layout id of the R-point (1st index) and
    !R-vector coords. relative to the direct lattice
    !basis vectors (2nd index).
    integer, allocatable          :: deg_R_point(:)                           !Degeneracy of the R-point specified by its
    !memory layout id.
    complex(kind=dp), allocatable :: real_space_hamiltonian_elements(:, :, :) !Hamiltonian matrix elements (1st and 2nd indexes)
    !and memory layout id of the R-point (3rd index) in eV.
    complex(kind=dp), allocatable :: real_space_position_elements(:, :, :, :) !Position operator matrix elements
    !(1st and 2nd indexes),
    !cartesian coordinate (3rd index) and memory layout id
    !of the R-point (4th index) in A.
    real(kind=dp)                 :: e_fermi = 0.0_dp                         !Fermi energy.
    real(kind=dp)                 :: deg_thr = 1.0E-4_dp                      !Degeneracy threshold in eV.
    real(kind=dp)                 :: deg_offset = 0.04_dp                     !Offset for regularization in case of
    !degeneracies, in eV.
  end type SsTC_sys

  type SsTC_external_vars
    real(kind=dp), allocatable :: data(:) !External variable data array.
  end type SsTC_external_vars

  type SsTC_local_k_data
    character(len=120)                                :: name
    integer, allocatable                              :: integer_indices(:)               !Each entry contains the range
    !of each of the integer indices.
    complex(kind=dp), allocatable                     :: k_data(:)                        !Data local for each k with
    !integer index in memory array,
    procedure(SsTC_local_calculator), pointer, nopass :: local_calculator => null()       !Pointer to the local calculator.
    integer                                           :: particular_integer_component = 0 !Specification of some integer component.
  end type SsTC_local_k_data

  type, extends(SsTC_local_k_data) :: SsTC_global_k_data
    integer, allocatable                               :: continuous_indices(:)           !Each entry contains the range
    !of each continuous indices.
    type(SsTC_external_vars), allocatable              :: ext_var_data(:)                 !External variable data.
    procedure(SsTC_global_calculator), pointer, nopass :: global_calculator => null()     !Pointer to the global calculator.
    integer, allocatable                               :: iterables(:, :)                 !Iterable dictionary.
  end type SsTC_global_k_data

  !Interfaces for the generic functions returning k dependent quantities.
  abstract interface !nvfortran: remove abstract to avoid error.
    function SsTC_local_calculator(k_data, system, k, error) result(u)
      import :: SsTC_local_k_data, SsTC_sys, SsTC_external_vars, dp

      class(SsTC_local_k_data), intent(in) :: k_data
      type(SsTC_sys), intent(in)           :: system
      real(kind=dp), intent(in)            :: k(3)
      logical, intent(inout)               :: error

      complex(kind=dp) :: u(product(k_data%integer_indices))
    end function SsTC_local_calculator
  end interface

  abstract interface !nvfortran: remove abstract to avoid error.
    function SsTC_global_calculator(task, system, k, error) result(u)
      import :: SsTC_global_k_data, SsTC_sys, SsTC_external_vars, dp

      class(SsTC_global_k_data), intent(in) :: task
      type(SsTC_sys), intent(in)            :: system
      real(kind=dp), intent(in)             :: k(3)
      logical, intent(inout)                :: error

      complex(kind=dp) :: u(product(task%integer_indices), product(task%continuous_indices))
    end function SsTC_global_calculator
  end interface

  public :: SsTC_local_k_data
  public :: SsTC_global_k_data

  public :: SsTC_local_calculator
  public :: SsTC_global_calculator

  public :: SsTC_sys
  public :: SsTC_sys_constructor

  public :: SsTC_external_vars
  public :: SsTC_external_variable_constructor

  public :: SsTC_integer_array_element_to_memory_element
  public :: SsTC_integer_memory_element_to_array_element
  public :: SsTC_continuous_array_element_to_memory_element
  public :: SsTC_continuous_memory_element_to_array_element
  public :: SsTC_construct_iterable

contains

  function SsTC_sys_constructor(name, path_to_tb_file, &
                                efermi, deg_thr, deg_offset) &
    result(system)

    implicit none

    character(len=*), intent(in)        :: name
    character(len=*), intent(in)        :: path_to_tb_file
    real(kind=dp), optional, intent(in) :: efermi, &
                                           deg_thr, deg_offset

    type(SsTC_sys) :: system

    character(len=400) :: filename
    integer            :: num_bands, nrpts, &
                          division, remainder, &
                          irpts, i, j, &
                          dummy1, dummy2, dummy3

    real(kind=dp), allocatable :: dummyR(:)

    system%name = name

    if (present(efermi)) system%e_fermi = efermi
    if (present(deg_thr)) system%deg_thr = deg_thr
    if (present(deg_offset)) system%deg_offset = deg_offset

    filename = trim(path_to_tb_file)//trim(name)//"_tb.dat"
    filename = trim(filename)

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a, a, a)") "          Initializing system "//trim(name)//"."
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a, a, a)") "          Reading file "//trim(filename)//"."

    open (newunit=stdin, action="read", file=filename)
    read (unit=stdin, fmt=*)
    do i = 1, 3
      read (unit=stdin, fmt=*) (system%direct_lattice_basis(i, j), j=1, 3)
    enddo

    do i = 1, 3
      do j = 1, 3
        system%metric_tensor(i, j) = dot_product(system%direct_lattice_basis(i, :), system%direct_lattice_basis(j, :))
      enddo
    enddo
    system%cell_volume = sqrt(system%metric_tensor(1, 1)*system%metric_tensor(2, 2)*system%metric_tensor(3, 3) + &
                              system%metric_tensor(2, 1)*system%metric_tensor(3, 2)*system%metric_tensor(1, 3) + &
                              system%metric_tensor(1, 2)*system%metric_tensor(2, 3)*system%metric_tensor(3, 1) - &
                              system%metric_tensor(3, 1)*system%metric_tensor(2, 2)*system%metric_tensor(1, 3) - &
                              system%metric_tensor(2, 1)*system%metric_tensor(1, 2)*system%metric_tensor(3, 3) - &
                              system%metric_tensor(3, 2)*system%metric_tensor(2, 3)*system%metric_tensor(1, 1))

    read (unit=stdin, fmt=*) num_bands
    system%num_bands = num_bands
    read (unit=stdin, fmt=*) nrpts
    system%num_R_points = nrpts
    allocate (system%R_point(system%num_R_points, 3), &
              system%deg_R_point(system%num_R_points), &
              system%real_space_hamiltonian_elements(system%num_bands, system%num_bands, system%num_R_points), &
              system%real_space_position_elements(system%num_bands, system%num_bands, 3, system%num_R_points))

    division = nrpts/15
    remainder = modulo(nrpts, 15)

    if (remainder == 0) then
      do i = 1, division
        read (unit=stdin, fmt=*) (system%deg_R_point(15*(i - 1) + j), j=1, 15)
      enddo
    else
      do i = 1, division
        read (unit=stdin, fmt=*) (system%deg_R_point(15*(i - 1) + j), j=1, 15)
      enddo
      read (unit=stdin, fmt=*) (system%deg_R_point(15*(i - 1) + j), j=1, remainder)
    endif

    read (unit=stdin, fmt=*)
    allocate (dummyR(2))
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Reading Hamiltonian..."

    system%real_space_hamiltonian_elements = cmplx_0
    do irpts = 1, nrpts
      read (unit=stdin, fmt=*) (system%R_point(irpts, i), i=1, 3)
      do i = 1, num_bands
        do j = 1, num_bands
          read (unit=stdin, fmt=*) dummy1, dummy2, dummyR(1), dummyR(2)
          !As pointed out in the W90 source code v3.1.0, get_oper.F90, ln 147, addition is required instead of equality.
          system%real_space_hamiltonian_elements(i, j, irpts) = &
            system%real_space_hamiltonian_elements(i, j, irpts) + cmplx(dummyR(1), dummyR(2), dp)
        enddo
      enddo
      read (unit=stdin, fmt=*)
    enddo

    deallocate (dummyR)
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Done."

    allocate (dummyR(6))
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Reading position operator..."

    system%real_space_position_elements = cmplx_0
    do irpts = 1, nrpts
      read (unit=stdin, fmt=*) dummy1, dummy2, dummy3
      do i = 1, num_bands
        do j = 1, num_bands
          read (unit=stdin, fmt=*) dummy1, dummy2, dummyR(1), dummyR(2), dummyR(3), dummyR(4), dummyR(5), dummyR(6)
          !As pointed out in the W90 source code v3.1.0, get_oper.F90, ln 147, addition is required instead of equality.
          system%real_space_position_elements(i, j, 1, irpts) = &
            system%real_space_position_elements(i, j, 1, irpts) + cmplx(dummyR(1), dummyR(2), dp)
          system%real_space_position_elements(i, j, 2, irpts) = &
            system%real_space_position_elements(i, j, 2, irpts) + cmplx(dummyR(3), dummyR(4), dp)
          system%real_space_position_elements(i, j, 3, irpts) = &
            system%real_space_position_elements(i, j, 3, irpts) + cmplx(dummyR(5), dummyR(6), dp)
        enddo
      enddo
      if (irpts == nrpts) cycle
      read (unit=stdin, fmt=*)
    enddo

    deallocate (dummyR)
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          Done."

    close (unit=stdin)

    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") "          System loaded sucessfully."
    if ((rank == 0) .and. verbose) write (unit=stdout, fmt="(a)") ""

  end function SsTC_sys_constructor

  function SsTC_external_variable_constructor(start, end, steps) result(vars)

    implicit none

    !Function to set external variable data.
    real(kind=dp), intent(in) :: start, end
    integer, intent(in)       :: steps

    type(SsTC_external_vars) :: vars

    integer :: i

    allocate (vars%data(steps))

    if (steps == 1) then
      vars%data(1) = start
    else
      do i = 1, steps
        vars%data(i) = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
      enddo
    endif

  end function SsTC_external_variable_constructor

  !==UTILITY FUNCTIONS TO PASS "MEMORY LAYOUT" INDICES TO "ARRAY LAYOUT" AND VICE-VERSA==!
  !See discussion at https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays

  function SsTC_integer_array_element_to_memory_element(data_k, i_arr) result(i_mem)

    implicit none

    !Get integer indices from array layout to memory layout.
    class(SsTC_local_k_data), intent(in) :: data_k
    integer, intent(in)                  :: i_arr(size(data_k%integer_indices))

    integer :: i_mem
    integer :: counter

    i_mem = 1
    do counter = 0, size(data_k%integer_indices) - 1
      i_mem = (i_mem - 1)*data_k%integer_indices(counter + 1) + i_arr(counter + 1)
    enddo

  end function SsTC_integer_array_element_to_memory_element

  function SsTC_integer_memory_element_to_array_element(data_k, i_mem) result(i_arr)

    implicit none

    !Get integer indices from memory layout to array layout.
    class(SsTC_local_k_data), intent(in) :: data_k
    integer, intent(in)                  :: i_mem

    integer :: i_arr(size(data_k%integer_indices))
    integer :: counter, reduction

    reduction = i_mem
    do counter = size(data_k%integer_indices), 1, -1
      if (counter == 1) then
        i_arr(counter) = reduction
      else
        i_arr(counter) = modulo(reduction, data_k%integer_indices(counter))
        if (i_arr(counter) == 0) i_arr(counter) = data_k%integer_indices(counter)
        reduction = int((reduction - i_arr(counter))/data_k%integer_indices(counter)) + 1
      endif
    enddo

  end function SsTC_integer_memory_element_to_array_element

  function SsTC_continuous_array_element_to_memory_element(task, r_arr) result(r_mem)

    implicit none

    !Get continuous indices from array layout to memory layout.
    class(SsTC_global_k_data), intent(in) :: task
    integer, intent(in)                   :: r_arr(size(task%continuous_indices))

    integer :: r_mem
    integer :: counter

    r_mem = 1
    do counter = 0, size(task%continuous_indices) - 1
      r_mem = (r_mem - 1)*task%continuous_indices(counter + 1) + r_arr(counter + 1)
    enddo

  end function SsTC_continuous_array_element_to_memory_element

  function SsTC_continuous_memory_element_to_array_element(task, r_mem) result(r_arr)

    implicit none

    !Get continuous indices from memory layout to array layout.
    class(SsTC_global_k_data), intent(in) :: task
    integer, intent(in)                   :: r_mem

    integer :: r_arr(size(task%continuous_indices))
    integer :: counter, reduction

    reduction = r_mem
    do counter = size(task%continuous_indices), 1, -1
      if (counter == 1) then
        r_arr(counter) = reduction
      else
        r_arr(counter) = modulo(reduction, task%continuous_indices(counter))
        if (r_arr(counter) == 0) r_arr(counter) = task%continuous_indices(counter)
        reduction = int((reduction - r_arr(counter))/task%continuous_indices(counter)) + 1
      endif
    enddo

  end function SsTC_continuous_memory_element_to_array_element

  subroutine SsTC_construct_iterable(global, vars)

    implicit none

    !Creates a dictionaty with all the possible permutations of the
    !considered variation of the continuous variables specified
    !in the array elements of "vars".
    !E.g., if global has N continuous variables with size s_i for
    !i = [1, N] and vars is a 1d array of size j < N,
    !containing in each entry the label of the variable
    !that will be permutated i_j, the iterable
    !part of global will be set with a 2d dictionaty containing, in the
    !1st index the label of the permutation and in the 2nd index
    !the particular permutation, of size(global%continuous_indices),
    !of the variables i_j for all j in size(vars).
    class(SsTC_global_k_data), intent(inout) :: global
    integer, intent(in)                      :: vars(:)

    integer :: n, &
               i, j, &
               counter, reduction

    integer, allocatable :: temp_arr(:)

    n = 1
    do i = 1, size(vars)
      n = n*size(global%ext_var_data(vars(i))%data)
    enddo

    allocate (global%iterables(n, size(global%continuous_indices)), &
              temp_arr(size(global%ext_var_data)))

    global%iterables = 1
    temp_arr = 0

    do i = 1, n

      reduction = i
      do counter = size(vars), 1, -1
        if (counter == 1) then
          temp_arr(vars(counter)) = reduction
        else
          temp_arr(vars(counter)) = modulo(reduction, size(global%ext_var_data(vars(counter))%data))
          if (temp_arr(vars(counter)) == 0) temp_arr(vars(counter)) = size(global%ext_var_data(vars(counter))%data)
          reduction = int((reduction - temp_arr(vars(counter)))/size(global%ext_var_data(vars(counter))%data)) + 1
        endif
      enddo

      do j = 1, size(vars)
        global%iterables(i, vars(j)) = temp_arr(vars(j))
      enddo

    enddo

    deallocate (temp_arr)

  end subroutine SsTC_construct_iterable

end module SsTC_data_structures
