module data_structures

  use utility

  implicit none

  type sys
    character(len=120)            :: name
    integer                       :: num_bands                                !Number of bands.
    real(kind=dp)                 :: direct_lattice_basis(3, 3)               !Direct lattice basis vectors in A. 1st index is vector label, 2nd index is vector component.
    real(kind=dp)                 :: metric_tensor(3, 3)                      !Metric tensor of the direct lattice basis.
    real(kind=dp)                 :: cell_volume                              !Cell volume.
    integer                       :: num_R_points                             !Number of R points given in the *_tb.dat file.
    integer, allocatable          :: R_point(:, :)                            !Memory layout id of the R-point (1st index) and R-vector coords. relative to the direct lattice basis vectors (2nd index).
    integer, allocatable          :: deg_R_point(:)                           !Degeneracy of the R-point specified by its memory layout id.
    complex(kind=dp), allocatable :: real_space_hamiltonian_elements(:, :, :) !Hamiltonian matrix elements (1st and 2nd indexes) and memory layout id of the R-point (3rd index) in eV.
    complex(kind=dp), allocatable :: real_space_position_elements(:, :, :, :) !Position operator matrix elements (1st and 2nd indexes), cartesian coordinate (3rd index) and memory layout id of the R-point (4th index) in A.
    real(kind=dp)                 :: e_fermi = 0.0_dp                         !Fermi energy.
    real(kind=dp)                 :: deg_thr = 1.0E-4_dp                      !Degeneracy threshold in eV.
    real(kind=dp)                 :: deg_offset = 0.04_dp                     !Offset for regularization in case of deeneracies in eV.
    !Below: Floquet stuff.
    real(kind=dp)                 :: Nt = 65 !2^6 + 1 Discretization points of each period.
    real(kind=dp)                 :: Ns = 10 !Considered Harmonics.
    logical                       :: diag = .false. !If we only consider diagonal terms of the pos. operator.
    !Optical stuff.
    logical                       :: adpt_smearing = .true.
    real(kind=dp)                 :: smearing = 1.0_dp
  end type sys

  type local_k_data
    character(len=120)                            :: name
    integer, allocatable                          :: integer_indices(:)               !Each entry contains the range of each of the integer indices.
    complex(kind=dp), allocatable                 :: k_data(:)                        !Data local for each k with integer index in memory array,
    procedure (local_calculator), pointer, nopass :: local_calculator                 !Pointer to the local calculator.
    integer                                       :: particular_integer_component = 0 !Specification of some integer component.
  end type local_k_data

  type, extends(local_k_data) :: global_k_data
    integer,          allocatable                  :: continuous_indices(:) !Each entry contains the range of each continuous indices.
    type(external_vars), allocatable               :: ext_var_data(:)       !External variable data.
    procedure (global_calculator), pointer, nopass :: global_calculator     !Pointer to the global calculator.
    integer, allocatable                           :: iterables(:, :)
  end type global_k_data

  type external_vars
    real(kind=dp), allocatable :: data(:) !External variable data array.
  end type external_vars

  !Interfaces for the generic functions returning k dependent quantities.
  abstract interface
    function global_calculator(task, system, k, error) result(u)
      import :: global_k_data, sys, external_vars, dp
      class(global_k_data), intent(in) :: task
      type(sys),              intent(in) :: system
      real(kind=dp),          intent(in) :: k(3)
      logical, intent(inout) :: error

      complex(kind=dp)                   :: u(product(task%integer_indices), product(task%continuous_indices))
    end function global_calculator
  end interface

  abstract interface
    function local_calculator(k_data, system, k, error) result(u)
      import :: local_k_data, sys, external_vars, dp
      class(local_k_data), intent(in) :: k_data
      type(sys),          intent(in)  :: system
      real(kind=dp),      intent(in)  :: k(3)
      logical, intent(inout) :: error

      complex(kind=dp)                :: u(product(k_data%integer_indices))
    end function local_calculator
  end interface

  public :: local_k_data
  public :: global_k_data

  public :: sys
  public :: sys_constructor

  public :: external_vars
  public :: external_variable_constructor

  public :: integer_array_element_to_memory_element
  public :: integer_memory_element_to_array_element
  public :: continuous_array_element_to_memory_element
  public :: continuous_memory_element_to_array_element
  public :: construct_iterable

  contains

  function sys_constructor(name, path_to_tb_file, efermi, deg_thr, deg_offset, &
                           floq_Nt, floq_Ns, floq_diag, &
                           optical_smearing) result(system)

    character(len=*),        intent(in) :: name
    character(len=*),        intent(in) :: path_to_tb_file
    real(kind=dp), optional, intent(in) :: efermi, &
                                           deg_thr, deg_offset, &
                                           optical_smearing
    integer, optional, intent(in) :: floq_Nt, floq_NS
    logical, optional, intent(in) :: floq_diag

    type(sys) :: system

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
    !==FLOQUET==!
    if (present(floq_Nt)) system%Nt = floq_Nt
    if (present(floq_Ns)) system%Ns = floq_Ns
    if (present(floq_diag)) system%diag = floq_diag
    !==OPTICAL==!
    if (present(optical_smearing)) then
      system%adpt_smearing = .false.
      system%smearing = optical_smearing
    endif

    filename = trim(path_to_tb_file)//trim(name)//"_tb.dat"
    filename = trim(filename)

    write(unit=112, fmt = "(A)") "Reading file"//filename//"."

    open(unit=113, action="read", file = filename)
    read(unit=113, fmt = *)
    do i = 1, 3
      read(unit=113, fmt = *) (system%direct_lattice_basis(i, j), j = 1, 3)
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

    read(unit=113, fmt = *) num_bands
    system%num_bands = num_bands
    read(unit=113, fmt = *) nrpts
    system%num_R_points = nrpts
    allocate(system%R_point(system%num_R_points, 3), &
             system%deg_R_point(system%num_R_points), &
             system%real_space_hamiltonian_elements(system%num_bands, system%num_bands, system%num_R_points), &
             system%real_space_position_elements(system%num_bands, system%num_bands, 3, system%num_R_points))

    division = nrpts/15
    remainder = modulo(nrpts, 15)

    if (remainder == 0) then
      do i = 1, division
          read(unit=113, fmt = *) (system%deg_R_point(15*(i - 1) + j), j = 1, 15)
      enddo
    else
      do i = 1, division
        read(unit=113, fmt = *) (system%deg_R_point(15*(i - 1) + j), j = 1, 15)
      enddo
      read(unit=113, fmt = *) (system%deg_R_point(15*(i - 1) + j), j = 1, remainder)
    endif

    read(unit=113, fmt = *)
    allocate(dummyR(2))
    write(unit=112, fmt = "(A)") "Reading Hamiltonian."

    do irpts = 1, nrpts
      read(unit=113, fmt = *) (system%R_point(irpts, i), i = 1, 3)
      do i = 1, num_bands
        do j = 1, num_bands
          read(unit=113, fmt = *) dummy1, dummy2, dummyR(1), dummyR(2)
          system%real_space_hamiltonian_elements(i, j, irpts) = cmplx(dummyR(1), dummyR(2), dp)
        enddo
      enddo
      read(unit=113, fmt = *)
    enddo

    deallocate(dummyR)
    write(unit=112, fmt = "(A)") "Done."

    allocate(dummyR(6))
    write(unit=112, fmt = "(A)") "Reading Position Operator."

    do irpts = 1, nrpts - 1
      read(unit=113, fmt = *) dummy1, dummy2, dummy3
      do i = 1, num_bands
        do j = 1, num_bands
          read(unit=113, fmt = *) dummy1, dummy2, dummyR(1), dummyR(2), dummyR(3), dummyR(4), dummyR(5), dummyR(6)
          system%real_space_position_elements(i, j, 1, irpts) = cmplx(dummyR(1), dummyR(2), dp)
          system%real_space_position_elements(i, j, 2, irpts) = cmplx(dummyR(3), dummyR(4), dp)
          system%real_space_position_elements(i, j, 3, irpts) = cmplx(dummyR(5), dummyR(6), dp)
        enddo
      enddo
      read(unit=113, fmt = *)
    enddo

    read(unit=113, fmt = *) dummy1, dummy2, dummy3
    do i = 1, num_bands
      do j = 1, num_bands
        read(unit=113, fmt = *) dummy1, dummy2, dummyR(1), dummyR(2), dummyR(3), dummyR(4), dummyR(5), dummyR(6)
        system%real_space_position_elements(i, j, 1, nrpts) = cmplx(dummyR(1), dummyR(2), dp)
        system%real_space_position_elements(i, j, 2, nrpts) = cmplx(dummyR(3), dummyR(4), dp)
        system%real_space_position_elements(i, j, 3, nrpts) = cmplx(dummyR(5), dummyR(6), dp)
      enddo
    enddo

    deallocate(dummyR)
    write(unit=112, fmt = "(A)") "Done."

    close(unit=113)

  end function sys_constructor

  function external_variable_constructor(start, end, steps) result(vars)
    !Function to set external variable data.
    real(kind=dp), intent(in) :: start, end
    integer,       intent(in) :: steps

    type (external_vars) :: vars

    integer :: i

    allocate(vars%data(steps))

    if (steps == 1) then
      vars%data(1) = start
    else
      do i = 1, steps
        vars%data(i) = start + (end - start)*real(i - 1,dp)/real(steps - 1,dp)
      enddo
    endif

  end function external_variable_constructor

  !==UTILITY FUNCTIONS TO PASS "MEMORY LAYOUT" INDICES TO "ARRAY LAYOUT" AND VICE-VERSA==!
  !See discussion at https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays

  function integer_array_element_to_memory_element(data_k, i_arr) result (i_mem)
    !Get integer indices from array layout to memory layout.
    class(local_k_data), intent(in) :: data_k
    integer,                intent(in) :: i_arr(size(data_k%integer_indices))

    integer :: i_mem

    integer :: counter

    i_mem = 1
    do counter = 0, size(data_k%integer_indices) - 1
      i_mem = (i_mem - 1)*data_k%integer_indices(counter + 1) + i_arr(counter + 1)
    enddo

  end function integer_array_element_to_memory_element

  function integer_memory_element_to_array_element(data_k, i_mem) result (i_arr)
    !Get integer indices from memory layout to array layout.
    class(local_k_data), intent(in) :: data_k
    integer, intent(in)                :: i_mem

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

  end function integer_memory_element_to_array_element

  function continuous_array_element_to_memory_element(task, r_arr) result (r_mem)
    !Get continuous indices from array layout to memory layout.
    class(global_k_data), intent(in) :: task
    integer,                intent(in) :: r_arr(size(task%continuous_indices))

    integer :: r_mem

    integer :: counter

    r_mem = 1
    do counter = 0, size(task%continuous_indices) - 1
      r_mem = (r_mem - 1)*task%continuous_indices(counter + 1) + r_arr(counter + 1)
    enddo

  end function continuous_array_element_to_memory_element

  function continuous_memory_element_to_array_element(task, r_mem) result (r_arr)
    !Get continuous indices from memory layout to array layout.
    class(global_k_data), intent(in) :: task
    integer,                intent(in) :: r_mem

    integer :: r_arr(size(task%continuous_indices))

    integer :: counter, reduction

    reduction = r_mem
    do counter = size(task%continuous_indices), 1, -1
      if (counter == 1) then
        r_arr(counter) = reduction
      else
        r_arr(counter) = modulo(reduction, task%continuous_indices(counter))
        if (r_arr(counter)==0) r_arr(counter) = task%continuous_indices(counter)
        reduction = int((reduction - r_arr(counter))/task%continuous_indices(counter)) + 1
      endif
    enddo

  end function continuous_memory_element_to_array_element

  subroutine construct_iterable(global, vars)
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
    class(global_k_data) :: global
    integer, intent(in) :: vars(:)

    integer :: i, j, n, counter, reduction

    integer, allocatable :: temp_arr(:)

    n = 1
    do i = 1, size(vars)
      n = n * size(global%ext_var_data(vars(i))%data)
    enddo

    allocate(global%iterables(n, size(global%continuous_indices)), &
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

    deallocate(temp_arr)
    

  end subroutine construct_iterable

end module data_structures