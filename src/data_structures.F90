module data_structures

  use utility

  implicit none

  type sys
    character(len=120)            :: name
    integer                       :: num_bands !Number of bands.
    real(kind=dp)                 :: direct_lattice_basis(3, 3) !Direct lattice basis vectors in A. 1st index is vector label, 2nd index is vector component.
    integer                       :: num_R_points !Number of R points given in the *_tb.dat file.
    integer, allocatable          :: R_point(:, :) !Memory layout id of the R-point (1st index) and R-vector coords. relative to the direct lattice basis vectors (2nd index).
    integer, allocatable          :: deg_R_point(:) !Degeneracy of the R-point specified by its memory layout id.
    complex(kind=dp), allocatable :: real_space_hamiltonian_elements(:, :, :) !Hamiltonian matrix elements (1st and 2nd indexes) and memory layout id of the R-point (3rd index) in eV.
    complex(kind=dp), allocatable :: real_space_position_elements(:, :, :, :) !Position operator matrix elements (1st and 2nd indexes), cartesian coordinate (3rd index) and memory layout id of the R-point (4th index) in A.
    real(kind=dp)                 :: e_fermi = 0.0_dp !Fermi energy.
    real(kind=dp)                 :: deg_thr = 1.0E-4_dp !Degeneracy threshold in eV.
    real(kind=dp)                 :: deg_offset = 0.04_dp !Offset for regularization in case of deeneracies in eV.
  end type sys

  type local_k_data
    !Name of the local quantity or in the extended type case, the name of the task to integrate.
    character(len=120)                            :: name
    !Specification of the integer indices of the local quantity, or in the extended type case,
    !the integer indices which will be integrated.
    integer, allocatable                          :: integer_indices(:)
    !Local k-data.
    complex(kind=dp), allocatable                 :: k_data(:)
    !Pointer to calculator function.
    procedure (local_calculator), pointer, nopass :: local_calculator
    !If only some component out of the integer indices wants to be calculated this can be specified here in memory layout.
    integer                                        :: particular_integer_component = 0 
  end type local_k_data

  type, extends(local_k_data) :: global_k_data
    !Specification of the continuous indices of the quantity which will be integrated.
    integer,          allocatable                  :: continuous_indices(:)
    !External variable data.
    type(external_vars), allocatable               :: ext_var_data(:)
    !Result of the integration, contains the integer index and the continuous index in memory layout, respectively.
    complex(kind=dp), allocatable                  :: result(:, :)
    !Pointer to calculator function.
    procedure (global_calculator), pointer, nopass :: global_calculator
  end type global_k_data

  type, extends(global_k_data) :: BZ_integral_task
    !Integration method.
    character(len=120)                             :: method
    !Integration samples.
    integer                                        :: samples(3)
  end type BZ_integral_task

  type external_vars
    !External variable data array.
    real(kind=dp), allocatable :: data(:)
  end type external_vars

  !Interface for the generic function returning the integrand of the quantity to be integrated at point "k",
  !a.k.a. the "calculator".
  abstract interface
    function global_calculator(task, system, k) result(u)
      import :: global_k_data, sys, external_vars, dp
      class(global_k_data), intent(in) :: task
      type(sys),              intent(in) :: system
      real(kind=dp),          intent(in) :: k(3)

      complex(kind=dp)                   :: u(product(task%integer_indices), product(task%continuous_indices))
    end function global_calculator
  end interface

  abstract interface
    function local_calculator(k_data, system, k) result(u)
      import :: local_k_data, sys, external_vars, dp
      class(local_k_data), intent(in) :: k_data
      type(sys),          intent(in)  :: system
      real(kind=dp),      intent(in)  :: k(3)

      complex(kind=dp)                :: u(product(k_data%integer_indices))
    end function local_calculator
  end interface

  public  :: sys
  public  :: sys_constructor

  public  :: local_k_data

  public  :: global_k_data
  public  :: task_constructor

  public :: external_vars
  public :: external_variable_constructor

  public  :: print_task_result

  public  :: integer_array_element_to_memory_element
  public  :: integer_memory_element_to_array_element
  public  :: continuous_array_element_to_memory_element
  public  :: continuous_memory_element_to_array_element

  contains

  function sys_constructor(name, path_to_tb_file, efermi) result(system)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: path_to_tb_file
    real(kind=dp), optional, intent(in) :: efermi

    type(sys) :: system

    character(len=400) :: filename
    integer            :: num_bands, nrpts, &
                          division, remainder, &
                          irpts, i, j, &
                          dummy1, dummy2, dummy3

    real(kind=dp), allocatable :: dummyR(:)

    system%name = name

    if(present(efermi)) system%e_fermi = efermi

    filename = trim(path_to_tb_file)//trim(name)//"_tb.dat"
    filename = trim(filename)

    write(unit=112, fmt = "(A)") "Reading file"//filename//"."

    open(unit=113, action="read", file = filename)
    read(unit=113, fmt = *)
    do i = 1, 3
      read(unit=113, fmt = *) (system%direct_lattice_basis(i, j), j = 1, 3)
    enddo
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

  function task_constructor(name, calculator, &
                            N_int_ind, int_ind_range, &
                            N_ext_vars, ext_vars_start, ext_vars_end, ext_vars_steps, &
                            method, samples, &
                            part_int_comp) result(task)
    !Function to construct the data structure "task" containing all info for BN integration. 
    !This function sets integer indices, external variables, integration method 
    !and number of points in the BZ discretization. If only a particular integer component 
    !wants to be computed, this can also be set. The BZ integrand must be given by function
    !"calculator" with the interface calculator_C1M3.

    character(len=*),            intent(in) :: name
    procedure (global_calculator)           :: calculator
    integer, intent(in)                     :: N_int_ind
    integer, intent(in)                     :: int_ind_range(N_int_ind)
    integer, intent(in)                     :: N_ext_vars
    real(kind=dp),    optional,  intent(in) :: ext_vars_start(N_ext_vars), ext_vars_end(N_ext_vars)
    integer,          optional,  intent(in) :: ext_vars_steps(N_ext_vars)
    character(len=*), optional,  intent(in) :: method
    integer, intent(in)                     :: samples(3)
    integer,          optional,  intent(in) :: part_int_comp(N_int_ind)

    type(BZ_integral_task) :: task

    integer :: i

    !Set name.
    task%name = name

    !Set integer index data.
    allocate(task%integer_indices(N_int_ind))
    do i = 1, N_int_ind
      task%integer_indices(i) = int_ind_range(i)
    enddo

    !Set external variable data.
    if (((N_ext_vars).ge.1).and.(present(ext_vars_start)).and.(present(ext_vars_end)).and.(present(ext_vars_steps))) then
      allocate(task%continuous_indices(N_ext_vars), task%ext_var_data(N_ext_vars))
      do i = 1, N_ext_vars
        task%continuous_indices(i) = ext_vars_steps(i)
        allocate(task%ext_var_data(i)%data(ext_vars_steps(i)))
        task%ext_var_data(i) = external_variable_constructor(ext_vars_start(i), ext_vars_end(i), ext_vars_steps(i))
      enddo
    else
      allocate(task%continuous_indices(1))
      task%continuous_indices(1) = 1.0_dp
    endif

    !Set calculator pointer (function alias).
    task%global_calculator => calculator

    !Set result.
    allocate(task%result(product(task%integer_indices), product(task%continuous_indices)))

    !Set calculation of a particular integer component.
    if (present(part_int_comp)) task%particular_integer_component = integer_array_element_to_memory_element(task, part_int_comp)

    !Set integration method.
    if (present(method)) then
      if (method == "extrapolation") then
        task%method = "extrapolation"
        write(unit=112, fmt="(A)") "Warning: To employ the extrapolation method all the elements of the array 'samples' must be expressible as either 1 or 2^n + 1 for n = 0, 1, ..."
      elseif (method == "rectangle") then
        task%method = "rectangle"
      else
        print*, "Integration method not recognized. Setting up rectangle method"
        task%method = "rectangle"
      endif
    else
      task%method = "rectangle"
    endif

    !Set number of integration samples (discretization of BZ).
    task%samples = samples

  end function task_constructor

  subroutine print_task_result(task, system)
    !Subroutine to format and output files related to the result of the task "task".
    class(global_k_data), intent(in) :: task
    type(sys),              intent(in) :: system

    character(len=400) :: filename
    integer :: i_arr(size(task%integer_indices)), &
               r_arr(size(task%continuous_indices))
    integer :: i_mem, r_mem, count

    do i_mem = 1, product(task%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

      filename = trim(system%name)//'-'//trim(task%name)//'_'
      do count = 1, size(task%integer_indices)
        filename = trim(filename)//achar(48 + i_arr(count))
      enddo
      filename = trim(filename)//'.dat'

      open(unit=111, action="write", file=filename)

      do r_mem = 1, product(task%continuous_indices) !For each continuous index.

        r_arr = continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.

        write(unit=111, fmt=*) (task%ext_var_data(count)%data(r_arr(count)), count = 1, size(task%continuous_indices)),&
        real(task%result(i_mem, r_mem), dp), aimag(task%result(i_mem, r_mem))

      enddo

      close(unit=111)

    enddo

  end subroutine print_task_result

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

end module data_structures