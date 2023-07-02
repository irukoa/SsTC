module data_structures

  implicit none

  integer, parameter, private :: dp = 8

  type sys
    character(len=120)            :: name
    integer                       :: num_bands
    real(kind=dp)                 :: lattice_constant !In A.
    real(kind=dp)                 :: direct_lattice_basis(3, 3), &  !Relative to lattice_constant.
                                     reciprocal_lattice_basis(3, 3) !Relative to 1/lattice_constant.
    complex(kind=dp), allocatable :: real_space_hamiltonian_elements(:, :, :) !In eV.
    complex(kind=dp), allocatable :: real_space_position_elements(:, :, :, :) !In A.
  end type sys

  type BZ_integral_task
    !Label for the integration task.
    character(len=120)                           :: name
    !Specification of the integer indices of the quantity which will be integrated.
    integer,          allocatable                :: integer_indices(:)
    !Specification of the continuous indices of the quantity which will be integrated.
    integer,          allocatable                :: continuous_indices(:)
    !External variable data.
    type(external_vars), allocatable             :: ext_var_data(:)
    !Result of the integration, contains the integer index and the continuous index in memory layout, respectively.
    complex(kind=dp), allocatable                :: result(:, :)
    !Pointer to calculator function.
    procedure (calculator_C1M3), pointer, nopass :: calculator
    !Integration method.
    character(len=120)                           :: method
    !Integration samples.
    integer                                      :: samples(3)
    !If only some component out of the integer indices wants to be calculated this can be specified here in memory layout.
    integer                                      :: particular_integer_component = 0 
  end type BZ_integral_task

  type external_vars
    !External variable data array.
    real(kind=dp), allocatable :: data(:)
  end type external_vars

  !Interface for the generic function returning the integrand of the quantity to be integrated at point "k",
  !a.k.a. the "calculator".
  abstract interface
    function calculator_C1M3(task, system, k) result(u)
      import :: BZ_integral_task, sys, external_vars, dp
      type(BZ_integral_task), intent(in) :: task
      type(sys),                intent(in) :: system
      real(kind=dp),            intent(in) :: k(3)

      complex(kind=dp)                     :: u(product(task%integer_indices), product(task%continuous_indices))
    end function calculator_C1M3
  end interface

  public  :: sys

  public  :: BZ_integral_task
  public  :: task_constructor

  private :: external_vars
  private :: external_variable_constructor

  public  :: print_task_result

  public  :: integer_array_element_to_memory_element
  public  :: integer_memory_element_to_array_element
  public  :: continuous_array_element_to_memory_element
  public  :: continuous_memory_element_to_array_element

  contains

  !I will do the constructor for the system later on the line.

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
    procedure (calculator_C1M3)             :: calculator
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
    task%calculator => calculator

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
    type(BZ_integral_task), intent(in) :: task
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

  function integer_array_element_to_memory_element(task, i_arr) result (i_mem)
    !Get integer indices from array layout to memory layout.
    type(BZ_integral_task), intent(in) :: task
    integer,                intent(in) :: i_arr(size(task%integer_indices))

    integer :: i_mem

    integer :: counter

    i_mem = 1
    do counter = 0, size(task%integer_indices) - 1
      i_mem = (i_mem - 1)*task%integer_indices(counter + 1) + i_arr(counter+1)
    enddo

  end function integer_array_element_to_memory_element

  function integer_memory_element_to_array_element(task, i_mem) result (i_arr)
    !Get integer indices from memory layout to array layout.
    type(BZ_integral_task), intent(in) :: task
    integer, intent(in)                :: i_mem

    integer :: i_arr(size(task%integer_indices))

    integer :: counter, reduction

    reduction = i_mem
    do counter = size(task%integer_indices), 1, -1
      if (counter == 1) then
        i_arr(counter) = reduction
      else
        i_arr(counter) = modulo(reduction, task%integer_indices(counter))
        if (i_arr(counter) == 0) i_arr(counter) = task%integer_indices(counter)
        reduction = int((reduction - i_arr(counter))/task%integer_indices(counter)) + 1
      endif
    enddo

  end function integer_memory_element_to_array_element

  function continuous_array_element_to_memory_element(task, r_arr) result (r_mem)
    !Get continuous indices from array layout to memory layout.
    type(BZ_integral_task), intent(in) :: task
    integer,                intent(in) :: r_arr(size(task%continuous_indices))

    integer :: r_mem

    integer :: counter

    r_mem = 1
    do counter = 0, size(task%continuous_indices) - 1
      r_mem = (r_mem - 1)*task%continuous_indices(counter+1) + r_arr(counter + 1)
    enddo

  end function continuous_array_element_to_memory_element

  function continuous_memory_element_to_array_element(task, r_mem) result (r_arr)
    !Get continuous indices from memory layout to array layout.
    type(BZ_integral_task), intent(in) :: task
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