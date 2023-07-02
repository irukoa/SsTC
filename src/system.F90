module system

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

  type BZ_integrated_data
    !Label for the integration task.
    character(len=120)            :: task
    !Specification of the integer indices of the quantity which will be integrated.
    integer,          allocatable :: integer_indices(:)
    !Specification of the continuous indices of the quantity which will be integrated.
    integer,          allocatable :: continuous_indices(:)
    !Result of the integration, contains the integer index and the continuous index in memory layout, respectively.
    complex(kind=dp), allocatable :: res(:, :)
    !If only some component out of the integer indices wants to be calculated this can be specified here in memory layout.
    integer                       :: particular_integer_component = 0 
    !Pointer to calculator function.
    procedure (calculator_C1M3), pointer, nopass :: calculator
    !Integration method.
    character(len=120)            :: method
    !Integration samples.
    integer                       :: samples(3)
    
    type(external_vars), allocatable :: ext_var_data(:)
  end type BZ_integrated_data

  type external_vars
    !External variable data array.
    real(kind=dp), allocatable    :: data(:)
  end type external_vars


  abstract interface
    function calculator_C1M3(task, system, k) result(u)
      import :: BZ_integrated_data, sys, external_vars, dp
      type(BZ_integrated_data), intent(in) :: task
      type(sys),                intent(in) :: system
      real(kind=dp),            intent(in) :: k(3)

      complex(kind=dp)                     :: u(product(task%integer_indices), product(task%continuous_indices))
    end function calculator_C1M3
  end interface

  public :: sys
  public :: BZ_integrated_data
  public :: external_vars

  public :: print_task_result

  public :: integer_array_element_to_memory_element
  public :: integer_memory_element_to_array_element
  public :: continuous_array_element_to_memory_element
  public :: continuous_memory_element_to_array_element

  contains

  !I will do the constructor for the system later on the line.

  function external_variable_constructor(start, end, steps) result(vars)
    real(kind=dp), intent(in) :: start, end
    integer, intent(in) :: steps

    type (external_vars) :: vars

    integer :: i

    allocate(vars%data(steps))

    if (steps == 1) then
      vars%data(1) = start
    else
      do i = 1, steps
        vars%data(i) = start + (end - start)*real(i-1,dp)/real(steps-1,dp)
      enddo
    endif

  end function external_variable_constructor

  function task_constructor(name, Nint, int_range, Next_vars, ext_vars_start, ext_vars_end, ext_vars_steps, calculator, part_int_comp, method, samples) result(task)
    character(len=*), intent(in)  :: name
    integer, intent(in)           :: Nint
    integer, intent(in)           :: int_range(Nint)
    integer, intent(in) :: Next_vars
    real(kind=dp), optional, intent(in) :: ext_vars_start(Next_vars), ext_vars_end(Next_vars)
    integer, optional, intent(in) :: ext_vars_steps(Next_vars)
    procedure (calculator_C1M3)   :: calculator
    integer, optional, intent(in) :: part_int_comp(Nint)
    character(len=*), optional, intent(in)  :: method
    integer, intent(in) :: samples(3)

    type(BZ_integrated_data) :: task

    integer :: i

    task%task = name

    allocate(task%integer_indices(Nint))
    do i = 1, Nint
      task%integer_indices(i) = int_range(i)
    enddo

    if (((Next_vars).ge.1).and.(present(ext_vars_start)).and.(present(ext_vars_end)).and.(present(ext_vars_steps))) then
      allocate(task%continuous_indices(Next_vars), task%ext_var_data(Next_vars))
      do i = 1, Next_vars
        task%continuous_indices(i) = ext_vars_steps(i)
        allocate(task%ext_var_data(i)%data(ext_vars_steps(i)))
        task%ext_var_data(i) = external_variable_constructor(ext_vars_start(i), ext_vars_end(i), ext_vars_steps(i))
      enddo
    else
      allocate(task%continuous_indices(1))
      task%continuous_indices(1) = 1.0_dp
    endif

    task%calculator => calculator

    allocate(task%res(product(task%integer_indices), product(task%continuous_indices)))

    if (present(part_int_comp)) task%particular_integer_component = integer_array_element_to_memory_element(task, part_int_comp)

    if (present(method)) then
      if (method == "extrapolation") then
        task%method = "extrapolation"
      elseif (method == "rectangle") then
        task%method = "rectangle"
      else
        print*, "Integration method not recognized. Setting up rectangle method"
        task%method = "rectangle"
      endif
    else
      task%method = "rectangle"
    endif

    task%samples = samples

  end function task_constructor

  subroutine print_task_result(task, system)
    !Subroutine to format and output files.
    type(BZ_integrated_data), intent(in) :: task
    type(sys),                intent(in) :: system

    character(len=400) :: filename
    integer :: i_arr(size(task%integer_indices)), r_arr(size(task%continuous_indices))
    integer :: i_mem, r_mem, count

    do i_mem = 1, product(task%integer_indices) !For each integer index.

      i_arr = integer_memory_element_to_array_element(task, i_mem) !Pass to array layout.

      filename = trim(system%name)//'-'//trim(task%task)//'_'
      do count = 1, size(task%integer_indices)
        filename = trim(filename)//achar(48 + i_arr(count))
      enddo
      filename = trim(filename)//'.dat'

      open(unit=111, action="write", file=filename)

      do r_mem = 1, product(task%continuous_indices) !For each continuous index.

        r_arr = continuous_memory_element_to_array_element(task, r_mem) !Pass to array layout.

        write(unit=111, fmt=*) (task%ext_var_data(count)%data(r_arr(count)), count = 1, size(task%continuous_indices)),&
        real(task%res(i_mem, r_mem), dp), aimag(task%res(i_mem, r_mem))

      enddo

      close(unit=111)

    enddo

  end subroutine print_task_result

  function integer_array_element_to_memory_element(task, i_arr) result (i_mem)
    !Get integer indices from array layout to memory layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: i_arr(size(task%integer_indices))
    integer :: i_mem

    integer :: counter

    i_mem = 1
    do counter = 0, size(task%integer_indices) - 1
      i_mem = (i_mem-1)*task%integer_indices(counter+1) + i_arr(counter+1)
    enddo

  end function integer_array_element_to_memory_element

  function integer_memory_element_to_array_element(task, i_mem) result (i_arr)
    !Get integer indices from memory layout to array layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: i_mem
    integer :: i_arr(size(task%integer_indices))

    integer :: counter, reduction

    reduction = i_mem
    do counter = size(task%integer_indices), 1, -1
      if (counter == 1) then
        i_arr(counter) = reduction
      else
        i_arr(counter) = modulo(reduction, task%integer_indices(counter))
        if (i_arr(counter)==0) i_arr(counter) = task%integer_indices(counter)
        reduction = int((reduction-i_arr(counter))/task%integer_indices(counter)) + 1
      endif
    enddo

  end function integer_memory_element_to_array_element

  function continuous_array_element_to_memory_element(task, r_arr) result (r_mem)
    !Get continuous indices from array layout to memory layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: r_arr(size(task%continuous_indices))
    integer :: r_mem

    integer :: counter

    r_mem = 1
    do counter = 0, size(task%continuous_indices) - 1
      r_mem = (r_mem-1)*task%continuous_indices(counter+1) + r_arr(counter+1)
    enddo

  end function continuous_array_element_to_memory_element

  function continuous_memory_element_to_array_element(task, r_mem) result (r_arr)
    !Get continuous indices from memory layout to array layout.
    type(BZ_integrated_data), intent(in) :: task
    integer, intent(in) :: r_mem
    integer :: r_arr(size(task%continuous_indices))

    integer :: counter, reduction

    reduction = r_mem
    do counter = size(task%continuous_indices), 1, -1
      if (counter == 1) then
        r_arr(counter) = reduction
      else
        r_arr(counter) = modulo(reduction, task%continuous_indices(counter))
        if (r_arr(counter)==0) r_arr(counter) = task%continuous_indices(counter)
        reduction = int((reduction-r_arr(counter))/task%continuous_indices(counter)) + 1
      endif
    enddo

  end function continuous_memory_element_to_array_element

end module system