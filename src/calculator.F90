module calculator

  use utility
  use system

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
  end type BZ_integrated_data

  type external_vars
    real(kind=dp), allocatable :: data(:)
  end type external_vars

  !=======EXAMPLE=======!
  !Let's say that we want to integrate the complex optical conductivity \sigma^{ab}(\omega).
  !In this case, a, b are integer indices ranging each from 1 to 3 and \omega is a continuous
  !variable which we suppose has been discretized in N points in some range as 
  !\omega_i = \omega_1 + (\omega_N - \omega_1)*(i-1)/(N-1).
  !To be able to compute this quantity, the following things are needed:
  !i) A specification of the derived type associated to the task of calculating \sigma^{ab}(\omega).
  !For example, 
  !
  !   type(BZ_integrated_data) :: complex_optical_conductivity
  !
  !   complex_optical_conductivity%task = "opt_cond"                      <- Alias for the task.
  !   allocate(complex_optical_conductivity%integer_indices(2))           <- Meaning that we consider indices a and b.
  !   complex_optical_conductivity%integer_indices(1) = 3                 <- Meaning that a takes ONLY the values 1, 2, 3.
  !   complex_optical_conductivity%integer_indices(2) = 3                 <- Meaning that b takes ONLY the values 1, 2, 3.
  !   allocate(complex_optical_conductivity%continuous_indices(1))        <- Meaning that we consider the continuous variable \omega.
  !   complex_optical_conductivity%continuous_indices(1) = N              <- Meaning that \omega_i can ONLY take the N values 1, 2, ..., N.
  !   allocate(complex_optical_conductivity%res( &                        <- Allocate result.
  !            product(complex_optical_conductivity%integer_indices),&
  !            product(complex_optical_conductivity%continuous_indices)))
  !
  !ii) An if clause on the calculator_dict subroutine as follows:
  !
  !   subroutine calculator_dict(task, i_arr, r_arr, k, data_k)
  !   .
  !   .
  !     if (task%task == "opt_cond") then
  !       data_k = calculator_opt_cond(k, i_arr, r_arr)
  !     elseif
  !   .
  !   .
  !   end subroutine calculator_dict
  !iii) A function, calculator_opt_cond defined as public in some module in the scope of program 
  !floquet_tight_binding (main.F90) and this module, which takes as inputs
  !a, b, \omega_i and k and returns a scalar which contains the contribution to the integral
  !of a given kpt.

  public :: BZ_integrated_data
  public :: external_vars
  public :: print_task_result
  public :: calculator_dict
  public :: integer_array_element_to_memory_element
  public :: integer_memory_element_to_array_element
  public :: continuous_array_element_to_memory_element
  public :: continuous_memory_element_to_array_element

  contains

  subroutine print_task_result(task, system, external_variables)
    !Subroutine to format and output files.
    type(BZ_integrated_data), intent(in) :: task
    type(sys),                intent(in) :: system
    type(external_vars),      intent(in) :: external_variables(size(task%continuous_indices))

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

        write(unit=111, fmt=*) (external_variables(count)%data(r_arr(count)), count = 1, size(task%continuous_indices)),&
        real(task%res(i_mem, r_mem), dp), aimag(task%res(i_mem, r_mem))

      enddo

      close(unit=111)

    enddo

  end subroutine print_task_result

  subroutine calculator_dict(task, i_arr, r_arr, k, data_k)
    !General calculator dictionary. Specify the task and refer to a 
    !function that returns the contribution of the given kpt and indices.

    type(BZ_integrated_data), intent(in)    :: task
    real(kind=dp),            intent(in)    :: k(3)
    integer,                  intent(in)    :: i_arr(:), r_arr(:)
    complex(kind=dp),         intent(inout) :: data_k

    if (task%task == "test") then
      !data_k = calculator_test(k, i_arr, r_arr)
    else
      print*, "Task label not recognized."
      stop
    endif
  end subroutine calculator_dict

  !function calculator_test(k, i_arr, r_arr) result(u)
    !TESTING PURPOSES CALCULATOR.
  !  real(kind=dp), intent(in) :: k(3)
  !  integer,       intent(in) :: i_arr(:), r_arr(:)

  !  complex(kind=dp) :: u

  !  u = r_arr(1)*i_arr(1)*(k(1)**2)*exp(sin(10*k(1)))

  !end function calculator_test

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

end module calculator