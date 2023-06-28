module calculator

  use utility

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
  !   complex_optical_conductivity%continuous_indices(1) = 1              <- Meaning that \omega_i can ONLY take the N values 1, 2, ..., N.
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
  public :: calculator_dict

  contains

  subroutine calculator_dict(task, i_arr, r_arr, k, data_k)
    !General calculator dictionary. Specify the task and refer to a 
    !function that returns the contribution of the given kpt and indices.

    type(BZ_integrated_data), intent(in)    :: task
    real(kind=dp),            intent(in)    :: k(3)
    integer,                  intent(in)    :: i_arr(:), r_arr(:)
    complex(kind=dp),         intent(inout) :: data_k

    if (task%task == "test") then
      data_k = calculator_test(k, i_arr, r_arr)
    else
      print*, "Task label not recognized."
      stop
    endif
  end subroutine calculator_dict

  function calculator_test(k, i_arr, r_arr) result(u)
    !TESTING PURPOSES CALCULATOR.
    real(kind=dp), intent(in) :: k(3)
    integer,       intent(in) :: i_arr(:), r_arr(:)

    complex(kind=dp) :: u

    u = r_arr(1)*i_arr(1)*(k(1)**2)*exp(sin(10*k(1)))

  end function calculator_test

end module calculator