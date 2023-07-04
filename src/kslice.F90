module kslice

  use utility
  use data_structures

  implicit none

  type, extends(local_k_data) :: k_slice !TODO: MAKE CONSTRUCTOR.
    real(kind=dp) :: corner(3), vector(2, 3)
    integer       :: mesh(2)
  end type k_slice

end module kslice