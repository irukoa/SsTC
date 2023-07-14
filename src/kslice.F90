module kslice

  use utility
  use data_structures

  implicit none

  type, extends(global_k_data) :: k_slice_task !TODO: MAKE CONSTRUCTOR, SAMPLER AND PRINTER BASED ON KPATH.
    real(kind=dp) :: corner(3), vector(2, 3)
    integer       :: mesh(2)
    complex(kind=dp), allocatable :: kslice_data(:, :, :)!Integer index, continuous index and kpt index respectively.
  end type k_slice_task

end module kslice
