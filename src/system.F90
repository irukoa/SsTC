module system

  use utility

  type sys
    character(len=120) :: name
    integer :: num_bands
    real(kind=dp) :: lattice_constant !In A.
    real(kind=dp) :: direct_lattice_basis(3, 3), reciprocal_lattice_basis(3, 3) !Relative to a and 1/a, respectively.
    complex(kind=dp), allocatable :: real_space_hamiltonian_elements(:, :, :) !In eV.
    complex(kind=dp), allocatable :: real_space_position_elements(:, :, :, :) !In A.
  end type sys

  public :: sys

  !Subroutines in the form of WAN HAM.

end module system