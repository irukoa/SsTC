module system

  use utility

  type sys
    character(len=120) :: name
    integer :: num_bands
    real(kind=dp) :: lattice_constant !In A.
    ! real(kind=dp) :: direct_lattice_basis(3, 3), reciprocal_lattice_basis(3, 3) !Relative to a and 1/a, respectively.
    complex(kind=dp), allocatable :: real_space_hamiltonian_elements(:, :, :) !In eV.
    complex(kind=dp), allocatable :: real_space_position_elements(:, :, :, :) !In A.
  end type sys

  type external_vars
    real(kind=dp) :: w_s, w_e !In eV: first \omega and last \omega being considered, respectively.
    integer :: Nw !Number of \omega values being considered.
    real(kind=dp), allocatable :: w !In eV: List of \omega values.
    real(kind=dp) :: t_s, t_e !In s: first t and last t being considered, respectively.
    integer :: Nt !Number of t values being considered.
    real(kind=dp), allocatable :: t !In s: List of t values.
  end type external_vars

  public :: sys
  public :: external_vars

  !Subroutines in the form of WAN HAM.

end module system