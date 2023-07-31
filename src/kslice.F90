module SsTC_kslice

  use SsTC_utility
  use SsTC_data_structures

  implicit none

  private

  type, extends(SsTC_global_k_data) :: SsTC_kslice_task !TODO: MAKE CONSTRUCTOR, SAMPLER AND PRINTER BASED ON KPATH.
    real(kind=dp) :: corner(3), vector(2, 3)
    integer       :: mesh(2)
    complex(kind=dp), allocatable :: kslice_data(:, :, :)!Integer index, continuous index and kpt index respectively.
  end type SsTC_kslice_task

  public :: SsTC_kslice_task

end module SsTC_kslice
