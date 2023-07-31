module SsTC_calculators_general

  use SsTC_utility
  use SsTC_data_structures
  use SsTC_local_k_quantities
  use SsTC_kpath

  implicit none

  private

  public :: SsTC_bands_kpath_task_constructor
  public :: SsTC_bands

contains

  !====DEFAULT CALCULATORS====!

  !====BANDS CALCULATOR AND CONSTRUCTOR====!
  function SsTC_bands(k_data, system, k, error) result(u)
    class(SsTC_local_k_data), intent(in) :: k_data
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp)                :: u(k_data%integer_indices(1))

    complex(kind=dp) :: hamiltonian(system%num_bands, system%num_bands), &
                        rot(system%num_bands, system%num_bands)
    real(kind=dp) :: eig(system%num_bands)

    hamiltonian = SsTC_wannier_hamiltonian(system, k)
    call SsTC_utility_diagonalize(hamiltonian, system%num_bands, eig, rot, error)
    if (error) then
      write (unit=stderr, fmt="(a)") "Error in function bands when computing the eigenvalues of the Hamiltonian."
      return
    endif
    u = eig
  end function SsTC_bands
  !==========DEFAULT BANDS KPATH TASK==========!
  subroutine SsTC_bands_kpath_task_constructor(task, system, Nvec, vec_coord, nkpts)

    type(SsTC_sys), intent(in)  :: system
    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)

    type(SsTC_kpath_task), intent(out) :: task

    call SsTC_kpath_constructor(task=task, name="def_bands", &
                                l_calculator=SsTC_bands, &
                                Nvec=Nvec, vec_coord=vec_coord, nkpts=nkpts, &
                                N_int_ind=1, int_ind_range=(/system%num_bands/), &
                                N_ext_vars=1)
  end subroutine SsTC_bands_kpath_task_constructor
  !====END BANDS CALCULATOR AND CONSTRUCTOR====!

end module SsTC_calculators_general
