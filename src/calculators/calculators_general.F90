module calculators_general

  use utility
  use data_structures
  use local_k_quantities
  use kpath

  implicit none

  public :: bands
  public :: bands_kpath_task_constructor

contains

  !====DEFAULT CALCULATORS====!

  !====BANDS CALCULATOR AND CONSTRUCTOR====!
  function bands(k_data, system, k, error) result(u)
    class(local_k_data), intent(in) :: k_data
    type(sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp)                :: u(k_data%integer_indices(1))

    complex(kind=dp) :: hamiltonian(system%num_bands, system%num_bands), &
                        rot(system%num_bands, system%num_bands)
    real(kind=dp) :: eig(system%num_bands)

    hamiltonian = wannier_hamiltonian(system, k)
    call utility_diagonalize(hamiltonian, system%num_bands, eig, rot, error)
    if (error) then
      write (unit=113, fmt="(a, i3, a)") "Error in function bands when computing the eigenvalues of the Hamiltonian."
      return
    endif
    u = eig
  end function bands
  !==========DEFAULT BANDS KPATH TASK==========!
  subroutine bands_kpath_task_constructor(task, system, Nvec, vec_coord, nkpts)

    type(sys), intent(in)  :: system
    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)

    type(k_path_task), intent(out) :: task

    call kpath_constructor(task=task, name="def_bands", &
                           l_calculator=bands, &
                           Nvec=Nvec, vec_coord=vec_coord, nkpts=nkpts, &
                           N_int_ind=1, int_ind_range=(/system%num_bands/), &
                           N_ext_vars=1)
  end subroutine bands_kpath_task_constructor
  !====END BANDS CALCULATOR AND CONSTRUCTOR====!

end module calculators_general
