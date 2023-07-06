module calculators_general

  use utility
  use data_structures
  use local_k_quantities
  use kpath

  implicit none

  public :: bands
  public :: bands_kpath_task_constructor
  public :: wannier_hamiltonian
  public :: wannier_berry_connection
  public :: wannier_dhamiltonian_dk
  public :: wannier_d2hamiltonian_dk2
  public :: hamiltonian_occ_matrix
  public :: non_abelian_d

  contains

  !====DEFAULT CALCULATORS====!

  !====BANDS CALCULATOR AND CONSTRUCTOR====!
  function bands(k_data, system, k) result(u)
    class(local_k_data), intent(in) :: k_data
    type(sys),           intent(in) :: system
    real(kind=dp),       intent(in) :: k(3)

    complex(kind=dp)                :: u(k_data%integer_indices(1))

    complex(kind=dp) :: hamiltonian(system%num_bands, system%num_bands),&
                        rot(system%num_bands, system%num_bands)
    real(kind=dp) :: eig(system%num_bands)

    logical :: error
    character(len=120) :: errormsg

    hamiltonian = wannier_hamiltonian(system, k)
    call utility_diagonalize(hamiltonian, system%num_bands, eig, rot, error, errormsg)
    u = eig
  end function bands
  !==========DEFAULT BANDS KPATH TASK==========!
  function bands_kpath_task_constructor(system, Nvec, vec_coord, nkpts) result(default_bands_kpath_task)

    type(sys), intent(in)  :: system
    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)

    type(k_path_task) :: default_bands_kpath_task

    default_bands_kpath_task = kpath_constructor(name = "def_bands", &
                                                 l_calculator = bands, &
                                                 Nvec = Nvec, vec_coord = vec_coord, nkpts = nkpts, &
                                                 N_int_ind = 1, int_ind_range = (/system%num_bands/), &
                                                 N_ext_vars = 1)
  end function bands_kpath_task_constructor
  !====END BANDS CALCULATOR AND CONSTRUCTOR====!


  !====AUXILIARY ROUTINES====!

  function wannier_hamiltonian(system, k) result(HW)
    !Output: Wannier basis Hamiltonian.
    !1st and 2nd indexes: bands.
    type(sys),          intent(in)  :: system
    real(kind=dp),      intent(in)  :: k(3)

    complex(kind=dp) :: HW(system%num_bands, system%num_bands)

    type(local_k_data), allocatable :: H_get(:)
    integer, allocatable :: i_arr(:)
    integer :: i_mem

    !Get Hamiltonian in memory layout.
    call get_hamiltonian(system, k, H_get)

    allocate(i_arr(size(H_get(1)%integer_indices)))
    do i_mem = 1, product(H_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = integer_memory_element_to_array_element(H_get(1), i_mem)
      HW(i_arr(1), i_arr(2)) = H_get(1)%k_data(i_mem)
    enddo
    deallocate(H_get, i_arr)

  end function wannier_hamiltonian

  function wannier_berry_connection(system, k) result(AW)
    !Output: Wannier basis Berry connection.
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    type(sys),          intent(in)  :: system
    real(kind=dp),      intent(in)  :: k(3)

    complex(kind=dp) :: AW(system%num_bands, system%num_bands, 3)

    type(local_k_data), allocatable :: A_get(:)
    integer, allocatable :: i_arr(:)
    integer :: i_mem

    !Get Berry connection in memory layout.
    call get_position(system, k, A_get)

    allocate(i_arr(size(A_get(1)%integer_indices)))
    do i_mem = 1, product(A_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = integer_memory_element_to_array_element(A_get(1), i_mem)
      AW(i_arr(1), i_arr(2), i_arr(3)) = A_get(1)%k_data(i_mem)
    enddo
    deallocate(A_get, i_arr)

  end function wannier_berry_connection

  function wannier_dhamiltonian_dk(system, k) result(DHW)
    !Output: 1st k-derivative of the Wannier Hamiltonian.
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    type(sys),          intent(in)  :: system
    real(kind=dp),      intent(in)  :: k(3)

    complex(kind=dp) :: DHW(system%num_bands, system%num_bands, 3)

    type(local_k_data), allocatable :: H_get(:)
    integer, allocatable :: i_arr(:)
    integer :: i_mem

    !Get Hamiltonian in memory layout.
    call get_hamiltonian(system, k, H_get, Nder_i = 1)

    allocate(i_arr(size(H_get(1)%integer_indices)))
    do i_mem = 1, product(H_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = integer_memory_element_to_array_element(H_get(1), i_mem)
      DHW(i_arr(1), i_arr(2), i_arr(3)) = H_get(1)%k_data(i_mem)
    enddo
    deallocate(H_get, i_arr)

  end function wannier_dhamiltonian_dk

  function wannier_d2hamiltonian_dk2(system, k) result(DDHW)
    !Output: 2nd k-derivative of the Wannier Hamiltonian.
    !1st and 2nd indexes: bands, 3rd and 4th indexes: cartesian comp.
    type(sys),          intent(in)  :: system
    real(kind=dp),      intent(in)  :: k(3)

    complex(kind=dp) :: DDHW(system%num_bands, system%num_bands, 3, 3)

    type(local_k_data), allocatable :: H_get(:)
    integer, allocatable :: i_arr(:)
    integer :: i_mem

    !Get Hamiltonian in memory layout.
    call get_hamiltonian(system, k, H_get, Nder_i = 2)

    allocate(i_arr(size(H_get(1)%integer_indices)))
    do i_mem = 1, product(H_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = integer_memory_element_to_array_element(H_get(1), i_mem)
      DDHW(i_arr(1), i_arr(2), i_arr(3), i_arr(4)) = H_get(1)%k_data(i_mem)
    enddo
    deallocate(H_get, i_arr)

  end function wannier_d2hamiltonian_dk2

  function hamiltonian_occ_matrix(system, eig) result(rho)
    !Output: Fermi occupation matrix in the Hamiltonian basis.
    !1st and 2nd indexes: bands.
    type(sys),          intent(in)  :: system
    real(kind=dp),      intent(in)  :: eig(system%num_bands)

    complex(kind=dp) :: rho(system%num_bands, system%num_bands)

    integer :: ibnd

    rho = cmplx_0

    do ibnd = 1, system%num_bands
      if (eig(ibnd) .le. system%e_fermi) rho(ibnd, ibnd) = 1.0_dp
    enddo

  end function hamiltonian_occ_matrix

  function non_abelian_d(system, eig, rot, HW_a) result(DH)
    !Output: Vector valued matrix D in Eq. (32) in 
    !10.1103/PhysRevB.75.195121 .
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    type(sys),          intent(in) :: system
    real(kind=dp),      intent(in) :: eig(system%num_bands)
    complex(kind=dp),   intent(in) :: rot(system%num_bands, system%num_bands)
    complex(kind=dp),   intent(in) :: HW_a(system%num_bands, system%num_bands, 3)

    complex(kind=dp) :: DH(system%num_bands, system%num_bands, 3)

    integer :: i, n, m

    do i = 1, 3
      DH(:, :, i) = matmul(matmul(transpose(conjg(rot)), HW_a(:, :, i)), rot)
      do n = 1, system%num_bands
        do m = 1, system%num_bands
          if (abs(eig(n)-eig(m)) < system%deg_thr) then
            DH(n, m, i) = cmplx_0
          else
            DH(n, m, i) = DH(n, m, i)*((eig(m) - eig(n))/((eig(m) - eig(n))**2 + (system%deg_offset)**2))
          endif
        enddo
      enddo
    enddo
  end function non_abelian_d

end module calculators_general