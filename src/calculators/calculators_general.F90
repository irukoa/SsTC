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
  function bands(k_data, system, k, error) result(u)
    class(local_k_data), intent(in) :: k_data
    type(sys),           intent(in) :: system
    real(kind=dp),       intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp)                :: u(k_data%integer_indices(1))

    complex(kind=dp) :: hamiltonian(system%num_bands, system%num_bands),&
                        rot(system%num_bands, system%num_bands)
    real(kind=dp) :: eig(system%num_bands)

    hamiltonian = wannier_hamiltonian(system, k)
    call utility_diagonalize(hamiltonian, system%num_bands, eig, rot, error)
    if (error) then
      write(unit=113, fmt="(a, i3, a)") "Error in function bands when computing the eigenvalues of the Hamiltonian."
      return
    endif
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

  function velocities(system, HW_a, eig, rot, error) result(v)
    type(sys),          intent(in) :: system
    real(kind=dp),      intent(in) :: eig(system%num_bands)
    complex(kind=dp),   intent(in) :: rot(system%num_bands, system%num_bands)
    complex(kind=dp),   intent(in) :: HW_a(system%num_bands, system%num_bands, 3)

    logical, intent(inout) :: error

    complex(kind=dp) :: v(system%num_bands, system%num_bands, 3)

    integer :: i, n, deg(system%num_bands)
    complex(kind=dp), allocatable :: dummy_rot(:, :)
    real(kind=dp) :: degen_vels(system%num_bands)

    !Get degeneracies.
    deg = utility_get_degen(eig, system%deg_thr)

    !Get velocity matrix as in the case of nondegenerate kpt.
    do i = 1, 3
      v(:, :, i) = matmul(matmul(transpose(conjg(rot)), HW_a(:, :, i)), rot)
    enddo

    !In the case of degenerate bands,
    if (maxval(deg) > 1) then
      !for each coordinate,
      do i = 1, 3
        !initialize degen_vels,
        forall (n=1:system%num_bands) degen_vels(n) = v(n, n, i)
        !and for each eigenvalue,
        do n = 1, system%num_bands
          !check degeneracy.
          if (deg(n) > 1) then
            !In the case of a deg(n) dimensional degenerate subspace,
            allocate (dummy_rot(deg(n), deg(n)))
            !diagonalize the subspace and overwrite to degen_vels.
            call utility_diagonalize(v(n:n + deg(n) - 1, n:n + deg(n) - 1, i), &
                                     deg(n), degen_vels(n:n + deg(n) - 1), dummy_rot, error)
            if (error) then
              write(unit=113, fmt="(a)") "Error in function velocities when computing the eigenvalues of the degenerate subspace."
              return
            endif
            deallocate(dummy_rot)
          endif
        enddo
        !Overwrite to v.
        forall (n=1:system%num_bands) v(n, n, i) = degen_vels(n)
      enddo
    endif

  end function velocities

  function inverse_effective_mass(system, HW_a_b, HW_a, eig, rot, error) result(mu)
    type(sys),          intent(in) :: system
    real(kind=dp),      intent(in) :: eig(system%num_bands)
    complex(kind=dp),   intent(in) :: rot(system%num_bands, system%num_bands)
    complex(kind=dp),   intent(in) :: HW_a_b(system%num_bands, system%num_bands, 3, 3)
    complex(kind=dp),   intent(in) :: HW_a(system%num_bands, system%num_bands, 3)

    logical, intent(inout) :: error

    complex(kind=dp) :: mu(system%num_bands, system%num_bands, 3, 3), &
                        rot_H_a(system%num_bands, system%num_bands, 3), &
                        D_a(system%num_bands, system%num_bands, 3), &
                        rot_H_a_b(system%num_bands, system%num_bands, 3, 3)

    integer :: i, j, ij, n, deg(system%num_bands)
    complex(kind=dp), allocatable :: dummy_rot(:, :)
    real(kind=dp) :: degen_mass(system%num_bands)

    !Get degeneracies.
    deg = utility_get_degen(eig, system%deg_thr)

    !Get non-abelian D matrix.
    D_a = non_abelian_d(system, eig, rot, HW_a)

    !Get rotated quantities.
    do i = 1, 3
      rot_H_a(:, :, i) = matmul(matmul(transpose(conjg(rot)), HW_a(:, :, i)), rot)
      do j = 1, 3
        rot_H_a_b(:, :, i, j) = matmul(matmul(transpose(conjg(rot)), HW_a_b(:, :, i, j)), rot)
      enddo
    enddo

    !Get inverse effective mass matrix as in the case of nondegenerate kpt.
    do ij = 1, 6
      i = alpha_S(ij)
      j = beta_S(ij)
      !Eq. 28 of 10.1103/PhysRevB.75.195121
      mu(:, :, i, j) = matmul(rot_H_a(:, :, i), D_a(:, :, j))
      mu(:, :, i, j) = rot_H_a_b(:, :, i, j) + mu(:, :, i, j) + transpose(conjg(mu(:, :, i, j)))
      !Symmetrization.
      mu(:, :, j, i) = mu(:, :, i, j)
    enddo

    !In the case of degenerate bands,
    if (maxval(deg) > 1) then
      !for each coordinate,
      do ij = 1, 6
        i = alpha_S(ij)
        j = beta_S(ij)
        !initialize degen_mass,
        forall (n=1:system%num_bands) degen_mass(n) = mu(n, n, i, j)
        !and for each eigenvalue,
        do n = 1, system%num_bands
          !check degeneracy.
          if (deg(n) > 1) then
            !In the case of a deg(n) dimensional degenerate subspace,
            allocate (dummy_rot(deg(n), deg(n)))
            !diagonalize the subspace and overwrite to degen_mass.
            call utility_diagonalize(mu(n:n + deg(n) - 1, n:n + deg(n) - 1, i, j), &
                                     deg(n), degen_mass(n:n + deg(n) - 1), dummy_rot, error)
            if (error) then
              write(unit=113, fmt="(a)") "Error in function inverse_effective_mass when computing the eigenvalues of the degenerate subspace."
              return
            endif
            deallocate(dummy_rot)
          endif
        enddo
        !Overwrite to mu.
        forall (n=1:system%num_bands) mu(n, n, i, j) = degen_mass(n)
        !Symmetrization.
        mu(:, :, j, i) = mu(:, :, i, j)
      enddo
    endif

  end function inverse_effective_mass

end module calculators_general