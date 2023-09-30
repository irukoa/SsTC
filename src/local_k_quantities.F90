#include "cond_comp.h"
module SsTC_local_k_quantities

  USE OMP_LIB
  USE MPI_F08

  use SsTC_utility
  use SsTC_mpi_comms
  use SsTC_data_structures

  implicit none

  private

  public :: SsTC_wannier_hamiltonian
  public :: SsTC_wannier_berry_connection
  public :: SsTC_wannier_dhamiltonian_dk
  public :: SsTC_wannier_d2hamiltonian_dk2
  public :: SsTC_wannier_dberry_connection_dk
  public :: SsTC_wannier_momentum
  public :: SsTC_hamiltonian_occ_matrix
  public :: SsTC_non_abelian_d

  public :: SsTC_deleig
  public :: SsTC_inverse_effective_mass

  public :: SsTC_get_hamiltonian
  public :: SsTC_get_position

contains

  !====AUXILIARY ROUTINES====!

  function SsTC_wannier_hamiltonian(system, k) result(HW)

    implicit none

    !Output: Wannier basis Hamiltonian.
    !1st and 2nd indexes: bands.
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: k(3)

    complex(kind=dp) :: HW(system%num_bands, system%num_bands)

    type(SsTC_local_k_data), allocatable :: H_get(:)
    integer, allocatable                 :: i_arr(:)
    integer                              :: i_mem

    !Get Hamiltonian in memory layout.
    call SsTC_get_hamiltonian(system, k, H_get)

    allocate (i_arr(size(H_get(1)%integer_indices)))
    do i_mem = 1, product(H_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = SsTC_integer_memory_element_to_array_element(H_get(1), i_mem)
      HW(i_arr(1), i_arr(2)) = H_get(1)%k_data(i_mem)
    enddo
    deallocate (H_get, i_arr)

  end function SsTC_wannier_hamiltonian

  function SsTC_wannier_berry_connection(system, k) result(AW)

    implicit none

    !Output: Wannier basis Berry connection.
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: k(3)

    complex(kind=dp) :: AW(system%num_bands, system%num_bands, 3)

    type(SsTC_local_k_data), allocatable :: A_get(:)
    integer, allocatable                 :: i_arr(:)
    integer                              :: i_mem

    !Get Berry connection in memory layout.
    call SsTC_get_position(system, k, A_get)

    allocate (i_arr(size(A_get(1)%integer_indices)))
    do i_mem = 1, product(A_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = SsTC_integer_memory_element_to_array_element(A_get(1), i_mem)
      AW(i_arr(1), i_arr(2), i_arr(3)) = A_get(1)%k_data(i_mem)
    enddo
    deallocate (A_get, i_arr)

  end function SsTC_wannier_berry_connection

  function SsTC_wannier_dhamiltonian_dk(system, k) result(DHW)

    implicit none

    !Output: 1st k-derivative of the Wannier Hamiltonian.
    !1st and 2nd indexes: bands, 3rd index: derivative comp.
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: k(3)

    complex(kind=dp) :: DHW(system%num_bands, system%num_bands, 3)

    type(SsTC_local_k_data), allocatable :: H_get(:)
    integer, allocatable                 :: i_arr(:)
    integer                              :: i_mem

    !Get Hamiltonian in memory layout.
    call SsTC_get_hamiltonian(system, k, H_get, Nder_i=1)

    allocate (i_arr(size(H_get(1)%integer_indices)))
    do i_mem = 1, product(H_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = SsTC_integer_memory_element_to_array_element(H_get(1), i_mem)
      DHW(i_arr(1), i_arr(2), i_arr(3)) = H_get(1)%k_data(i_mem)
    enddo
    deallocate (H_get, i_arr)

  end function SsTC_wannier_dhamiltonian_dk

  function SsTC_wannier_d2hamiltonian_dk2(system, k) result(DDHW)

    implicit none

    !Output: 2nd k-derivative of the Wannier Hamiltonian.
    !1st and 2nd indexes: bands,
    !3rd and 4th indexes: derivative directions.
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: k(3)

    complex(kind=dp) :: DDHW(system%num_bands, system%num_bands, 3, 3)

    type(SsTC_local_k_data), allocatable :: H_get(:)
    integer, allocatable                 :: i_arr(:)
    integer                              :: i_mem

    !Get Hamiltonian in memory layout.
    call SsTC_get_hamiltonian(system, k, H_get, Nder_i=2)

    allocate (i_arr(size(H_get(1)%integer_indices)))
    do i_mem = 1, product(H_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = SsTC_integer_memory_element_to_array_element(H_get(1), i_mem)
      DDHW(i_arr(1), i_arr(2), i_arr(3), i_arr(4)) = H_get(1)%k_data(i_mem)
    enddo
    deallocate (H_get, i_arr)

  end function SsTC_wannier_d2hamiltonian_dk2

  function SsTC_wannier_dberry_connection_dk(system, k) result(DAW)

    implicit none

    !Output: 1st k-derivative of the Wannier Berry connection.
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    !4th index : derivative direction
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: k(3)

    complex(kind=dp) :: DAW(system%num_bands, system%num_bands, 3, 3)

    type(SsTC_local_k_data), allocatable :: A_get(:)
    integer, allocatable                 :: i_arr(:)
    integer                              :: i_mem

    !Get Hamiltonian in memory layout.
    call SsTC_get_position(system, k, A_get, Nder_i=1)

    allocate (i_arr(size(A_get(1)%integer_indices)))
    do i_mem = 1, product(A_get(1)%integer_indices)
      !Pass to array layout.
      i_arr = SsTC_integer_memory_element_to_array_element(A_get(1), i_mem)
      DAW(i_arr(1), i_arr(2), i_arr(3), i_arr(4)) = A_get(1)%k_data(i_mem)
    enddo
    deallocate (A_get, i_arr)

  end function SsTC_wannier_dberry_connection_dk

  function SsTC_wannier_momentum(system, k) result(PW)

    implicit none

    !Output: Full momentum operator p_nm^a(k) in the Wannier gauge:
    !p_nm^a(k) = \partial H_0_nm / \partial k^a - i*[A^a, H_0]_nm,
    !i.e., the action of the covariant derivative on H_0.
    !the units are eV*A.
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.

    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: k(3)

    complex(kind=dp) :: PW(system%num_bands, system%num_bands, 3)

    complex(kind=dp) :: HW(system%num_bands, system%num_bands), &
                        DHW(system%num_bands, system%num_bands, 3), &
                        AW(system%num_bands, system%num_bands, 3), &
                        commutator(system%num_bands, system%num_bands, 3)

    integer :: i

    !Get Hamiltonian,
    HW = SsTC_wannier_hamiltonian(system, k)
    !it's k-derivative,
    DHW = SsTC_wannier_dhamiltonian_dk(system, k)
    !and the Berry connection.
    AW = SsTC_wannier_berry_connection(system, k)

    !Get commutator and set result.
    do i = 1, 3
      commutator(:, :, i) = matmul(AW(:, :, i), HW) - matmul(HW, AW(:, :, i))
      PW(:, :, i) = DHW(:, :, i) - cmplx_i*commutator(:, :, i)
    enddo

  end function SsTC_wannier_momentum

  function SsTC_hamiltonian_occ_matrix(system, eig) result(rho)

    implicit none

    !Output: Fermi occupation matrix in the Hamiltonian basis.
    !1st and 2nd indexes: bands.
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in)  :: eig(system%num_bands)

    complex(kind=dp) :: rho(system%num_bands, system%num_bands)

    integer :: ibnd

    rho = cmplx_0

    do ibnd = 1, system%num_bands
      if (eig(ibnd) .le. system%e_fermi) rho(ibnd, ibnd) = 1.0_dp
    enddo

  end function SsTC_hamiltonian_occ_matrix

  function SsTC_non_abelian_d(system, eig, rot, HW_a) result(DH)

    implicit none

    !Output: Vector valued matrix D in Eq. (32) in
    !10.1103/PhysRevB.75.195121 .
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    type(SsTC_sys), intent(in)   :: system
    real(kind=dp), intent(in)    :: eig(system%num_bands)
    complex(kind=dp), intent(in) :: rot(system%num_bands, system%num_bands)
    complex(kind=dp), intent(in) :: HW_a(system%num_bands, system%num_bands, 3)

    complex(kind=dp) :: DH(system%num_bands, system%num_bands, 3)

    integer :: i, n, m

    do i = 1, 3
      DH(:, :, i) = matmul(matmul(transpose(conjg(rot)), HW_a(:, :, i)), rot)
      do n = 1, system%num_bands
        do m = 1, system%num_bands
          if (abs(eig(n) - eig(m)) < system%deg_thr) then
            DH(n, m, i) = cmplx_0
          else
            DH(n, m, i) = DH(n, m, i)*((eig(m) - eig(n))/((eig(m) - eig(n))**2 + (system%deg_offset)**2))
          endif
        enddo
      enddo
    enddo
  end function SsTC_non_abelian_d

  function SsTC_deleig(system, HW_a, eig, rot, error) result(v)

    implicit none

    !Output: Velocities v_{nm, a} in Eq. (23) in
    !10.1103/PhysRevB.75.195121 .
    !1st and 2nd indexes: bands, 3rd index: cartesian comp.
    type(SsTC_sys), intent(in)   :: system
    real(kind=dp), intent(in)    :: eig(system%num_bands)
    complex(kind=dp), intent(in) :: rot(system%num_bands, system%num_bands)
    complex(kind=dp), intent(in) :: HW_a(system%num_bands, system%num_bands, 3)
    logical, intent(inout)       :: error

    complex(kind=dp) :: v(system%num_bands, system%num_bands, 3)

    integer :: i, n, deg(system%num_bands)
    complex(kind=dp), allocatable :: dummy_rot(:, :)
    real(kind=dp) :: degen_vels(system%num_bands)

    !Get degeneracies.
    deg = SsTC_utility_get_degen(eig, system%deg_thr)

    !Get velocity matrix as in the case of nondegenerate kpt.
    do i = 1, 3
      v(:, :, i) = matmul(matmul(transpose(conjg(rot)), HW_a(:, :, i)), rot)
    enddo

    !In the case of degenerate bands,
    if (maxval(deg) > 1) then
      !for each coordinate,
      do i = 1, 3
        !initialize degen_vels,
        do n = 1, system%num_bands
          degen_vels(n) = real(v(n, n, i), dp)
        enddo
        !and for each eigenvalue,
        do n = 1, system%num_bands
          !check degeneracy.
          if (deg(n) > 1) then
            !In the case of a deg(n) dimensional degenerate subspace,
            allocate (dummy_rot(deg(n), deg(n)))
            !diagonalize the subspace and overwrite to degen_vels.
            call SsTC_utility_diagonalize(v(n:n + deg(n) - 1, n:n + deg(n) - 1, i), &
                                          deg(n), degen_vels(n:n + deg(n) - 1), dummy_rot, error)
            if (error) return
            deallocate (dummy_rot)
          endif
        enddo
        !Overwrite to v.
        do n = 1, system%num_bands
          v(n, n, i) = degen_vels(n)
        enddo
      enddo
    endif

  end function SsTC_deleig

  function SsTC_inverse_effective_mass(system, HW_a_b, HW_a, eig, rot, error) result(mu)

    implicit none

    !Output: Inverse effective mass \mu_{nm, ab} in Eq. (24) in
    !10.1103/PhysRevB.75.195121 .
    !1st and 2nd indexes: bands, 3rd and 4th index: cartesian comp.
    type(SsTC_sys), intent(in)   :: system
    real(kind=dp), intent(in)    :: eig(system%num_bands)
    complex(kind=dp), intent(in) :: rot(system%num_bands, system%num_bands)
    complex(kind=dp), intent(in) :: HW_a_b(system%num_bands, system%num_bands, 3, 3)
    complex(kind=dp), intent(in) :: HW_a(system%num_bands, system%num_bands, 3)
    logical, intent(inout)       :: error

    complex(kind=dp) :: mu(system%num_bands, system%num_bands, 3, 3)

    complex(kind=dp) :: rot_H_a(system%num_bands, system%num_bands, 3), &
                        rot_H_a_b(system%num_bands, system%num_bands, 3, 3), &
                        D_a(system%num_bands, system%num_bands, 3)

    integer                       :: i, j, ij, &
                                     n, &
                                     deg(system%num_bands)
    complex(kind=dp), allocatable :: dummy_rot(:, :)
    real(kind=dp)                 :: degen_mass(system%num_bands)

    !Get degeneracies.
    deg = SsTC_utility_get_degen(eig, system%deg_thr)

    !Get non-abelian D matrix.
    D_a = SsTC_non_abelian_d(system, eig, rot, HW_a)

    !Get rotated quantities.
    do i = 1, 3
      rot_H_a(:, :, i) = matmul(matmul(transpose(conjg(rot)), HW_a(:, :, i)), rot)
      do j = 1, 3
        rot_H_a_b(:, :, i, j) = matmul(matmul(transpose(conjg(rot)), HW_a_b(:, :, i, j)), rot)
      enddo
    enddo

    !Get inverse effective mass matrix as in the case of nondegenerate kpt.
    do ij = 1, 6
      i = SsTC_alpha_S(ij)
      j = SsTC_beta_S(ij)
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
        i = SsTC_alpha_S(ij)
        j = SsTC_beta_S(ij)
        !initialize degen_mass,
        do n = 1, system%num_bands
          degen_mass(n) = real(mu(n, n, i, j), dp)
        enddo
        !and for each eigenvalue,
        do n = 1, system%num_bands
          !check degeneracy.
          if (deg(n) > 1) then
            !In the case of a deg(n) dimensional degenerate subspace,
            allocate (dummy_rot(deg(n), deg(n)))
            !diagonalize the subspace and overwrite to degen_mass.
            call SsTC_utility_diagonalize(mu(n:n + deg(n) - 1, n:n + deg(n) - 1, i, j), &
                                          deg(n), degen_mass(n:n + deg(n) - 1), dummy_rot, error)
            if (error) return
            deallocate (dummy_rot)
          endif
        enddo
        !Overwrite to mu.
        do n = 1, system%num_bands
          mu(n, n, i, j) = degen_mass(n)
        enddo
        !Symmetrization.
        mu(:, :, j, i) = mu(:, :, i, j)
      enddo
    endif

  end function SsTC_inverse_effective_mass

  !====CORE ROUTINES====!

  subroutine SsTC_get_hamiltonian(system, k, H, Nder_i, only_i)

    implicit none

    type(SsTC_sys), intent(in)                        :: system
    real(kind=dp), intent(in)                         :: k(3)
    type(SsTC_local_k_data), allocatable, intent(out) :: H(:)
    integer, optional                                 :: Nder_i !Order of the derivative. Nder = 0 means normal Hamiltonian,
    !Nder = 1, \partial H/\partial k^i...
    logical, optional                                 :: only_i !Determine, in the case of Nder > 0,
    !if all derivatives i with i < Nder are requested.

    integer :: Nder
    logical :: only

    integer              :: i, j, &
                            i_mem, &
                            irpts, ivec
    integer, allocatable ::  i_arr(:)
    real(kind=dp)        :: kdotr, prod_R, vec(3)
    complex(kind=dp)     :: temp_res

    if (.not. (present(Nder_i))) then
      Nder = 0
    else
      Nder = Nder_i
    endif

    if (.not. (present(only_i))) then
      only = .true.
    else
      only = only_i
    endif

    if (only) allocate (H(1)) !Only compute the Nth Hamiltonian derivative.
    if (.not. only) allocate (H(Nder + 1)) !Compute from the 0th to the Nth Hamiltonian derivative.

    if (.not. only) then
      do i = 1, Nder + 1
        allocate (H(i)%integer_indices(2 + i - 1)) !2 band indices and i - 1derivative indices.
        H(i)%integer_indices(1) = system%num_bands
        H(i)%integer_indices(2) = system%num_bands
        if (i /= 1) then
          do j = 1, i - 1
            H(i)%integer_indices(2 + j) = 3
          enddo
        endif
      enddo
    else
      allocate (H(1)%integer_indices(2 + Nder)) !2 band indices and Nder derivative indices.
      H(1)%integer_indices(1) = system%num_bands
      H(1)%integer_indices(2) = system%num_bands
      if (Nder /= 0) then
        do i = 1, Nder
          H(1)%integer_indices(2 + i) = 3
        enddo
      endif
    endif

    if (only) then !Case of only a requested derivative.

      allocate (i_arr(size(H(1)%integer_indices))) !Array indices storage.
      allocate (H(1)%k_data(product(H(1)%integer_indices))) !Data storage.
      H(1)%k_data = cmplx_0 !Initialize.

      do i_mem = 1, product(H(1)%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(H(1), i_mem) !Get array index.

        temp_res = cmplx_0
        !_OMPTGT_(PARALLEL PRIVATE (kdotr, vec, prod_R) REDUCTION (+: temp_res))
        !_OMPTGT_(DO)
        do irpts = 1, system%num_R_points !For each point in the Bravais lattice.

          !Compute factor appearing in the exponential (k is in coords relative to recip. lattice vectors).
          kdotr = 2.0_dp*pi*dot_product(system%R_point(irpts, :), k)

          !Compute Bravais lattice vector for label irpts.
          vec = 0.0_dp
          do ivec = 1, 3
            vec = vec + system%R_point(irpts, ivec)*system%direct_lattice_basis(ivec, :)
          enddo

          !In the case of derivatives compute the prefactor involving a product of Bravais lattice vector components.
          prod_R = 1.0_dp
          if (Nder /= 0) then
            do j = 1, Nder
              prod_R = prod_R*vec(i_arr(j + 2))
            enddo
          endif

          temp_res = temp_res + & !Compute sum.
                     ((cmplx_i)**Nder)*(prod_R)*exp(cmplx_i*kdotr)* &
                     system%real_space_hamiltonian_elements(i_arr(1), i_arr(2), irpts) !&
          !/real(system%deg_R_point(irpts), dp)
          !TODO: Recheck the deg_R_points division.
        enddo!irpts
        !_OMPTGT_(END DO)
        !_OMPTGT_(END PARALLEL)
        H(1)%k_data(i_mem) = temp_res
      enddo!i_mem

      deallocate (i_arr)

    else !Case of requestad all derivatives.
      do i = 1, Nder + 1

        allocate (i_arr(size(H(i)%integer_indices))) !Array indices storage.
        allocate (H(i)%k_data(product(H(i)%integer_indices))) !Data storage.
        H(i)%k_data = cmplx_0 !Initialize.

        do i_mem = 1, product(H(i)%integer_indices) !For each integer index.

          i_arr = SsTC_integer_memory_element_to_array_element(H(i), i_mem) !Get array index.

          temp_res = cmplx_0
          !_OMPTGT_(PARALLEL PRIVATE (kdotr, vec, prod_R) REDUCTION (+: temp_res))
          !_OMPTGT_(DO)
          do irpts = 1, system%num_R_points !For each point in the Bravais lattice.

            !Compute factor appearing in the exponential (k is in coords relative to recip. lattice vectors).
            kdotr = 2.0_dp*pi*dot_product(system%R_point(irpts, :), k)

            !Compute Bravais lattice vector for label irpts.
            vec = 0.0_dp
            do ivec = 1, 3
              vec = vec + system%R_point(irpts, ivec)*system%direct_lattice_basis(ivec, :)
            enddo

            !In the case of derivatives compute the prefactor involving a product of Bravais lattice vector components.
            prod_R = 1.0_dp
            if (i /= 1) then
              do j = 1, i - 1
                prod_R = prod_R*vec(i_arr(j + 2))
              enddo
            endif

            temp_res = temp_res + & !Compute sum.
                       ((cmplx_i)**(i - 1))*(prod_R)* &
                       exp(cmplx_i*kdotr)*system%real_space_hamiltonian_elements(i_arr(1), i_arr(2), irpts) !&
            !/real(system%deg_R_point(irpts), dp)
            !TODO: Recheck the deg_R_points division.
          enddo!irpts
          !_OMPTGT_(END DO)
          !_OMPTGT_(END PARALLEL)
          H(i)%k_data(i_mem) = temp_res
        enddo!i_mem

        deallocate (i_arr)

      enddo
    endif

  end subroutine SsTC_get_hamiltonian

  subroutine SsTC_get_position(system, k, A, Nder_i, only_i)

    implicit none

    type(SsTC_sys), intent(in)                        :: system
    real(kind=dp), intent(in)                         :: k(3)
    type(SsTC_local_k_data), allocatable, intent(out) :: A(:)
    integer, optional                                 :: Nder_i !Order of the derivative. Nder = 0 means normal Berry connection,
    !Nder = 1, \partial A^i/\partial k^j...
    logical, optional                                 :: only_i !Determine, in the case of Nder > 0,
    !if all derivatives i with i < Nder are requested.

    integer :: Nder
    logical :: only

    integer              :: i, j, &
                            i_mem, &
                            irpts, ivec
    integer, allocatable ::  i_arr(:)
    real(kind=dp)        :: kdotr, prod_R, vec(3)
    complex(kind=dp)     :: temp_res

    if (.not. (present(Nder_i))) then
      Nder = 0
    else
      Nder = Nder_i
    endif

    if (.not. (present(only_i))) then
      only = .true.
    else
      only = only_i
    endif

    if (only) allocate (A(1)) !Only compute the Nth position operator's derivative.
    if (.not. only) allocate (A(Nder + 1)) !Compute from the 0th to the Nth position operator's derivative.

    if (.not. only) then
      do i = 1, Nder + 1
        allocate (A(i)%integer_indices(3 + i - 1)) !2 band indices, 1 component index and i - 1 derivative indices.
        A(i)%integer_indices(1) = system%num_bands
        A(i)%integer_indices(2) = system%num_bands
        A(i)%integer_indices(3) = 3
        if (i /= 1) then
          do j = 1, i - 1
            A(i)%integer_indices(3 + j) = 3
          enddo
        endif
      enddo
    else
      allocate (A(1)%integer_indices(3 + Nder)) !2 band indices, 1 component index and i - 1 derivative indices.
      A(1)%integer_indices(1) = system%num_bands
      A(1)%integer_indices(2) = system%num_bands
      A(1)%integer_indices(3) = 3
      if (Nder /= 0) then
        do i = 1, Nder
          A(1)%integer_indices(3 + i) = 3
        enddo
      endif
    endif

    if (only) then !Case of only a requested derivative.

      allocate (i_arr(size(A(1)%integer_indices))) !Array indices storage.
      allocate (A(1)%k_data(product(A(1)%integer_indices))) !Data storage.
      A(1)%k_data = cmplx_0 !Initialize.

      do i_mem = 1, product(A(1)%integer_indices) !For each integer index.

        i_arr = SsTC_integer_memory_element_to_array_element(A(1), i_mem) !Get array index.

        temp_res = cmplx_0
        !_OMPTGT_(PARALLEL PRIVATE (kdotr, vec, prod_R) REDUCTION (+: temp_res))
        !_OMPTGT_(DO)
        do irpts = 1, system%num_R_points !For each point in the Bravais lattice.

          !Compute factor appearing in the exponential (k is in coords relative to recip. lattice vectors).
          kdotr = 2.0_dp*pi*dot_product(system%R_point(irpts, :), k)

          !Compute Bravais lattice vector for label irpts.
          vec = 0.0_dp
          do ivec = 1, 3
            vec = vec + system%R_point(irpts, ivec)*system%direct_lattice_basis(ivec, :)
          enddo

          !In the case of derivatives compute the prefactor involving a product of Bravais lattice vector components.
          prod_R = 1.0_dp
          if (Nder /= 0) then
            do j = 1, Nder
              prod_R = prod_R*vec(i_arr(j + 3))
            enddo
          endif

          temp_res = temp_res + & !Compute sum.
                     ((cmplx_i)**Nder)*(prod_R)* &
                     exp(cmplx_i*kdotr)*system%real_space_position_elements(i_arr(1), i_arr(2), i_arr(3), irpts) !&
          !/real(system%deg_R_point(irpts), dp)
          !TODO: Recheck the deg_R_points division.
        enddo!irpts
        !_OMPTGT_(END DO)
        !_OMPTGT_(END PARALLEL)
        A(1)%k_data(i_mem) = temp_res
      enddo!i_mem

      deallocate (i_arr)

    else !Case of requestad all derivatives.
      do i = 1, Nder + 1

        allocate (i_arr(size(A(i)%integer_indices))) !Array indices storage.
        allocate (A(i)%k_data(product(A(i)%integer_indices))) !Data storage.
        A(i)%k_data = cmplx_0 !Initialize.

        do i_mem = 1, product(A(i)%integer_indices) !For each integer index.

          i_arr = SsTC_integer_memory_element_to_array_element(A(i), i_mem) !Get array index.

          temp_res = cmplx_0
          !_OMPTGT_(PARALLEL PRIVATE (kdotr, vec, prod_R) REDUCTION (+: temp_res))
          !_OMPTGT_(DO)
          do irpts = 1, system%num_R_points !For each point in the Bravais lattice.

            !Compute factor appearing in the exponential (k is in coords relative to recip. lattice vectors).
            kdotr = 2.0_dp*pi*dot_product(system%R_point(irpts, :), k)

            !Compute Bravais lattice vector for label irpts.
            vec = 0.0_dp
            do ivec = 1, 3
              vec = vec + system%R_point(irpts, ivec)*system%direct_lattice_basis(ivec, :)
            enddo

            !In the case of derivatives compute the prefactor involving a product of Bravais lattice vector components.
            prod_R = 1.0_dp
            if (i /= 1) then
              do j = 1, i - 1
                prod_R = prod_R*vec(i_arr(j + 3))
              enddo
            endif

            temp_res = temp_res + & !Compute sum.
                       ((cmplx_i)**(i - 1))*(prod_R)* &
                       exp(cmplx_i*kdotr)*system%real_space_position_elements(i_arr(1), i_arr(2), i_arr(3), irpts) !&
            !/real(system%deg_R_point(irpts), dp)
            !TODO: Recheck the deg_R_points division.
          enddo!irpts
          !_OMPTGT_(END DO)
          !_OMPTGT_(END PARALLEL)
          A(i)%k_data(i_mem) = temp_res
        enddo!i_mem

        deallocate (i_arr)

      enddo
    endif

  end subroutine SsTC_get_position

end module SsTC_local_k_quantities
