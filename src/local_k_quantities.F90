module local_k_quantities

  use utility
  use data_structures

  implicit none

  contains

  subroutine get_hamiltonian(system, k, H, Nder_i, only_i)
    type(sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    integer, optional :: Nder_i !Order of the derivative. Nder = 0 means normal Hamiltonian, Nder = 1, \partial H/\partial k^i...
    logical, optional :: only_i !

    type(local_k_data), allocatable, intent(out) :: H(:)

    integer :: Nder
    logical :: only

    integer :: i, j, i_mem, irpts, ivec
    integer, allocatable ::  i_arr(:)

    real(kind=dp) :: kdotr, prod_R, vec(3)

    if (.not.(present(Nder_i))) then
      Nder = 0
    else
      Nder = Nder_i
    endif

    if (.not.(present(only_i))) then
      only = .true.
    else
      only = only_i
    endif

    if ((only .EQV. .true.)) allocate(H(1)) !Only compute the Nth Hamiltonian derivative.
    if ((only .EQV. .false.)) allocate(H(Nder + 1)) !Compute from the 0th to the Nth Hamiltonian derivative.

    if ((only .EQV. .false.)) then
      do i = 1, Nder + 1
        allocate(H(i)%integer_indices(2 + i - 1)) !2 band indices and i - 1derivative indices.
        H(i)%integer_indices(1) = system%num_bands
        H(i)%integer_indices(2) = system%num_bands
        if (i/=1) then
          do j = 1, i - 1
            H(i)%integer_indices(2 + j) = 3
          enddo
        endif
      enddo
    else
      allocate(H(1)%integer_indices(2 + Nder)) !2 band indices and Nder derivative indices.
      H(1)%integer_indices(1) = system%num_bands
      H(1)%integer_indices(2) = system%num_bands
      if (Nder/=0) then
        do i = 1, Nder
          H(1)%integer_indices(2 + i) = 3
        enddo
      endif
    endif

    if ((only .EQV. .true.)) then !Case of only a requested derivative.

      allocate(i_arr(size(H(1)%integer_indices))) !Array indices storage.
      allocate(H(1)%k_data(product(H(1)%integer_indices))) !Data storage.
      H(1)%k_data = cmplx_0 !Initialize.

      do i_mem = 1, product(H(1)%integer_indices) !For each integer index.

        i_arr = integer_memory_element_to_array_element(H(1), i_mem) !Get array index.

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
          if (Nder/=0) then
            do j = 1, Nder
              prod_R = prod_R*vec(i_arr(j + 2))
            enddo
          endif

          H(1)%k_data(i_mem) = H(1)%k_data(i_mem) + & !Compute sum.
          ((cmplx_i)**Nder)*(prod_R)*exp(cmplx_i*kdotr)*system%real_space_hamiltonian_elements(i_arr(1), i_arr(2), irpts)

        enddo!irpts
      enddo!i_mem

      deallocate(i_arr)

    else !Case of requestad all derivatives.
      do i = 1, Nder + 1



        allocate(i_arr(size(H(i)%integer_indices))) !Array indices storage.
        allocate(H(i)%k_data(product(H(i)%integer_indices))) !Data storage.
        H(i)%k_data = cmplx_0 !Initialize.
  
        do i_mem = 1, product(H(i)%integer_indices) !For each integer index.
  
          i_arr = integer_memory_element_to_array_element(H(i), i_mem) !Get array index.
  
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
            if (i/=1) then
              do j = 1, i - 1
                prod_R = prod_R*vec(i_arr(j + 2))
              enddo
            endif
  
            H(i)%k_data(i_mem) = H(i)%k_data(i_mem) + & !Compute sum.
            ((cmplx_i)**(i - 1))*(prod_R)*exp(cmplx_i*kdotr)*system%real_space_hamiltonian_elements(i_arr(1), i_arr(2), irpts)
  
          enddo!irpts
        enddo!i_mem
  
        deallocate(i_arr)




      enddo
    endif

  end subroutine get_hamiltonian

end module local_k_quantities