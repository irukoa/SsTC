module calculators_floquet

  use utility
  use data_structures
  use local_k_quantities

  implicit none

  contains

  function wannier_tdep_hamiltonian(system, q, k, diag) result(H_TK)

    type(sys), intent(in)     :: system

    real(kind=dp), intent(in) :: q(3), & !In A^-1.
                                 k(3)    !In crystal coordineates.

    logical, intent(in)       :: diag !If only WF centres are used.

    complex(kind=dp)          :: H_TK(system%num_bands, system%num_bands)

    integer          :: n, m, l, j, irpts, ivec
    complex(kind=dp) :: modulation, vecpos(3), temp_res, H_nm_TR
    real(kind=dp)    :: vecbrav(3), kdotr
    logical          :: qis0

    qis0 = (sqrt(sum(q*q)).lt.1.0E-6_dp)

    do n = 1, system%num_bands
      do m = 1, system%num_bands

        temp_res = cmplx_0
        !$OMP PARALLEL PRIVATE (H_nm_TR, modulation, vecpos, kdotr, vecbrav) REDUCTION (+: temp_res)
        !$OMP DO
        do irpts = 1, system%num_R_points

          H_nm_TR = system%real_space_hamiltonian_elements(n, m, irpts)
          modulation = cmplx(1.0_dp, 0.0_dp)
          if (qis0) goto 1 !Ignore modulation, we never needed the rotating frame transformation.


          if (diag) then !If only diagonal elements of the pos. operator exist we should not count for bands l, j.

            vecpos =  system%real_space_position_elements(m, m, :, irpts) - system%real_space_position_elements(n, n, :, irpts)
            modulation = exp(cmplx_i*sum(q*vecpos))

          else

            do l = 1, system%num_bands
              do j = 1, system%num_bands

                vecpos =  system%real_space_position_elements(m, l, :, irpts) - system%real_space_position_elements(j, n, :, irpts)

                modulation = modulation + &
                exp(cmplx_i*sum(q*vecpos))

              enddo!j
            enddo!l

          endif

          1 continue
          H_nm_TR = modulation*H_nm_TR
          !At this point H_nm_TR stores the real lattice (Wannier) resolved time
          !dependent Hamiltonian for the force modulation q and R-point irpts.

          !Now we compute its Fourier transform.

          !Compute factor appearing in the exponential (k is in coords relative to recip. lattice vectors).
          kdotr = 2.0_dp*pi*dot_product(system%R_point(irpts, :), k)

          !Compute Bravais lattice vector for label irpts.
          vecbrav = 0.0_dp
          do ivec = 1, 3
            vecbrav = vecbrav + system%R_point(irpts, ivec)*system%direct_lattice_basis(ivec, :)
          enddo

          temp_res = temp_res + & !Compute sum.
          exp(cmplx_i*kdotr)*H_nm_TR &
          /real(system%deg_R_point(irpts), dp)

        enddo!irpts
        !$OMP END DO
        !$OMP END PARALLEL
        H_TK(n, m) = temp_res
      enddo!m
    enddo!n

  end function wannier_tdep_hamiltonian

!  !function tdep_hamiltonian(system, omega, t, amplitudes)
!
!  function driving_field(amplitudes, phases, omega, t, t0) result(u)
!    real(kind=dp), intent(in) :: amplitudes(:, :), phases(:, :), &
!                                 omega, t, t0
!
!    complex(kind=dp) :: u(3)
!
!    integer :: icoord, iharm
!    real(kind=dp) :: n
!
!    u = cmplx_0
!
!    do iharm = 1, size(amplitudes(:, 1))
!      n = real(iharm, dp) - 0.5_dp*real(size(amplitudes(:, 1)) + 1, dp)
!      do icoord = 1, 3
!        u(icoord) = u(icoord) + &
!        amplitudes(iharm, icoord)*exp(cmplx_i*phases(iharm, icoord))*&
!        exp(cmplx_i*n*omega*(t-t0))
!      enddo
!    enddo
!
!  end function driving_field
!
!  function int_driving_field(amplitudes, phases, omega, t, t0) result(u)
!    real(kind=dp), intent(in) :: amplitudes(:, :), phases(:, :), &
!                                 omega, t, t0
!
!    complex(kind=dp) :: u(3)
!
!    integer :: icoord, iharm
!    real(kind=dp) :: n
!
!    u = cmplx_0
!
!    do iharm = 1, size(amplitudes(:, 1))
!      n = real(iharm, dp) - 0.5_dp*real(size(amplitudes(:, 1)) + 1, dp)
!      if (n == 0) cycle
!      do icoord = 1, 3
!        u(icoord) = u(icoord) + &
!        amplitudes(iharm, icoord)*exp(cmplx_i*phases(iharm, icoord))*&
!        exp(cmplx_i*n*omega*(t-t0))/&
!        (cmplx_i*n*omega)
!      enddo
!    enddo
!
!  end function int_driving_field
!
!  function floq_current(system, floquet_curr, k) result(j)
!    class(global_k_data), intent(in) :: floquet_curr
!    type(sys),           intent(in) :: system
!    real(kind=dp),       intent(in) :: k(3)
!
!    complex(kind=dp)                :: j(product(floquet_curr%integer_indices), product(floquet_curr%continuous_indices))
!
!    integer :: ie, iw, it
!
!    real(kind=dp) :: e, omega, t
!
!    count = 1
!
!    do iharm = 1, (size(floquet_curr%continuous_indices)-2)/6 !For each harmonic index.
!      do icomp = 1, 3 !For each cartesian index.
!
!        !do iamp = 1, size(floq_curr%continuous_indices(count)%data) !For each considered change in amplitude with id iharm, icomp.
!          !do iphs = 1, size(floq_curr%continuous_indices(count + 1)%data) !For each considered change in phase with id iharm, icomp.
!            amplitudes(iharm, icomp)%data = floq_curr%continuous_indices(count)%data!(iamp)
!            phases(iharm, icomp)%data = floq_curr%continuous_indices(count + 1)%data!(iphs)
!          !enddo
!        !enddo
!        count = count + 2
!
!      enddo
!    enddo
!
!    do iamp = 1, size(amplitudes(iharm, icomp)%data)
!
!
!
!
!
!
!   BASICAMENTE ESTO ES LO QUE MAS CONVIENE.
!    PARA CADA CASO QUE OCNSIDEREMOS CONVIENE CREAR UN CAMPO ELECTRICO 
!    QUE CORRESPONDA A UNA INTERFACE GENERAL. LO QUE SE ME OCURRE ES
!    AMPLITUDA VARIABLE (HASTA 2 AMPLITUDES), FRECUENCIA VARIABLE, POLARIZACION VARIABLE
!    POLARIZACION Y FRECUENCIA VARIABLES.  
 !    !do ie = 1, floquet_curr%continuous_indices(1)
!    !  e = floquet_curr%ext_var_data(1)%data(ie)
!
!    !  do iw = 1, floquet_curr%continuous_indices(2)
!    !    omega = floquet_curr%ext_var_data(2)%data(iw)
!
!    !    do it = 1, floquet_curr%continuous_indices(3)
!    !      t = floquet_curr%ext_var_data(3)%data(it)
!
!    !      rarr = (/ie, iw, it/)
!    !      rmem = ...r_arr
!    !      j(, r_mem)
!    !    enddo
!
!    !  enddo
!
!    !enddo
!
!  end function floq_current
!
end module calculators_floquet