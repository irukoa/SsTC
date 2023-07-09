module calculators_floquet

  use utility
  use data_structures
  use local_k_quantities
  use kpath

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
    logical          :: qis0!, qisNaN

    qis0 = (sqrt(sum(q*q)).lt.1.0E-6_dp)
    !qisNaN = (sqrt(sum(q*q)).gt.1.0E6_dp)

    !if (qisNaN) then
    !  H_TK = cmplx(0.0_dp, 0.0_dp)
    !  return !Raise an error?
    !endif

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
          H_nm_TR = modulation*H_nm_TR !TODO: SOME CHECK ON MODULATION AND ISSUE A WARNING?
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

  function quasienergy_kpath_task_constructor(system, Nvec, vec_coord, nkpts, &
                                              Nharm, &
                                              axstart, axend, axsteps, &
                                              pxstart, pxend, pxsteps, &
                                              aystart, ayend, aysteps, &
                                              pystart, pyend, pysteps, &
                                              azstart, azend, azsteps, &
                                              pzstart, pzend, pzsteps, &
                                              omegastart, omegaend, omegasteps, &
                                              t0start, t0end, t0steps) result(default_quasienergy_kpath_task)

    type(sys), intent(in)  :: system
    integer, intent(in) :: Nvec
    real(kind=dp), intent(in) :: vec_coord(Nvec, 3)
    integer, intent(in) :: nkpts(Nvec - 1)

    integer, intent(in) :: Nharm
    real(kind=dp), dimension(Nharm), intent(in) :: axstart, axend, &
                                                   pxstart, pxend, &
                                                   aystart, ayend, &
                                                   pystart, pyend, &
                                                   azstart, azend, &
                                                   pzstart, pzend

    integer, dimension(Nharm), intent(in) :: axsteps, pxsteps, aysteps, pysteps, azsteps, pzsteps

    real(kind=dp), intent(in) :: omegastart, omegaend, &
                                 t0start, t0end

    integer, intent(in) :: omegasteps, t0steps

    real(kind=dp), dimension(6*Nharm + 2) :: start, end
    integer, dimension(6*Nharm + 2) :: steps
    integer :: iharm

    type(k_path_task) :: default_quasienergy_kpath_task

    do iharm = 1, Nharm
      start(6*(iharm - 1) + 1) = axstart(iharm)
      end(6*(iharm - 1) + 1) = axend(iharm)
      steps(6*(iharm - 1) + 1) = axsteps(iharm)

      start(6*(iharm - 1) + 2) = pxstart(iharm)
      end(6*(iharm - 1) + 2) = pxend(iharm)
      steps(6*(iharm - 1) + 2) = pxsteps(iharm)

      start(6*(iharm - 1) + 3) = aystart(iharm)
      end(6*(iharm - 1) + 3) = ayend(iharm)
      steps(6*(iharm - 1) + 3) = aysteps(iharm)

      start(6*(iharm - 1) + 4) = pystart(iharm)
      end(6*(iharm - 1) + 4) = pyend(iharm)
      steps(6*(iharm - 1) + 4) = pysteps(iharm)

      start(6*(iharm - 1) + 5) = azstart(iharm)
      end(6*(iharm - 1) + 5) = azend(iharm)
      steps(6*(iharm - 1) + 5) = azsteps(iharm)

      start(6*(iharm - 1) + 6) = pzstart(iharm)
      end(6*(iharm - 1) + 6) = pzend(iharm)
      steps(6*(iharm - 1) + 6) = pzsteps(iharm)
    enddo

    start(6*Nharm + 1) = t0start
    end(6*Nharm + 1) = t0end
    steps(6*Nharm + 1) = t0steps

    start(6*Nharm + 2) = omegastart
    end(6*Nharm + 2) = omegaend
    steps(6*Nharm + 2) = omegasteps

    default_quasienergy_kpath_task = kpath_constructor(name = "quasienergy", &
                                                 g_calculator = quasienergy, &
                                                 Nvec = Nvec, vec_coord = vec_coord, nkpts = nkpts, &
                                                 N_int_ind = 1, int_ind_range = (/system%num_bands/), &
                                                 N_ext_vars     = 6*Nharm + 2, &
                                                 ext_vars_start = start, &
                                                 ext_vars_end   = end, &
                                                 ext_vars_steps = steps)
  end function quasienergy_kpath_task_constructor

  function quasienergy(floquet_task, system, k, error) result(u)
    class(global_k_data), intent(in) :: floquet_task
    type(sys), intent(in) :: system
    real(kind=dp),          intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp) :: u(product(floquet_task%integer_indices), product(floquet_task%continuous_indices))

    integer :: r_mem, r_arr(size(floquet_task%continuous_indices)), &
               i_mem

    real(kind=dp) :: omega, t0, tper, q(3), dt, &
    quasi(system%num_bands)
    real(kind=dp), allocatable :: amplitudes(:, :), &
                                  phases(:, :)
    integer :: Nharm, iharm, it, i

    complex(kind=dp) :: H_TK(system%num_bands, system%num_bands), &
                        expH_TK(system%num_bands, system%num_bands), &
                        tev(system%num_bands, system%num_bands), &
                        hf(system%num_bands, system%num_bands), &
                        rot(system%num_bands, system%num_bands)

    Nharm = (size(floquet_task%continuous_indices)-2)/6
    allocate(amplitudes(Nharm, 3), phases(Nharm, 3))

    u = cmplx_0

    do r_mem = 1, product(floquet_task%continuous_indices)

      r_arr = continuous_memory_element_to_array_element(floquet_task, r_mem)

      omega = floquet_task%ext_var_data(size(floquet_task%continuous_indices))&
      %data(r_arr(size(floquet_task%continuous_indices)))

      t0 = floquet_task%ext_var_data(size(floquet_task%continuous_indices)-1)&
      %data(r_arr(size(floquet_task%continuous_indices)-1))

      do iharm = 1, Nharm
        amplitudes(iharm, 1) = floquet_task%ext_var_data(6*(iharm - 1) + 1)&
        %data(r_arr(6*(iharm - 1) + 1))
        phases(iharm, 1) = floquet_task%ext_var_data(6*(iharm - 1) + 2)&
        %data(r_arr(6*(iharm - 1) + 2))
        amplitudes(iharm, 2) = floquet_task%ext_var_data(6*(iharm - 1) + 3)&
        %data(r_arr(6*(iharm - 1) + 3))
        phases(iharm, 2) = floquet_task%ext_var_data(6*(iharm - 1) + 4)&
        %data(r_arr(6*(iharm - 1) + 4))
        amplitudes(iharm, 3) = floquet_task%ext_var_data(6*(iharm - 1) + 5)&
        %data(r_arr(6*(iharm - 1) + 5))
        phases(iharm, 3) = floquet_task%ext_var_data(6*(iharm - 1) + 6)&
        %data(r_arr(6*(iharm - 1) + 6))
      enddo

      dt = (2*pi/omega)/real(system%Nt-1, dp)
      tev = cmplx_0
      forall (i = 1:system%num_bands) tev(i, i) = cmplx(1.0_dp, 0.0_dp)
      do it = 1, system%Nt
        tper = t0 + dt*real(it-1, dp) !In eV^-1.
        q = int_driving_field(amplitudes, phases, omega, tper) !In A^-1
        H_TK = wannier_tdep_hamiltonian(system, q, k, system%diag) !In eV.

        expH_TK = utility_exphs(-cmplx_i*dt*H_TK, system%num_bands, .true., error)
        if (error) then
          write(unit=113, fmt="(a, i3, a)") "Error in function quasienergy at t-step, ", it, &
          "when computing matrix exponential for modulation vector q = "
          write(unit=113, fmt="(3E18.8E3, a)") q, "A^-1."
          return
        endif

        tev = matmul(tev, expH_TK)

      enddo
      hf = cmplx_i*omega*utility_logu(tev, system%num_bands, error)/(2*pi)
      if (error) then
        write(unit=113, fmt="(a, i3, a)") "Error in function quasienergy when computing matrix log of the one-petiod time evolution operator."
        return
      endif

      call utility_diagonalize(hf, system%num_bands, quasi, rot, error)
      if (error) then
        write(unit=113, fmt="(a, i3, a)") "Error in function quasienergy when computing the quasienergy spectrum."
        return
      endif
      do i_mem = 1, product(floquet_task%integer_indices)
        u(i_mem, r_mem) = quasi(i_mem)
      enddo
    enddo

    deallocate(amplitudes, phases)

  end function quasienergy

  function int_driving_field(amplitudes, phases, omega, t) result(q)

    real(kind=dp), intent(in) :: amplitudes(:, :), phases(:, :), &
                                 omega, t

    real(kind=dp) :: q(3)

    integer :: icoord, iharm

    q = 0.0_dp

    do iharm = 1, size(amplitudes(:, 1))
      do icoord = 1, 3
        q(icoord) = q(icoord) + &
        amplitudes(iharm, icoord)*sin(iharm*omega*t - phases(iharm, icoord))/&
        (iharm*omega)
      enddo
    enddo

    !q, the integral of the driving electric field, is now in V/(m*eV). 
    !we have to multiply by e to obatin the integral of the driving force field in J/(m*eV)
    !the divide by |e| to obain it in 1/m. Lastly pass to 1/A by multiplying by 1^-10.
    q = q*1.0E-10_dp

  end function int_driving_field

end module calculators_floquet