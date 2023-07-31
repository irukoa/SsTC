module SsTC_calculators_floquet

  use SsTC_utility
  use SsTC_extrapolation_integration
  use SsTC_data_structures
  use SsTC_local_k_quantities
  use SsTC_kpath
  use SsTC_integrator

  implicit none

  private

  type, extends(SsTC_BZ_integral_task) :: SsTC_floq_BZ_integral_task
    integer              :: Nt = 65 !2^6 + 1 Discretization points of each period.
    integer            :: Ns = 10 !Considered Harmonics.
    logical                       :: diag = .false. !If we only consider diagonal terms of the pos. operator.
  end type SsTC_floq_BZ_integral_task

  type, extends(SsTC_kpath_task) :: SsTC_floq_kpath_task
    integer              :: Nt = 65 !2^6 + 1 Discretization points of each period.
    integer              :: Ns = 10 !Considered Harmonics.
    logical                       :: diag = .false. !If we only consider diagonal terms of the pos. operator.
  end type SsTC_floq_kpath_task

  public :: SsTC_floq_BZ_integral_task
  public :: SsTC_floq_kpath_task

  public :: SsTC_quasienergy_kpath_task_constructor
  public :: SsTC_quasienergy

  public :: SsTC_floq_curr_BZ_integral_constructor
  public :: SsTC_floq_curr

contains

  !====DEFAULT CALCULATORS====!

  !====QUASIENERGY CALCULATOR AND CONSTRUCTOR====!
  subroutine SsTC_quasienergy_kpath_task_constructor(floq_task, system, Nvec, vec_coord, nkpts, &
                                                     Nharm, &
                                                     axstart, axend, axsteps, &
                                                     pxstart, pxend, pxsteps, &
                                                     aystart, ayend, aysteps, &
                                                     pystart, pyend, pysteps, &
                                                     azstart, azend, azsteps, &
                                                     pzstart, pzend, pzsteps, &
                                                     omegastart, omegaend, omegasteps, &
                                                     t0start, t0end, t0steps, &
                                                     floq_Nt, floq_NS, floq_diag)

    type(SsTC_sys), intent(in)  :: system
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

    integer, optional, intent(in) :: floq_Nt, floq_NS
    logical, optional, intent(in) :: floq_diag

    type(SsTC_floq_kpath_task), intent(out) :: floq_task

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

    call SsTC_kpath_constructor(task=floq_task, name="quasienergy", &
                                g_calculator=SsTC_quasienergy, &
                                Nvec=Nvec, vec_coord=vec_coord, nkpts=nkpts, &
                                N_int_ind=1, int_ind_range=(/system%num_bands/), &
                                N_ext_vars=6*Nharm + 2, &
                                ext_vars_start=start, &
                                ext_vars_end=end, &
                                ext_vars_steps=steps)

    if (present(floq_Nt)) floq_task%Nt = floq_Nt
    if (present(floq_Ns)) floq_task%Ns = floq_Ns
    if (present(floq_diag)) floq_task%diag = floq_diag

  end subroutine SsTC_quasienergy_kpath_task_constructor
  !==========DEFAULT QUASIENERGY KPATH TASK==========!
  function SsTC_quasienergy(floquet_task, system, k, error) result(u)
    class(SsTC_global_k_data), intent(in) :: floquet_task
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
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

    select type (floquet_task)
    type is (SsTC_floq_kpath_task)

      Nharm = (size(floquet_task%continuous_indices) - 2)/6
      allocate (amplitudes(Nharm, 3), phases(Nharm, 3))

      u = cmplx_0

      do r_mem = 1, product(floquet_task%continuous_indices)

        r_arr = SsTC_continuous_memory_element_to_array_element(floquet_task, r_mem)

        omega = floquet_task%ext_var_data(size(floquet_task%continuous_indices)) &
                %data(r_arr(size(floquet_task%continuous_indices)))

        t0 = floquet_task%ext_var_data(size(floquet_task%continuous_indices) - 1) &
             %data(r_arr(size(floquet_task%continuous_indices) - 1))

        do iharm = 1, Nharm
          amplitudes(iharm, 1) = floquet_task%ext_var_data(6*(iharm - 1) + 1) &
                                 %data(r_arr(6*(iharm - 1) + 1))
          phases(iharm, 1) = floquet_task%ext_var_data(6*(iharm - 1) + 2) &
                             %data(r_arr(6*(iharm - 1) + 2))
          amplitudes(iharm, 2) = floquet_task%ext_var_data(6*(iharm - 1) + 3) &
                                 %data(r_arr(6*(iharm - 1) + 3))
          phases(iharm, 2) = floquet_task%ext_var_data(6*(iharm - 1) + 4) &
                             %data(r_arr(6*(iharm - 1) + 4))
          amplitudes(iharm, 3) = floquet_task%ext_var_data(6*(iharm - 1) + 5) &
                                 %data(r_arr(6*(iharm - 1) + 5))
          phases(iharm, 3) = floquet_task%ext_var_data(6*(iharm - 1) + 6) &
                             %data(r_arr(6*(iharm - 1) + 6))
        enddo

        dt = (2*pi/omega)/real(floquet_task%Nt - 1, dp)
        tev = cmplx_0
        forall (i=1:system%num_bands) tev(i, i) = cmplx(1.0_dp, 0.0_dp)

        do it = 2, floquet_task%Nt

          tper = t0 + dt*real(it - 1, dp) !In eV^-1.

          q = SsTC_int_driving_field(amplitudes, phases, omega, tper) !In A^-1

          H_TK = SsTC_wannier_tdep_hamiltonian(system, q, k, floquet_task%diag, error) !In eV.
          if (error) then
            write (unit=stderr, fmt="(a, i4, a)") "Error in function quasienergy at t-step, ", it, &
              "when computing time-dependent Hamiltonian for modulation vector q = "
            write (unit=stderr, fmt="(3e18.8e3, a)") q, "A^-1."
            return
          endif

          expH_TK = SsTC_utility_exphs(-cmplx_i*dt*H_TK, system%num_bands, .true., error)
          if (error) then
            write (unit=stderr, fmt="(a, i4, a)") "Error in function quasienergy at t-step, ", it, &
              "when computing matrix exponential for modulation vector q = "
            write (unit=stderr, fmt="(3e18.8e3, a)") q, "A^-1."
            return
          endif

          tev = matmul(tev, expH_TK)

        enddo

        hf = cmplx_i*omega*SsTC_utility_logu(tev, system%num_bands, error)/(2*pi)
        if (error) then
          write (unit=stderr, fmt="(a)") "Error in function quasienergy when computing &
          & matrix log of the one-petiod time evolution operator."
          return
        endif

        call SsTC_utility_diagonalize(hf, system%num_bands, quasi, rot, error)
        if (error) then
          write (unit=stderr, fmt="(a)") "Error in function quasienergy when computing the quasienergy spectrum."
          return
        endif

        do i_mem = 1, product(floquet_task%integer_indices)
          u(i_mem, r_mem) = quasi(i_mem)
        enddo

      enddo

      deallocate (amplitudes, phases)

    end select

  end function SsTC_quasienergy
  !====END QUASIENERGY CALCULATOR AND CONSTRUCTOR====!

  !====FLOQUET CURRENT CALCULATOR AND CONSTRUCTOR====!
  subroutine SsTC_floq_curr_BZ_integral_constructor(floq_task, method, samples, &
                                                    Nharm, &
                                                    axstart, axend, axsteps, &
                                                    pxstart, pxend, pxsteps, &
                                                    aystart, ayend, aysteps, &
                                                    pystart, pyend, pysteps, &
                                                    azstart, azend, azsteps, &
                                                    pzstart, pzend, pzsteps, &
                                                    omegastart, omegaend, omegasteps, &
                                                    t0start, t0end, t0steps, &
                                                    tstart, tend, tsteps, &
                                                    floq_Nt, floq_NS, floq_diag, &
                                                    particular_integer_component)

    character(len=*), optional, intent(in) :: method
    integer, optional, intent(in) :: samples(3)

    integer, intent(in) :: Nharm
    real(kind=dp), dimension(Nharm), intent(in) :: axstart, axend, &
                                                   pxstart, pxend, &
                                                   aystart, ayend, &
                                                   pystart, pyend, &
                                                   azstart, azend, &
                                                   pzstart, pzend

    integer, dimension(Nharm), intent(in) :: axsteps, pxsteps, aysteps, pysteps, azsteps, pzsteps

    real(kind=dp), intent(in) :: omegastart, omegaend, &
                                 t0start, t0end, &
                                 tstart, tend

    integer, intent(in) :: omegasteps, t0steps, tsteps

    real(kind=dp), dimension(6*Nharm + 3) :: start, end
    integer, dimension(6*Nharm + 3) :: steps
    integer :: i, iharm, iterable_vars(6*Nharm)

    integer, optional, intent(in) :: floq_Nt, floq_NS
    logical, optional, intent(in) :: floq_diag

    integer, optional, intent(in) :: particular_integer_component(:)

    type(SsTC_floq_BZ_integral_task), intent(out) :: floq_task

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

    start(6*Nharm + 3) = tstart
    end(6*Nharm + 3) = tend
    steps(6*Nharm + 3) = tsteps

    call SsTC_BZ_integral_task_constructor(task=floq_task, name="floq_curr", &
                                           g_calculator=SsTC_floq_curr, &
                                           method=method, samples=samples, &
                                           N_int_ind=1, int_ind_range=(/3/), &
                                           N_ext_vars=6*Nharm + 3, &
                                           ext_vars_start=start, &
                                           ext_vars_end=end, &
                                           ext_vars_steps=steps, &
                                           part_int_comp=particular_integer_component)

    !Function floq_curr assumes that all the information on the driving field is stored in the iterable.
    !Information on omega, t0 and t is passed via the last 3 indices of start(:), end(:) and steps(:) arrays.
    forall (i=1:6*Nharm) iterable_vars(i) = i
    call SsTC_construct_iterable(floq_task, iterable_vars)

    if (present(floq_Nt)) floq_task%Nt = floq_Nt
    if (present(floq_Ns)) floq_task%Ns = floq_Ns
    if (present(floq_diag)) floq_task%diag = floq_diag

  end subroutine SsTC_floq_curr_BZ_integral_constructor
  !==========DEFAULT FLOQUET CURRENT INTEGRAL TASK==========!
  function SsTC_floq_curr(floquet_task, system, k, error) result(u)
    class(SsTC_global_k_data), intent(in) :: floquet_task
    type(SsTC_sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp) :: u(product(floquet_task%integer_indices), product(floquet_task%continuous_indices))

    integer :: r_mem, r_arr(size(floquet_task%continuous_indices)), &
               i_mem

    real(kind=dp) :: omega, t0, q(3), dt, tper, t, quasi(system%num_bands), eig(system%num_bands)
    real(kind=dp), allocatable :: amplitudes(:, :), &
                                  phases(:, :)
    integer :: Nharm, iharm, i, iterable, it0, iomega, it, itper, info, &
               n, m, l, p, ir, is

    complex(kind=dp) :: H_TK(system%num_bands, system%num_bands), &
                        expH_TK(system%num_bands, system%num_bands), &
                        hf(system%num_bands, system%num_bands), &
                        rotWF(system%num_bands, system%num_bands), &
                        w_hamiltonian(system%num_bands, system%num_bands), &
                        rotWH(system%num_bands, system%num_bands), &
                        rho(system%num_bands, system%num_bands), &
                        vels(system%num_bands, system%num_bands, 3), &
                        rhoF(system%num_bands, system%num_bands), &
                        velsF(system%num_bands, system%num_bands, 3)

    complex(kind=dp), allocatable :: tev(:, :, :), pt(:, :, :), shrinkqs(:), qs(:, :, :), integrand(:, :, :), shrink_integrand(:, :)

    select type (floquet_task)
    type is (SsTC_floq_BZ_integral_task)

      !Gather required quantities in the Wannier basis.
      !Time-independent Hamiltonian.
      w_hamiltonian = SsTC_wannier_hamiltonian(system, k)
      !Get eigenvalues and rotation from Wannier to Hamiltonian basis.
      call SsTC_utility_diagonalize(w_hamiltonian, system%num_bands, eig, rotWH, error)
      if (error) then
        write (unit=stderr, fmt="(a)") "Error in function floq_curr when computing&
        & the eigenvalues of the time-independent Hamiltonian."
        return
      endif
      !Get occupations (a matrix in the Hamiltonian basis).
      rho = SsTC_hamiltonian_occ_matrix(system, eig)
      !Rotate them back to the Wannier basis (notice W\rho W^\dagger instead of W^\dagger\rho W).
      rho = matmul(matmul(rotWH, rho), transpose(conjg(rotWH)))
      !Get velocities (Hamiltonian basis). TODO: CHECK IF WE ALSO NEED NONDIAGONAL ELEMENTS OR ONLY DIAGONAL ONES.
      vels = SsTC_velocities(system, SsTC_wannier_dhamiltonian_dk(system, k), eig, rotWH, error)
      if (error) then
        write (unit=stderr, fmt="(a)") "Error in function floq_curr when computing the velocities in the degenerate subspace."
        return
      endif
      !Rotate back to Wannier basis.
      do i = 1, 3
        vels(:, :, i) = matmul(matmul(rotWH, vels(:, :, i)), transpose(conjg(rotWH)))
      enddo

      Nharm = (size(floquet_task%continuous_indices) - 3)/6
      allocate (amplitudes(Nharm, 3), phases(Nharm, 3))

      u = cmplx_0

      allocate (tev(floquet_task%Nt, system%num_bands, system%num_bands), &
                pt(floquet_task%Nt, system%num_bands, system%num_bands), &
                qs(-floquet_task%Ns:floquet_task%Ns, system%num_bands, system%num_bands), &
                shrinkqs(system%num_bands*system%num_bands), &
                integrand(floquet_task%Nt, system%num_bands, system%num_bands), &
                shrink_integrand(floquet_task%Nt, system%num_bands*system%num_bands))

      do iterable = 1, size(floquet_task%iterables(:, 1)) !Loop on: all ext. variables except t0, omega, t

        r_arr = floquet_task%iterables(iterable, :) !This list contains the particular permutation involving variation of driving field params.
        !The indices corresponding to t0, omega and t (the last 3) are all 1.

        do iharm = 1, Nharm
          amplitudes(iharm, 1) = floquet_task%ext_var_data(6*(iharm - 1) + 1) &
                                 %data(r_arr(6*(iharm - 1) + 1))
          phases(iharm, 1) = floquet_task%ext_var_data(6*(iharm - 1) + 2) &
                             %data(r_arr(6*(iharm - 1) + 2))
          amplitudes(iharm, 2) = floquet_task%ext_var_data(6*(iharm - 1) + 3) &
                                 %data(r_arr(6*(iharm - 1) + 3))
          phases(iharm, 2) = floquet_task%ext_var_data(6*(iharm - 1) + 4) &
                             %data(r_arr(6*(iharm - 1) + 4))
          amplitudes(iharm, 3) = floquet_task%ext_var_data(6*(iharm - 1) + 5) &
                                 %data(r_arr(6*(iharm - 1) + 5))
          phases(iharm, 3) = floquet_task%ext_var_data(6*(iharm - 1) + 6) &
                             %data(r_arr(6*(iharm - 1) + 6))
        enddo

        do it0 = 1, floquet_task%continuous_indices(6*Nharm + 1) !t0

          r_arr(6*Nharm + 1) = it0

          t0 = floquet_task%ext_var_data(6*Nharm + 1) &
               %data(it0)

          do iomega = 1, floquet_task%continuous_indices(6*Nharm + 2) !w

            r_arr(6*Nharm + 2) = iomega

            omega = floquet_task%ext_var_data(6*Nharm + 2) &
                    %data(iomega)

            dt = (2*pi/omega)/real(floquet_task%Nt - 1, dp)

            !This calculates the time-evolution operator in the Wannier basis for each time-instant within a period.

            !Initialize array and 1st time instant.
            tev = cmplx_0
            forall (i=1:system%num_bands) tev(1, i, i) = cmplx(1.0_dp, 0.0_dp)

            do itper = 2, floquet_task%Nt !For all remaining time instants.

              tper = t0 + dt*real(itper - 1, dp) !In eV^-1.

              q = SsTC_int_driving_field(amplitudes, phases, omega, tper) !In A^-1

              H_TK = SsTC_wannier_tdep_hamiltonian(system, q, k, floquet_task%diag, error) !In eV.
              if (error) then
                write (unit=stderr, fmt="(a, i4, a)") "Error in function floq_curr at t-step, ", it, &
                  "when computing time-dependent Hamiltonian for modulation vector q = "
                write (unit=stderr, fmt="(3e18.8e3, a)") q, "A^-1."
                return
              endif

              expH_TK = SsTC_utility_exphs(-cmplx_i*dt*H_TK, system%num_bands, .true., error)
              if (error) then
                write (unit=stderr, fmt="(a, i4, a)") "Error in function floq_curr at t-step, ", it, &
                  "when computing matrix exponential for modulation vector q = "
                write (unit=stderr, fmt="(3e18.8e3, a)") q, "A^-1."
                return
              endif

              tev(itper, :, :) = matmul(tev(itper - 1, :, :), expH_TK)

            enddo !itper

            !At this point, the time-evolution operator for each time-instant within a period has been calculated in the Wannier basis.

            !Get effective Floquet Hamiltonian in the Wannier basis.
            hf = cmplx_i*omega*SsTC_utility_logu(tev(floquet_task%Nt, :, :), system%num_bands, error)/(2*pi)
            if (error) then
              write (unit=stderr, fmt="(a)") "Error in function floq_curr when computing &
              & matrix log of the one-petiod time evolution operator."
              return
            endif

            !Diagonalize it to obtain the quasienergy spectra and the rotation matrix passing from Wannier to Floquet basis.
            call SsTC_utility_diagonalize(hf, system%num_bands, quasi, rotWF, error)
            if (error) then
              write (unit=stderr, fmt="(a)") "Error in function floq_curr when computing the quasienergy spectrum."
              return
            endif

            !Get Q_s operators in the Floquet basis.
            !For each omega multiple,
            do is = -floquet_task%Ns, floquet_task%Ns
              !for each time instant within the period,
              do itper = 1, floquet_task%Nt

                tper = t0 + dt*real(itper - 1, dp) !In eV^-1.

                !compute the contribution for each time-instant.
                integrand(itper, :, :) = matmul(TEV(itper, :, :), &
                                                SsTC_utility_exphs(cmplx_i*hf*tper, system%num_bands, .true., error)) &
                                         *exp(-cmplx_i*is*omega*tper) &
                                         /real(floquet_task%Nt - 1, dp)
                !Adimensional. Units are OK.

                !Shrink band indices.
                call SsTC_shrink_array(integrand(itper, :, :), shrink_integrand(itper, :), info) !TODO: if info=...???

              enddo

              !Extrapolate.
              call SsTC_integral_extrapolation(shrink_integrand, (/floquet_task%Nt/), (/t0, t0 + 2*pi/omega/), shrinkqs, info) !TODO: if info=...???
              !Expand to Q_s.
              call SsTC_expand_array(shrinkqs, qs(is, :, :), info) !TODO: if info=...???
              !At this point, the Q_s are in the Wannier basis, we pass them to the Floquet basis.
              qs(is, :, :) = matmul(matmul(transpose(conjg(rotWF)), qs(is, :, :)), rotWF)
            enddo

            !Right now, we have Q_s in the Floquet basis, the quasienergy spectra and the rotation matrix from Wannier to Floquet basis.
            !We also have computed the occupation matrix and velocity matrix in the Wannier basis (rho and vels respectively).
            !We create a local copy of both objects in the Floquet basis,
            rhoF = matmul(matmul(transpose(conjg(rotWF)), rho), rotWF)
            do i = 1, 3
              velsF(:, :, i) = matmul(matmul(transpose(conjg(rotWF)), vels(:, :, i)), rotWF)
            enddo
            !so now everything is resolved in the Floquet basis.

            !For each "external" time instant.
            do it = 1, floquet_task%continuous_indices(6*Nharm + 3) !t

              r_arr(6*Nharm + 3) = it

              t = floquet_task%ext_var_data(6*Nharm + 3) &
                  %data(it)!In s.
              !We pass to eV^-1 by dividing by hbar_over_e
              t = t/hbar_over_e

              r_mem = SsTC_continuous_array_element_to_memory_element(floquet_task, r_arr)

              do i_mem = 1, product(floquet_task%integer_indices)

                do n = 1, system%num_bands
                do m = 1, system%num_bands
                do l = 1, system%num_bands
                do p = 1, system%num_bands
                do ir = -floquet_task%Ns, floquet_task%Ns
                do is = -floquet_task%Ns, floquet_task%Ns

                  u(i_mem, r_mem) = u(i_mem, r_mem) + &
                                    rhoF(n, m)*velsF(l, p, i_mem)*qs(ir, p, n)*conjg(qs(is, l, m))* &
                                    exp(-cmplx_i*t*(quasi(n) - quasi(m) + real(ir - is, dp)*omega))!Units: eV*Angstrom.

                enddo
                enddo
                enddo
                enddo
                enddo
                enddo

              enddo !i_mem

            enddo !t
          enddo !w
        enddo !t0

      enddo!iterable

      deallocate (amplitudes, phases, tev, pt, qs, shrink_integrand, shrinkqs)
    end select

    u = u*1.0_dp !Factor in PostW90: fac = (1.0E20_dp*physics%elem_charge_SI**2)/(physics%hbar_SI*cell_volume)

    !stop !TODO: Testing purposes.

  end function SsTC_floq_curr
  !====END FLOQUET CURRENT CALCULATOR AND CONSTRUCTOR====!

  !====CORE PROCEDURES====!
  function SsTC_wannier_tdep_hamiltonian(system, q, k, diag, error) result(H_TK)

    type(SsTC_sys), intent(in)     :: system

    real(kind=dp), intent(in) :: q(3), & !In A^-1.
                                 k(3)    !In crystal coordineates.

    logical, intent(in)       :: diag !If only WF centres are used.
    logical, intent(inout) ::  error

    complex(kind=dp)          :: H_TK(system%num_bands, system%num_bands)

    integer          :: n, m, l, j, irpts, ivec
    complex(kind=dp) :: modulation, vecpos(3), temp_res, H_nm_TR
    real(kind=dp)    :: vecbrav(3), kdotr
    logical          :: qis0

    qis0 = (sqrt(sum(q*q)) .lt. 1.0E-6_dp)

    do n = 1, system%num_bands
      do m = 1, system%num_bands

        temp_res = cmplx_0
!$OMP         PARALLEL PRIVATE (H_nm_TR, modulation, vecpos, kdotr, vecbrav) REDUCTION (+: temp_res)
!$OMP         DO
        do irpts = 1, system%num_R_points

          H_nm_TR = system%real_space_hamiltonian_elements(n, m, irpts)
          modulation = cmplx(1.0_dp, 0.0_dp)
          if (qis0) goto 1 !Ignore modulation, we never needed the rotating frame transformation.

          if (diag) then !If only diagonal elements of the pos. operator exist we should not count for bands l, j.

            vecpos = system%real_space_position_elements(m, m, :, irpts) - system%real_space_position_elements(n, n, :, irpts)
            modulation = exp(cmplx_i*sum(q*vecpos))

          else

            do l = 1, system%num_bands
              do j = 1, system%num_bands

                vecpos = system%real_space_position_elements(m, l, :, irpts) - system%real_space_position_elements(j, n, :, irpts)

                modulation = modulation + &
                             exp(cmplx_i*sum(q*vecpos))

              enddo!j
            enddo!l

          endif

1         continue

          if (modulation /= modulation) then
            !In this case, modulation is NaN.
            error = .true.
          endif
          if (error) cycle

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
!$OMP         END DO
!$OMP         END PARALLEL
        H_TK(n, m) = temp_res
      enddo!m
    enddo!n
    if (error) then
      write (unit=stderr, fmt="(a)") "Error in function wannier_tdep_hamiltonian: &
      & a large q has made the modulation factor become NaN and as a consequence H(t) is undetermined."
    endif

  end function SsTC_wannier_tdep_hamiltonian

  function SsTC_int_driving_field(amplitudes, phases, omega, t) result(q)

    real(kind=dp), intent(in) :: amplitudes(:, :), phases(:, :), &
                                 omega, t

    real(kind=dp) :: q(3)

    integer :: icoord, iharm

    q = 0.0_dp

    do iharm = 1, size(amplitudes(:, 1))
      do icoord = 1, 3
        q(icoord) = q(icoord) + &
                    amplitudes(iharm, icoord)*sin(iharm*omega*t - phases(iharm, icoord))/ &
                    (iharm*omega)
      enddo
    enddo

    !q, the integral of the driving electric field, is now in V/(m*eV).
    !we have to multiply by e to obatin the integral of the driving force field in J/(m*eV)
    !the divide by |e| to obain it in 1/m. Lastly pass to 1/A by multiplying by 1^-10.
    q = q*1.0E-10_dp

  end function SsTC_int_driving_field
  !=================================!

end module SsTC_calculators_floquet
