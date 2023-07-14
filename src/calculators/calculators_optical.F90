module calculators_optical

  use utility
  use data_structures
  use local_k_quantities
  use integrator

  implicit none

  type, extends(BZ_integral_task) :: optical_BZ_integral_task
    logical                       :: adpt_smearing = .true.
    real(kind=dp)                 :: smearing = 1.0_dp
  end type

contains

  subroutine default_optical_conductivity_constructor(optical_task, method, samples, &
                                                      omegastart, omegaend, omegasteps, &
                                                      particular_integer_component, &
                                                      optical_smearing)

    character(len=*), optional, intent(in) :: method
    integer, optional, intent(in) :: samples(3)

    real(kind=dp), intent(in) :: omegastart, omegaend
    integer, intent(in) :: omegasteps

    integer, optional, intent(in) :: particular_integer_component(:)
    real(kind=dp), optional, intent(in) ::  optical_smearing

    type(optical_BZ_integral_task), intent(out) :: optical_task

    call BZ_integral_task_constructor(task=optical_task, name="opt_cond", &
                                      g_calculator=optical_conductivity, &
                                      method=method, samples=samples, &
                                      N_int_ind=2, int_ind_range=(/3, 3/), &
                                      N_ext_vars=1, &
                                      ext_vars_start=(/omegastart/), &
                                      ext_vars_end=(/omegaend/), &
                                      ext_vars_steps=(/omegasteps/), &
                                      part_int_comp=particular_integer_component)

    if (present(optical_smearing)) then
      optical_task%adpt_smearing = .false.
      optical_task%smearing = optical_smearing
    endif

  end subroutine default_optical_conductivity_constructor

  function optical_conductivity(task, system, k, error) result(u)
    class(global_k_data), intent(in) :: task
    type(sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp)                   :: u(product(task%integer_indices), product(task%continuous_indices))

    integer :: i_mem, i_arr(product(task%integer_indices)), &
               r_mem
    integer :: i, j, n, m

    real(kind=dp) :: omega, eig(system%num_bands), smearing, arg, delta, bpart, spacing, dk
    complex(kind=dp), dimension(system%num_bands, system%num_bands) :: &
      w_hamiltonian, &
      rho, &
      rot
    complex(kind=dp), dimension(system%num_bands, system%num_bands, 3) :: &
      w_dk_hamiltonian, &
      w_connection, connection, &
      na_d, vels

    select type (task)
    type is (optical_BZ_integral_task)

      u = cmplx(0.0_dp, 0.0_dp)

      !Gather required quantities in the Wannier basis.
      w_hamiltonian = wannier_hamiltonian(system, k)
      w_connection = wannier_berry_connection(system, k)
      w_dk_hamiltonian = wannier_dhamiltonian_dk(system, k)

      !Get eigenvalues and rotation.
      call utility_diagonalize(w_hamiltonian, system%num_bands, eig, rot, error)
      if (error) then
        write (unit=113, fmt="(a)") "Error in function optical_conductivity when computing the eigenvalues of the Hamiltonian."
        return
      endif

      !Get occupations.
      rho = hamiltonian_occ_matrix(system, eig)
      !Get non-abelian D.
      na_d = non_abelian_d(system, eig, rot, w_dk_hamiltonian)

      !Get connection in the Hamiltonian basis.
      do i = 1, 3
        connection(:, :, i) = matmul(matmul(transpose(conjg(rot)), w_connection(:, :, i)), rot) + &
                              cmplx_i*na_d(:, :, i)
      enddo

      !Set default smearing and dk.
      smearing = task%smearing
      dk = 0.0_dp

      !In the case of adaptive smearing and an integral task, get velocities and typical spacing.
      if (task%adpt_smearing) then
        vels = velocities(system, w_dk_hamiltonian, eig, rot, error)
        if (error) then
          write (unit=113, fmt="(a)") "Error in function optical_conductivity when computing the velocities."
          return
        endif
        dk = (1.0_dp/(system%cell_volume*product(task%samples)))**(1.0_dp/3.0_dp) !Typical spacing = (inverse cell volume/samples)^(1/3).
      endif

      do i_mem = 1, product(task%integer_indices)

        if ((task%particular_integer_component .ne. 0) .and. (i_mem .ne. task%particular_integer_component)) cycle

        i_arr = integer_memory_element_to_array_element(task, i_mem)
        i = i_arr(1)
        j = i_arr(2)

        do n = 1, system%num_bands
          do m = 1, system%num_bands

            if (n == m) cycle
            if (eig(m) > maxval(task%ext_var_data(1)%data(:)) .or. eig(n) > maxval(task%ext_var_data(1)%data(:))) cycle

            bpart = (rho(n, n) - rho(m, m))*(eig(n) - eig(m))* & !TODO: CHECK THIS FORMULA.
                    connection(m, n, i)*connection(n, m, j)

            if (task%adpt_smearing) then
              spacing = sqrt(2.0_dp)*sqrt(sum(real(vels(n, n, :) - vels(m, m, :), dp)*real(vels(n, n, :) - vels(m, m, :), dp)))*dk
              smearing = min(spacing, task%smearing)
            endif

            do r_mem = 1, product(task%continuous_indices)

              omega = task%ext_var_data(1)%data(r_mem)
              arg = (eig(m) - eig(n) - omega)/smearing
              if (abs(arg) > 10.0_dp) cycle
              delta = utility_delta(arg)/smearing

              u(i_mem, r_mem) = u(i_mem, r_mem) - pi*bpart*delta

            enddo!omega

          enddo!m
        enddo!n

      enddo!ij
      !At this point the units of u are A^2. We have to
      !i) divide by the cell volume in A^3,
      !ii) multiply by 10^10 to pass to 1/m,
      !iii) multiply by e^2/\hbar to get the result in Amperes/(Volt*Meter)
      u = u*1.0E10*e_charge/(hbar_over_e*system%cell_volume)

    end select
  end function optical_conductivity

  subroutine default_jdos_constructor(optical_task, method, samples, &
                                      omegastart, omegaend, omegasteps, &
                                      optical_smearing)

    character(len=*), optional, intent(in) :: method
    integer, optional, intent(in) :: samples(3)

    real(kind=dp), intent(in) :: omegastart, omegaend
    integer, intent(in) :: omegasteps

    real(kind=dp), optional, intent(in) :: optical_smearing

    type(optical_BZ_integral_task), intent(out) :: optical_task

    call BZ_integral_task_constructor(task=optical_task, name="jdos", &
                                      g_calculator=jdos, &
                                      method=method, samples=samples, &
                                      N_int_ind=1, int_ind_range=(/1/), &
                                      N_ext_vars=1, &
                                      ext_vars_start=(/omegastart/), &
                                      ext_vars_end=(/omegaend/), &
                                      ext_vars_steps=(/omegasteps/))

    if (present(optical_smearing)) then
      optical_task%adpt_smearing = .false.
      optical_task%smearing = optical_smearing
    endif

  end subroutine default_jdos_constructor

  function jdos(task, system, k, error) result(u)
    class(global_k_data), intent(in) :: task
    type(sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp)                   :: u(product(task%integer_indices), product(task%continuous_indices))

    integer :: r_mem
    integer :: n, m

    real(kind=dp) :: omega, eig(system%num_bands), smearing, arg, delta, bpart, spacing, dk
    complex(kind=dp), dimension(system%num_bands, system%num_bands) :: &
      w_hamiltonian, &
      rho, &
      rot
    complex(kind=dp), dimension(system%num_bands, system%num_bands, 3) :: w_dk_hamiltonian, vels

    select type (task)
    type is (optical_BZ_integral_task)

      u = cmplx(0.0_dp, 0.0_dp)

      !Gather required quantities in the Wannier basis.
      w_hamiltonian = wannier_hamiltonian(system, k)
      w_dk_hamiltonian = wannier_dhamiltonian_dk(system, k)

      !Get eigenvalues and rotation.
      call utility_diagonalize(w_hamiltonian, system%num_bands, eig, rot, error)
      if (error) then
        write (unit=113, fmt="(a)") "Error in function jdos when computing the eigenvalues of the Hamiltonian."
        return
      endif

      !Get occupations.
      rho = hamiltonian_occ_matrix(system, eig)

      !Set default smearing and dk.
      smearing = task%smearing
      dk = 0.0_dp

      !In the case of adaptive smearing and an integral task, get velocities and typical spacing.
      if (task%adpt_smearing) then
        vels = velocities(system, w_dk_hamiltonian, eig, rot, error)
        if (error) then
          write (unit=113, fmt="(a)") "Error in function jdos when computing the velocities."
          return
        endif
        dk = (1.0_dp/(system%cell_volume*product(task%samples)))**(1.0_dp/3.0_dp) !Typical spacing = (inverse cell volume/samples)^(1/3).
      endif

      do n = 1, system%num_bands
        do m = 1, system%num_bands

          if (n == m) cycle
          if (eig(m) > maxval(task%ext_var_data(1)%data(:)) .or. eig(n) > maxval(task%ext_var_data(1)%data(:))) cycle

          bpart = rho(n, n)*(1.0_dp - rho(m, m))

          if (task%adpt_smearing) then
            spacing = sqrt(2.0_dp)*sqrt(sum(real(vels(n, n, :) - vels(m, m, :), dp)*real(vels(n, n, :) - vels(m, m, :), dp)))*dk
            smearing = min(spacing, task%smearing)
          endif

          do r_mem = 1, product(task%continuous_indices)

            omega = task%ext_var_data(1)%data(r_mem)
            arg = (eig(m) - eig(n) - omega)/smearing
            if (abs(arg) > 10.0_dp) cycle
            delta = utility_delta(arg)/smearing

            u(1, r_mem) = u(1, r_mem) + bpart*delta

          enddo!omega

        enddo!m
      enddo!n

      !At this point the units of u are (states)/eV.
    end select
  end function jdos

  subroutine default_shift_current_constructor(optical_task, method, samples, &
                                               omegastart, omegaend, omegasteps, &
                                               particular_integer_component, &
                                               optical_smearing)

    character(len=*), optional, intent(in) :: method
    integer, optional, intent(in) :: samples(3)

    real(kind=dp), intent(in) :: omegastart, omegaend
    integer, intent(in) :: omegasteps

    integer, optional, intent(in) :: particular_integer_component(:)
    real(kind=dp), optional, intent(in) ::  optical_smearing

    type(optical_BZ_integral_task), intent(out) :: optical_task

    call BZ_integral_task_constructor(task=optical_task, name="sc", &
                                      g_calculator=shift_current, &
                                      method=method, samples=samples, &
                                      N_int_ind=3, int_ind_range=(/3, 3, 3/), &
                                      N_ext_vars=1, &
                                      ext_vars_start=(/omegastart/), &
                                      ext_vars_end=(/omegaend/), &
                                      ext_vars_steps=(/omegasteps/), &
                                      part_int_comp=particular_integer_component)

    if (present(optical_smearing)) then
      optical_task%adpt_smearing = .false.
      optical_task%smearing = optical_smearing
    endif

  end subroutine default_shift_current_constructor

  function shift_current(task, system, k, error) result(u)
    class(global_k_data), intent(in) :: task
    type(sys), intent(in) :: system
    real(kind=dp), intent(in) :: k(3)
    logical, intent(inout) :: error

    complex(kind=dp)                   :: u(product(task%integer_indices), product(task%continuous_indices))

    integer :: i_mem, i_arr(product(task%integer_indices)), &
               r_mem
    integer :: i, j, l, n, m

    real(kind=dp) :: omega, eig(system%num_bands), smearing, arg1, arg2, delta1, delta2, bpart, spacing, dk
    complex(kind=dp), dimension(system%num_bands, system%num_bands) :: &
      w_hamiltonian, &
      rho, &
      rot
    complex(kind=dp), dimension(system%num_bands, system%num_bands, 3) :: &
      w_dk_hamiltonian, &
      w_connection, connection, &
      na_d, vels
    complex(kind=dp), dimension(system%num_bands, system%num_bands, 3, 3) :: &
      w_d2k_hamiltonian, &
      w_dk_connection, &
      gen_r

    select type (task)
    type is (optical_BZ_integral_task)

      u = cmplx(0.0_dp, 0.0_dp)

      !Gather required quantities in the Wannier basis.
      w_hamiltonian = wannier_hamiltonian(system, k)
      w_connection = wannier_berry_connection(system, k)
      w_dk_hamiltonian = wannier_dhamiltonian_dk(system, k)
      w_dk_connection = wannier_dberry_connection_dk(system, k)
      w_d2k_hamiltonian = wannier_d2hamiltonian_dk2(system, k)

      !Get eigenvalues and rotation.
      call utility_diagonalize(w_hamiltonian, system%num_bands, eig, rot, error)
      if (error) then
        write (unit=113, fmt="(a)") "Error in function shift_current when computing the eigenvalues of the Hamiltonian."
        return
      endif

      !Get occupations.
      rho = hamiltonian_occ_matrix(system, eig)
      !Get non-abelian D.
      na_d = non_abelian_d(system, eig, rot, w_dk_hamiltonian)

      !Get connection in the Hamiltonian basis.
      do i = 1, 3
        connection(:, :, i) = matmul(matmul(transpose(conjg(rot)), w_connection(:, :, i)), rot) + &
                              cmplx_i*na_d(:, :, i)
      enddo
      !Get covariant derivative of the dipoles in the Hamilonian basis.
      gen_r = cov_deriv_of_dipole(system, w_d2k_hamiltonian, w_dk_hamiltonian, w_connection, w_dk_connection, eig, rot, error)
      if (error) then
        write (unit=113, fmt="(a)") "Error in function shift_current when computing the generalized derivative of the dipoles."
        return
      endif

      !Set default smearing and dk.
      smearing = task%smearing
      dk = 0.0_dp

      !In the case of adaptive smearing and an integral task, get velocities and typical spacing.
      if (task%adpt_smearing) then
        vels = velocities(system, w_dk_hamiltonian, eig, rot, error)
        if (error) then
          write (unit=113, fmt="(a)") "Error in function shift_current when computing the velocities."
          return
        endif
        dk = (1.0_dp/(system%cell_volume*product(task%samples)))**(1.0_dp/3.0_dp) !Typical spacing = (inverse cell volume/samples)^(1/3).
      endif

      do i_mem = 1, product(task%integer_indices)

        if ((task%particular_integer_component .ne. 0) .and. (i_mem .ne. task%particular_integer_component)) cycle

        i_arr = integer_memory_element_to_array_element(task, i_mem)!TODO: ACCOUNT FOR SYMMETRIZATION.
        i = i_arr(1)
        j = i_arr(2)
        l = i_arr(3)

        do n = 1, system%num_bands
          do m = 1, system%num_bands

            if (n == m) cycle
            !if (eig(m) > maxval(task%ext_var_data(1)%data(:)) .or. eig(n) > maxval(task%ext_var_data(1)%data(:))) cycle

            bpart = (rho(n, n) - rho(m, m))* &
                    (connection(m, n, j)*gen_r(n, m, i, l) + connection(m, n, l)*gen_r(n, m, i, j))

            if (task%adpt_smearing) then
              spacing = sqrt(2.0_dp)*sqrt(sum(real(vels(n, n, :) - vels(m, m, :), dp)*real(vels(n, n, :) - vels(m, m, :), dp)))*dk
              smearing = min(spacing, task%smearing)
            endif

            do r_mem = 1, product(task%continuous_indices)

              omega = task%ext_var_data(1)%data(r_mem)

              arg1 = (eig(m) - eig(n) - omega)/smearing
              delta1 = utility_delta(arg1)/smearing

              arg2 = (eig(n) - eig(m) - omega)/smearing
              delta2 = utility_delta(arg2)/smearing

              u(i_mem, r_mem) = u(i_mem, r_mem) - cmplx_i*0.25_dp*pi*bpart*(delta1 + delta2)

            enddo!omega

          enddo!m
        enddo!n

      enddo!ijl
      !At this point the units of u are A^3/eV. We have to
      !i) divide by the cell volume in A^3,
      !ii) multiply by hbar/e to pass from 1/eV to s,
      !iii) multiply by e^3/\hbar^2 = e/(hbar_over_e^2) to get the result in Amperes/(Volt^2)
      !u = u*e_charge/((hbar_over_e)*system%cell_volume) !TODO: test these things, it appears to be some underflow

    end select
  end function shift_current

end module calculators_optical
