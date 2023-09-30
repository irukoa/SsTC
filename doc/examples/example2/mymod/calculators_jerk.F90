module mymod

  use SsTC_utility
  use SsTC_data_structures
  use SsTC_local_k_quantities
  use SsTC_integrator

  implicit none

  private

  type, extends(SsTC_BZ_integral_task) :: optical_BZ_integral_task
    logical       :: adpt_smearing = .true.
    real(kind=dp) :: smearing = 1.0_dp
  end type optical_BZ_integral_task

  public :: optical_BZ_integral_task

  public :: jerk_current_constructor
  public :: jerk_current

contains

  subroutine jerk_current_constructor(optical_task, method, samples, &
                                                   omegastart, omegaend, omegasteps, &
                                                   particular_integer_component, &
                                                   optical_smearing)

    character(len=*), optional, intent(in) :: method
    integer, optional, intent(in)          :: samples(3)

    real(kind=dp), intent(in) :: omegastart, omegaend
    integer, intent(in)       :: omegasteps

    integer, optional, intent(in)       :: particular_integer_component(:)
    real(kind=dp), optional, intent(in) ::  optical_smearing

    type(optical_BZ_integral_task), intent(out) :: optical_task

    call SsTC_BZ_integral_task_constructor(task=optical_task, name="jc", &
                                           g_calculator=jerk_current, &
                                           method=method, samples=samples, &
                                           N_int_ind=4, int_ind_range=(/3, 3, 3, 3/), &
                                           N_ext_vars=1, &
                                           ext_vars_start=(/omegastart/), &
                                           ext_vars_end=(/omegaend/), &
                                           ext_vars_steps=(/omegasteps/), &
                                           part_int_comp=particular_integer_component)

    if (present(optical_smearing)) then
      optical_task%adpt_smearing = .false.
      optical_task%smearing = optical_smearing
    endif

  end subroutine jerk_current_constructor

  function jerk_current(task, system, k, error) result(u)
    class(SsTC_global_k_data), intent(in) :: task
    type(SsTC_sys), intent(in)            :: system
    real(kind=dp), intent(in)             :: k(3)
    logical, intent(inout)                :: error

    complex(kind=dp) :: u(product(task%integer_indices), product(task%continuous_indices))

    integer :: i_mem, i_arr(size(task%integer_indices)), &
               r_mem
    integer :: i, j, l, q, &
               iq, jl, &
               n, m

    real(kind=dp) :: omega, smearing, &
                     spacing, dk, &
                     arg, delta, &
                     bpart, &
                     eig(system%num_bands)

    complex(kind=dp), dimension(system%num_bands, system%num_bands) :: w_hamiltonian, &
                                                                       rho, &
                                                                       rot

    complex(kind=dp), dimension(system%num_bands, system%num_bands, 3)    :: w_dk_hamiltonian, &
                                                                             w_connection, connection, &
                                                                             na_d, vels

    complex(kind=dp), dimension(system%num_bands, system%num_bands, 3, 3) :: mu

    complex(kind=dp), dimension(product(task%continuous_indices)) :: symcmp

    u = cmplx_0

    select type (task)
    type is (optical_BZ_integral_task)

      !Gather required quantities in the Wannier basis.
      w_hamiltonian = SsTC_wannier_hamiltonian(system, k)
      w_connection = SsTC_wannier_berry_connection(system, k)
      w_dk_hamiltonian = SsTC_wannier_dhamiltonian_dk(system, k)

      !Get eigenvalues and rotation.
      call SsTC_utility_diagonalize(w_hamiltonian, system%num_bands, eig, rot, error)

      !Get occupations.
      rho = SsTC_hamiltonian_occ_matrix(system, eig)
      !Get non-abelian D.
      na_d = SsTC_non_abelian_d(system, eig, rot, w_dk_hamiltonian)
      !Get inverse effective mass.
      mu = SsTC_inverse_effective_mass(system, SsTC_wannier_d2hamiltonian_dk2(system, k), &
                                  w_dk_hamiltonian, eig, rot, error)

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
        vels = SsTC_deleig(system, w_dk_hamiltonian, eig, rot, error)
        !Typical spacing = (inverse of cell volume/samples)^(1/3).
        dk = (1.0_dp/(system%cell_volume*product(task%samples)))**(1.0_dp/3.0_dp)
      endif

      do iq = 1, 6
        do jl = 1, 6

          i = SsTC_alpha_S(iq)
          q = SsTC_beta_S(iq)
          j = SsTC_alpha_S(jl)
          l = SsTC_beta_S(jl)

          i_arr = (/i, j, l, q/)
          i_mem = SsTC_integer_array_element_to_memory_element(task, i_arr)

          if ((task%particular_integer_component .ne. 0) .and. (i_mem .ne. task%particular_integer_component)) cycle

          do n = 1, system%num_bands
            do m = 1, system%num_bands

              if (n == m) cycle

              bpart = (rho(n, n) - rho(m, m))* &
                      (connection(n, m, j)*connection(m, n, l))* &
                      (mu(n, n, i, q) - mu(m, m, i, q))

              if (task%adpt_smearing) then
                spacing = sqrt(2.0_dp)*sqrt(sum(real(vels(n, n, :) - vels(m, m, :), dp)*real(vels(n, n, :) - vels(m, m, :), dp)))*dk
                smearing = min(spacing, task%smearing)
              endif

              do r_mem = 1, product(task%continuous_indices)

                omega = task%ext_var_data(1)%data(r_mem)

                arg = (eig(m) - eig(n) - omega)/smearing
                delta = SsTC_utility_delta(arg)/smearing

                u(i_mem, r_mem) = u(i_mem, r_mem) + 2.0_dp*pi*bpart*delta

              enddo!omega

            enddo!m
          enddo!n

          symcmp = u(i_mem, :)

          i_arr = (/q, l, j, i/)
          i_mem = SsTC_integer_array_element_to_memory_element(task, i_arr)

          u(i_mem, :) = symcmp

        enddo!jl
      enddo!iq
    end select
    !At this point the units of u are A^4. We have to:
    !i) divide by the cell volume in A^3,
    !ii) multiply bu 10^-10 to get u in Meters,
    !iii) multiply by e^4/\hbar^3 = e/(hbar_over_e^3) to get the result in Amperes*Meters/(Seconds^2*Volt^3).
    u = u*1.0E-10_dp*e_charge/((hbar_over_e**3)*system%cell_volume)
  end function jerk_current

end module mymod
