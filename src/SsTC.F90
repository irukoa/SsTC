#include "cond_comp.h"
module SsTC

  USE OMP_LIB
  USE MPI_F08

  use SsTC_utility

  use SsTC_mpi_comms

  use extrapolation_integration, only: SsTC_integral_extrapolation => integral_extrapolation, &
    SsTC_shrink_array => shrink_array, &
    SsTC_expand_array => expand_array

  use SsTC_data_structures

  use SsTC_kpath
  use SsTC_kslice
  use SsTC_sampler
  use SsTC_integrator

  use SsTC_local_k_quantities

  !###MOD HEADERS

  implicit none

  private

  !Core Procedures:
  !Initialization.
  public :: SsTC_init

  !Utility.
  public :: SsTC_alpha_S
  public :: SsTC_beta_S

  public :: SsTC_alpha_A
  public :: SsTC_beta_A

  public :: SsTC_utility_delta
  public :: SsTC_utility_delta_vec
  public :: SsTC_utility_get_degen
  public :: SsTC_utility_diagonalize
  public :: SsTC_utility_schur
  public :: SsTC_utility_SVD
  public :: SsTC_utility_exphs
  public :: SsTC_utility_logu

  !Extrapolation integration.
  public :: SsTC_integral_extrapolation
  public :: SsTC_shrink_array
  public :: SsTC_expand_array

  !Data structures.
  public :: SsTC_local_k_data
  public :: SsTC_global_k_data

  public :: SsTC_local_calculator
  public :: SsTC_global_calculator

  public :: SsTC_sys
  public :: SsTC_sys_constructor

  public :: SsTC_external_vars
  public :: SsTC_external_variable_constructor

  public :: SsTC_integer_array_element_to_memory_element
  public :: SsTC_integer_memory_element_to_array_element
  public :: SsTC_continuous_array_element_to_memory_element
  public :: SsTC_continuous_memory_element_to_array_element
  public :: SsTC_construct_iterable

  !Kpath.
  public :: SsTC_kpath_task

  public :: SsTC_kpath_constructor
  public :: SsTC_kpath_sampler
  public :: SsTC_print_kpath

  !Kslice.
  public :: SsTC_kslice_task

  public :: SsTC_kslice_task_constructor
  public :: SsTC_sample_kslice_task
  public :: SsTC_print_kslice

  !BZ Sampler.
  public :: SsTC_sampling_task

  public :: SsTC_sampling_task_constructor
  public :: SsTC_sample_sampling_task
  public :: SsTC_print_sampling

  !Integrator.
  public :: SsTC_BZ_integral_task

  public :: SsTC_BZ_integral_task_constructor
  public :: SsTC_sample_and_integrate_BZ_integral_task
  public :: SsTC_print_BZ_integral_task

  !Local k-dependent Quantities.
  public :: SsTC_wannier_hamiltonian
  public :: SsTC_wannier_berry_connection
  public :: SsTC_wannier_dhamiltonian_dk
  public :: SsTC_wannier_d2hamiltonian_dk2
  public :: SsTC_wannier_dberry_connection_dk
  public :: SsTC_hamiltonian_occ_matrix
  public :: SsTC_non_abelian_d
  public :: SsTC_velocities
  public :: SsTC_inverse_effective_mass

  public :: SsTC_get_hamiltonian
  public :: SsTC_get_position

  !Non-Core Procedures:
  !###MOD PROCEDURES

contains

  subroutine SsTC_init(nThreads, nNested, exec_label)

    implicit none

    integer, intent(in), optional :: nThreads
    integer, intent(in), optional :: nNested
    character(len=*), intent(in), optional :: exec_label

    integer :: selThreads

    call MPI_INITIALIZED(is_mpi_initialized, ierror)
    call MPI_FINALIZED(is_mpi_finalized, ierror)

    if (.not. ((is_mpi_initialized) .and. ((.not. is_mpi_finalized)))) then
      write (unit=stdout, fmt="(a)") "          SsTC: MPI has not bee initialized or has been finalized. Stopping..."
      stop
    endif

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    if (present(exec_label)) then
      if (rank == 0) open (newunit=stdout, action="write", file=trim(exec_label//".out"))
      if (rank == 0) open (newunit=stderr, action="write", file=trim(exec_label//".err"))
    else
      if (rank == 0) open (newunit=stdout, action="write", file="SsTC_exec.out")
      if (rank == 0) open (newunit=stderr, action="write", file="SsTC_exec.err")
    endif

    call date_and_time(values=timing)

    if (rank == 0) write (unit=stdout, fmt="(a)") "            /$$$$$$            /$$$$$$$$  /$$$$$$  "
    if (rank == 0) write (unit=stdout, fmt="(a)") "           /$$__  $$          |__  $$__/ /$$__  $$ "
    if (rank == 0) write (unit=stdout, fmt="(a)") "          | $$  \__/  /$$$$$$$   | $$   | $$  \__/ "
    if (rank == 0) write (unit=stdout, fmt="(a)") "          |  $$$$$$  /$$_____/   | $$   | $$       "
    if (rank == 0) write (unit=stdout, fmt="(a)") "           \____  $$|  $$$$$$    | $$   | $$       "
    if (rank == 0) write (unit=stdout, fmt="(a)") "           /$$  \ $$ \____  $$   | $$   | $$    $$ "
    if (rank == 0) write (unit=stdout, fmt="(a)") "          |  $$$$$$/ /$$$$$$$/   | $$   |  $$$$$$/ "
    if (rank == 0) write (unit=stdout, fmt="(a)") "           \______/ |_______/    |__/    \______/  "

    if (rank == 0) write (unit=stdout, fmt="(a)") ""
    if (rank == 0) write (unit=stdout, fmt="(a, a, a)") "          Version: ", trim(_VERSION)//"."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

    if (rank == 0) write (unit=stdout, &
                          fmt="(a, i2, a, i2, a, i4, a, i2, a, i2, a, i2, a)") &
      "          SsTC library initializing at ", &
      timing(2), "/", timing(3), "/", timing(1), &
      ", ", timing(5), ":", timing(6), &
      ":", timing(7), "."
    if (rank == 0) write (unit=stdout, fmt="(a, i5, a)") "          Running ", nProcs, " MPI processes."

    if (present(nThreads)) then
      if ((nThreads > 0)) then
        call OMP_SET_NUM_THREADS(nThreads)
        selThreads = nThreads
      endif
    else
      call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())
      selThreads = OMP_GET_MAX_THREADS()
    endif
    if (rank == 0) write (unit=stdout, fmt="(a, i5, a)") "          Paralell regions will run in ", selThreads, " threads."

    if (present(nNested)) then
      if ((nNested > 1)) then
        call OMP_SET_MAX_ACTIVE_LEVELS(nNested)
        if (rank == 0) write (unit=stdout, fmt="(a)") "          The number of nested active parallel regions"
        if (rank == 0) write (unit=stdout, fmt="(a, i2, a)") "          has been set to ", nNested, "."
      endif
    else
      call OMP_SET_MAX_ACTIVE_LEVELS(1)
      if (rank == 0) write (unit=stdout, fmt="(a)") "          The number of nested active parallel regions"
      if (rank == 0) write (unit=stdout, fmt="(a)") "          has been set to 1."
    endif

    if (rank == 0) write (unit=stdout, fmt="(a)") "          SsTC initialized."
    if (rank == 0) write (unit=stdout, fmt="(a)") ""

  end subroutine SsTC_init

end module SsTC
