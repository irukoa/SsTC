module SsTC

  USE OMP_LIB

  use utility
  use data_structures
  use extrapolation_integration
  use local_k_quantities
  use kpath
  use kslice
  use integrator
  use calculators_general
  use calculators_floquet
  use calculators_optical

contains

  subroutine init_SsTC(nThreads, nNested, exec_label)

    integer, intent(in), optional :: nThreads
    integer, intent(in), optional :: nNested
    character(len=*), intent(in), optional :: exec_label

    integer :: selThreads

    if (present(exec_label)) then
      open (newunit=stdout, action="write", file=trim(exec_label//".out"))
      open (newunit=stderr, action="write", file=trim(exec_label//".err"))
    else
      open (newunit=stdout, action="write", file="SsTC_exec.out")
      open (newunit=stderr, action="write", file="SsTC_exec.err")
    endif

    write (unit=stdout, fmt="(a)") "Initializing SsTC..."

    if (present(nThreads) .and. (nThreads > 0)) then
      call OMP_SET_NUM_THREADS(nThreads)
      selThreads = nThreads
    else
      call OMP_SET_NUM_THREADS(OMP_GET_MAX_THREADS())
      selThreads = OMP_GET_MAX_THREADS()
    endif

    write (unit=stdout, fmt="(a, i5, a)") "Paralell regions will run in ", selThreads, " threads."

    if (present(nNested) .and. (nNested > 1)) then
      call OMP_SET_MAX_ACTIVE_LEVELS(nNested)
      write (unit=stdout, fmt="(a, i2, a)") "Warning: the number of nested active parallel regions&
      & has been set to ", nNested, ". Beware of overhead."
    else
      call OMP_SET_MAX_ACTIVE_LEVELS(1)
      write (unit=stdout, fmt="(a)") "The number of nested active parallel regions has been set to 1."
    endif

    write (unit=stdout, fmt="(a)") "SsTC initialized."
    write (unit=stdout, fmt="(a)") ""

  end subroutine init_SsTC

end module SsTC
