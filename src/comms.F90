module SsTC_mpi_comms

  implicit none

  private

  logical, public :: is_mpi_initialized = .false., is_mpi_finalized = .false.
  integer, public :: rank = 0, nProcs = 1, ierror = 0 !MPI singleton defaults.

  public :: get_MPI_task_partition

contains

  subroutine get_MPI_task_partition(task_size, nodes, counts, displs)

    implicit none

    integer, intent(in) :: task_size, nodes
    integer, intent(out) :: counts(0:nodes - 1), &
                            displs(0:nodes - 1)

    integer :: ratio, remainder, i

    ratio = task_size/nodes
    remainder = modulo(task_size, nodes)

    do i = 0, nodes - 1
      if (i < remainder) then
        counts(i) = ratio + 1
        displs(i) = i*(ratio + 1)
      else
        counts(i) = ratio
        displs(i) = remainder*(ratio + 1) + (i - remainder)*ratio
      end if
    end do

  end subroutine get_MPI_task_partition

end module SsTC_mpi_comms
