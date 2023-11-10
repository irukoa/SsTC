program testing_driver

  use, intrinsic :: iso_fortran_env, only: error_unit

  use OMP_LIB
  use MPI_F08

  use SsTC

  use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
  & select_suite, run_selected, get_argument

  use suite_l_k_quantities, only: collect_suite_l_k_quantities
  use suite_integrator, only: collect_suite_integrator
  use suite_sampler, only: collect_suite_sampler
  use suite_kpath, only: collect_suite_kpath
  use suite_kslice, only: collect_suite_kslice

  implicit none

  integer :: ierror, nProcs, rank

  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  call MPI_INIT(ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
  if (rank == 0) write (error_unit, "(a, i5, a)") "Running on ", nProcs, " MPI processes."

  call OMP_SET_MAX_ACTIVE_LEVELS(1)

  call SsTC_init(verb=.true., nThreads=1, exec_label="s_test")

  stat = 0

  testsuites = [ &
               new_testsuite("Local k Quantities", collect_suite_l_k_quantities), &
               new_testsuite("Integrator", collect_suite_integrator), &
               new_testsuite("Sampler", collect_suite_sampler), &
               new_testsuite("Kpath", collect_suite_kpath), &
               new_testsuite("Kslice", collect_suite_kslice) &
               ]

  call get_argument(1, suite_name)
  call get_argument(2, test_name)

  if (allocated(suite_name)) then
    is = select_suite(testsuites, suite_name)
    if (is > 0 .and. is <= size(testsuites)) then
      if (allocated(test_name)) then
        write (error_unit, fmt) "Suite:", testsuites(is)%name
        call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
        if (stat < 0) then
          error stop 1
        end if
      else
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
      end if
    else
      write (error_unit, fmt) "Available testsuites"
      do is = 1, size(testsuites)
        write (error_unit, fmt) "-", testsuites(is)%name
      end do
      error stop 1
    end if
  else
    do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
    end do
  end if

  if (stat > 0) then
    write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

  if (nProcs == 1) then

    if (rank == 0) write (error_unit, "(a, i5, a)") "Repeating tests in multithread case..."
    call SsTC_init(verb=.true., nThreads=2, exec_label="p_test")

    call get_argument(1, suite_name)
    call get_argument(2, test_name)

    if (allocated(suite_name)) then
      is = select_suite(testsuites, suite_name)
      if (is > 0 .and. is <= size(testsuites)) then
        if (allocated(test_name)) then
          write (error_unit, fmt) "Suite:", testsuites(is)%name
          call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
          if (stat < 0) then
            error stop 1
          end if
        else
          write (error_unit, fmt) "Testing:", testsuites(is)%name
          call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
        end if
      else
        write (error_unit, fmt) "Available testsuites"
        do is = 1, size(testsuites)
          write (error_unit, fmt) "-", testsuites(is)%name
        end do
        error stop 1
      end if
    else
      do is = 1, size(testsuites)
        write (error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
      end do
    end if

    if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
    end if

  endif

  call MPI_FINALIZE(ierror)

end program testing_driver
