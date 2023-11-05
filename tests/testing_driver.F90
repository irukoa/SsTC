program testing_driver

  use, intrinsic :: iso_fortran_env, only: error_unit

  use OMP_LIB
  use MPI_F08

  use SsTC

  use testdrive, only: run_testsuite, new_testsuite, testsuite_type, &
  & select_suite, run_selected, get_argument

  use suite_l_k_quantities, only: collect_suite_l_k_quantities

  implicit none

  integer, parameter :: dp = 8

  real(dp), parameter :: pi = acos(-1.0_dp)

  integer :: ierror

  integer :: stat, is
  character(len=:), allocatable :: suite_name, test_name
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  call MPI_INIT(ierror)

  call OMP_SET_MAX_ACTIVE_LEVELS(1)

  call SsTC_init(verb=.false.)

  stat = 0

  testsuites = [ &
               new_testsuite("suite1", collect_suite_l_k_quantities) &
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
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
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
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do
  end if

  if (stat > 0) then
    write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop 1
  end if

  !call SsTC_init()

  call MPI_FINALIZE(ierror)

end program testing_driver
