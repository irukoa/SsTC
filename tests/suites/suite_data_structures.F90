module suite_data_structures

  use SsTC

  use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed

  implicit none

  integer, parameter, private :: dp = 8

  private

  public :: collect_suite_data_structures

contains

  subroutine collect_suite_data_structures(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("Test System Constructor", test_cons), &
                new_unittest("Test External Variable Constructor", test_ext), &
                new_unittest("Test Integer Array Layout to Memory Layout", test_iam), &
                new_unittest("Test Integer Memory Layout to Array Layout", test_ima), &
                new_unittest("Test Continuous Array Layout to Memory Layout", test_cam), &
                new_unittest("Test Continuous Memory Layout to Array Layout", test_cma), &
                new_unittest("Test Iterable Dictionary", test_iterable) &
                ]

  end subroutine collect_suite_data_structures

  subroutine test_cons(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    type(SsTC_sys) :: SsTC_toy_model

    SsTC_toy_model = SsTC_sys_constructor(name="toy_model", path_to_tb_file="./data/", efermi=0.0_dp)

    if (abs(SsTC_toy_model%cell_volume - 1.0_dp) > 1.0E-6_dp) cond = .false.
    if (abs(real(SsTC_toy_model%real_space_hamiltonian_elements(1, 1, 1), dp) + 0.5_dp) > 1.0E-6_dp) cond = .false.
    if (abs(real(SsTC_toy_model%real_space_hamiltonian_elements(1, 2, 1), dp) - 3.0_dp) > 1.0E-6_dp) cond = .false.
    if (abs(real(SsTC_toy_model%real_space_position_elements(1, 1, 1, 1), dp) + 0.5_dp) > 1.0E-6_dp) cond = .false.
    if (abs(real(SsTC_toy_model%real_space_position_elements(1, 2, 1, 1), dp) - 3.0_dp) > 1.0E-6_dp) cond = .false.

    if (.not. cond) then
      call test_failed(error, "Mismatch between the allocated data &
                               &and the reference.")
      return
    end if

  end subroutine test_cons

  subroutine test_ext(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    block

      type(SsTC_external_vars) :: ext

      allocate (ext%data(100))

      ext = SsTC_external_variable_constructor(0.0_dp, 1.0_dp, 100)

      if (abs(ext%data(25) - (24.0_dp/99.0_dp)) > 1.0E-6_dp) cond = .false.

    end block

    if (.not. cond) then
      call test_failed(error, "Mismatch between the allocated data &
                               &and the reference.")
      return
    end if

  end subroutine test_ext

  subroutine test_iam(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    block

      type(SsTC_local_k_data) :: ldata
      integer :: a, b, c, d, e, im, expl

      allocate (ldata%integer_indices(5))
      ldata%integer_indices = (/1, 2, 3, 4, 5/)

      expl = 0
      do a = 1, ldata%integer_indices(1)
        do b = 1, ldata%integer_indices(2)
          do c = 1, ldata%integer_indices(3)
            do d = 1, ldata%integer_indices(4)
              do e = 1, ldata%integer_indices(5)
                expl = expl + 1
                im = SsTC_integer_array_element_to_memory_element(ldata, (/a, b, c, d, e/))
                if (.not. (im == expl)) cond = .false.
              enddo
            enddo
          enddo
        enddo
      enddo

    end block

    if (.not. cond) then
      call test_failed(error, "Mismatch between some iteration of the  &
                               &jagged array and the reference.")
      return
    end if

  end subroutine test_iam

  subroutine test_ima(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    block

      type(SsTC_local_k_data) :: ldata
      integer :: a, b, c, d, e, iarr(5), expl(5), im

      allocate (ldata%integer_indices(5))
      ldata%integer_indices = (/1, 2, 3, 4, 5/)

      im = 0
      do a = 1, ldata%integer_indices(1)
        do b = 1, ldata%integer_indices(2)
          do c = 1, ldata%integer_indices(3)
            do d = 1, ldata%integer_indices(4)
              do e = 1, ldata%integer_indices(5)
                im = im + 1
                expl = (/a, b, c, d, e/)
                iarr = SsTC_integer_memory_element_to_array_element(ldata, im)
                if ((.not. (iarr(1) == expl(1))) .or. &
                    (.not. (iarr(2) == expl(2))) .or. &
                    (.not. (iarr(3) == expl(3))) .or. &
                    (.not. (iarr(4) == expl(4))) .or. &
                    (.not. (iarr(5) == expl(5)))) cond = .false.
              enddo
            enddo
          enddo
        enddo
      enddo

    end block

    if (.not. cond) then
      call test_failed(error, "Mismatch between some iteration of the  &
                               &jagged array and the reference.")
      return
    end if

  end subroutine test_ima

  subroutine test_cam(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    block

      type(SsTC_global_k_data) :: gdata

      integer :: a, b, rm, expl

      allocate (gdata%continuous_indices(2))
      gdata%continuous_indices = (/10, 10/)

      expl = 0
      do a = 1, gdata%continuous_indices(1)
        do b = 1, gdata%continuous_indices(2)
          expl = expl + 1
          rm = SsTC_continuous_array_element_to_memory_element(gdata, (/a, b/))
          if (.not. (rm == expl)) cond = .false.
        enddo
      enddo

    end block

    if (.not. cond) then
      call test_failed(error, "Mismatch between some iteration of the  &
                               &jagged array and the reference.")
      return
    end if

  end subroutine test_cam

  subroutine test_cma(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    block

      type(SsTC_global_k_data) :: gdata
      integer :: a, b, rarr(2), expl(2), rm

      allocate (gdata%continuous_indices(2))
      gdata%continuous_indices = (/10, 10/)

      rm = 0
      do a = 1, gdata%continuous_indices(1)
        do b = 1, gdata%continuous_indices(2)
          rm = rm + 1
          expl = (/a, b/)
          rarr = SsTC_continuous_memory_element_to_array_element(gdata, rm)
          if ((.not. (rarr(1) == expl(1))) .or. &
              (.not. (rarr(2) == expl(2)))) cond = .false.
        enddo
      enddo

    end block

    if (.not. cond) then
      call test_failed(error, "Mismatch between some iteration of the  &
                               &jagged array and the reference.")
      return
    end if

  end subroutine test_cma

  subroutine test_iterable(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: cond = .true.

    block

      type(SsTC_global_k_data) :: gdata

      allocate (gdata%continuous_indices(3))
      gdata%continuous_indices = (/5, 5, 5/)
      allocate (gdata%ext_var_data(3))
      gdata%ext_var_data(1) = SsTC_external_variable_constructor(0.0_dp, 1.0_dp, 5)
      gdata%ext_var_data(2) = SsTC_external_variable_constructor(0.0_dp, 1.0_dp, 5)
      gdata%ext_var_data(3) = SsTC_external_variable_constructor(0.0_dp, 1.0_dp, 5)

      call SsTC_construct_iterable(gdata, vars=(/1, 3/))

      if ((.not. (gdata%iterables(22, 1) == 5)) .or. &
          (.not. (gdata%iterables(22, 2) == 1)) .or. &
          (.not. (gdata%iterables(22, 3) == 2))) cond = .false.

    end block

    if (.not. cond) then
      call test_failed(error, "Mismatch between some iteration of the  &
                               &jagged array and the reference.")
      return
    end if

  end subroutine test_iterable

end module suite_data_structures
