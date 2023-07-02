module extrapolation_integration

  !V1.0.1(modified) By Álvaro R. Puente-Uriona.

  !Standalone module to compute integrals using the extrapolation
  !a.k.a. Romberg method.

  !The main routine, integral_extrapolation takes as inputs:
  !i) a real or complex, scalar* or vector* array in memory layout* "array", containing
  !the values of the integrand within a discretized mesh in d = 1, 2 or 3
  !dimensions.
  !ii) a size(d) integer array "sizes", specifying the number of points used
  !in each dimension for the discretization. Any entry of the array shall be expressible as
  !size(i) = 2^n + 1 for n = 0, 1, 2,... for extrapolation to be possible. The only exception
  !is size(i) = 1, for which there is no integration.
  !iii) a size(2*d) real or complex array "int_bounds" containing the integral bounds in
  !sequential order, such that the discretization of dimension i has been done as
  !x_i= int_bounds(2*i - 1) + [int_bounds(2*i) - int_bounds(2*i - 1)]/[size(i) - 1].
  !And takes as outputs:
  !i) a real or complex, scalar or vector array in memory layout "result", containing
  !the integral values in d = 1, 2 or 3 dimensions.
  !ii) an integer info, specifying the calculation result:
  ! ii.a) info =  1: sucess.
  ! ii.b) info =  0: not suitable for extrapolation, returning the rectangle method approximation.
  ! ii.c) info = -1: error.
  !
  !The routine shrink_array passes an array from arbitrary layout to memory layout and
  !the routine expand_array passes an array from memory layout to arbitrary layout.
  !
  !(*) See documentation.

  implicit none

  integer, parameter, private :: dp = 8

  public :: integral_extrapolation
  public :: shrink_array
  public :: expand_array

  interface integral_extrapolation
    module procedure scalar_integral_extrapolation_real
    module procedure vector_integral_extrapolation_real
    module procedure scalar_integral_extrapolation_complex
    module procedure vector_integral_extrapolation_complex
    module procedure scalar_integral_extrapolation_complex_array_real_bounds
    module procedure vector_integral_extrapolation_complex_array_real_bounds
  end interface integral_extrapolation

  interface shrink_array
    module procedure shrink_rarray1
    module procedure shrink_rarray2
    module procedure shrink_rarray3
    module procedure shrink_rarray4
    module procedure shrink_carray1
    module procedure shrink_carray2
    module procedure shrink_carray3
    module procedure shrink_carray4
  end interface shrink_array

  interface expand_array
    module procedure expand_rarray1
    module procedure expand_rarray2
    module procedure expand_rarray3
    module procedure expand_rarray4
    module procedure expand_carray1
    module procedure expand_carray2
    module procedure expand_carray3
    module procedure expand_carray4
  end interface expand_array

contains

  subroutine scalar_integral_extrapolation_real(array, sizes, int_bounds, result, info)

    real(kind=dp), intent(in)  :: array(:), int_bounds(:)
    integer, intent(in)  :: sizes(:)

    real(kind=dp), intent(out) :: result
    integer, intent(out) :: info

    integer                    :: n1, n2, n3, ep
    real(kind=dp)              :: nr1, nr2, nr3
    real(kind=dp), allocatable :: e3(:, :, :), &
                                  e2(:, :)

    result = 0.0_dp

    if (2*size(sizes) .ne. size(int_bounds)) then
      info = -1 !Error.
      return
    endif

    if (product(sizes) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    if (size(sizes) .eq. 3) then

      !Approximate to 0.
      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(4) - int_bounds(3)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(6) - int_bounds(5)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1 !Sucess.
      !  return

      !endif

      !Calculate n_i = log_2(N_i - 1)
      !Exceptional case: N_i = 1:
      !Assign an otherwise unachievable negative number and use it
      !as a flag for nscalar1_integrate_complex. That dimension will
      !not be integrated over.
      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif
      if (sizes(2) .eq. 1) then
        nr2 = -1.0_dp
        n2 = -1
      else
        nr2 = log(real(sizes(2) - 1, dp))/log(2.0_dp)
        n2 = nint(nr2)
      endif
      if (sizes(3) .eq. 1) then
        nr3 = -1.0_dp
        n3 = -1
      else
        nr3 = log(real(sizes(3) - 1, dp))/log(2.0_dp)
        n3 = nint(nr3)
      endif

      !Check if the mesh is suitable for extrapolation.
      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp) .or. &
          (abs(real(n2, dp) - nr2) .ge. 1.0E-6_dp) .or. &
          (abs(real(n3, dp) - nr3) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))* &
                             (int_bounds(4) - int_bounds(3))* &
                             (int_bounds(6) - int_bounds(5))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        allocate (e3(sizes(1), sizes(2), sizes(3)))
        !i) Pass from memory to arbitrary array,
        call expand_rarray3(array(:), e3(:, :, :), ep)
        !ii) integrate as if the data was scattered from [0, 1]x[0, 1]x[0, 1],
        result = nscalar3_integrate_real(e3(:, :, :), n1, n2, n3)
        deallocate (e3)
        !iii) multiply by integral bounds so the integral is properly normalized.
        result = result*((int_bounds(2) - int_bounds(1))* &
                         (int_bounds(4) - int_bounds(3))* &
                         (int_bounds(6) - int_bounds(5)))
        info = 1 !Sucess.

      endif

    elseif (size(sizes) .eq. 2) then

      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(4) - int_bounds(3)) .le. 1.0E-6_dp)) then

      !  result = 0.0_dp
      !  info = 1
      !  return

      !endif

      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif
      if (sizes(2) .eq. 1) then
        nr2 = -1.0_dp
        n2 = -1
      else
        nr2 = log(real(sizes(2) - 1, dp))/log(2.0_dp)
        n2 = nint(nr2)
      endif

      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp) .or. &
          (abs(real(n2, dp) - nr2) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))* &
                             (int_bounds(4) - int_bounds(3))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        allocate (e2(sizes(1), sizes(2)))
        call expand_rarray2(array(:), e2(:, :), ep)
        result = nscalar2_integrate_real(e2(:, :), n1, n2)
        deallocate (e2)
        result = result*((int_bounds(2) - int_bounds(1))* &
                         (int_bounds(4) - int_bounds(3)))
        info = 1

      endif

    elseif (size(sizes) .eq. 1) then

      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp)) then

      !  result = 0.0_dp
      !  info = 1
      !  return

      !endif

      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif

      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        result = nscalar1_integrate_real(array(:), n1)
        result = result*((int_bounds(2) - int_bounds(1)))
        info = 1

      endif

    else

      info = -1
      return

    endif

  end subroutine scalar_integral_extrapolation_real

  subroutine vector_integral_extrapolation_real(array, sizes, int_bounds, result, info)

    real(kind=dp), intent(in)  :: array(:, :), int_bounds(:)
    integer, intent(in)        :: sizes(:)

    real(kind=dp), intent(out) :: result(size(array(1, :)))
    integer, intent(out)       :: info

    integer                    :: i

    result = 0.0_dp

    if (2*size(sizes) .ne. size(int_bounds)) then
      info = -1 !Error.
      return
    endif

    if (product(sizes) .ne. size(array(:, 1))) then
      info = -1 !Error.
      return
    endif

    do i = 1, size(array(1, :))
      call scalar_integral_extrapolation_real(array(:, i), sizes, int_bounds, result(i), info)
    enddo

  end subroutine vector_integral_extrapolation_real

  subroutine scalar_integral_extrapolation_complex(array, sizes, int_bounds, result, info)

    complex(kind=dp), intent(in)  :: array(:), int_bounds(:)
    integer, intent(in)  :: sizes(:)

    complex(kind=dp), intent(out) :: result
    integer, intent(out) :: info

    integer                       :: n1, n2, n3, ep
    real(kind=dp)                 :: nr1, nr2, nr3
    complex(kind=dp), allocatable :: e3(:, :, :), &
                                     e2(:, :)

    result = cmplx(0.0_dp, 0.0_dp, dp)

    if (2*size(sizes) .ne. size(int_bounds)) then
      info = -1 !Error.
      return
    endif

    if (product(sizes) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    if (size(sizes) .eq. 3) then

      !Approximate to 0.
      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(4) - int_bounds(3)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(6) - int_bounds(5)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1 !Sucess.
      !  return

      !endif

      !Calculate n_i = log_2(N_i - 1)
      !Exceptional case: N_i = 1:
      !Assign an otherwise unachievable negative number and use it
      !as a flag for nscalar1_integrate_complex. That dimension will
      !not be integrated over.
      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif
      if (sizes(2) .eq. 1) then
        nr2 = -1.0_dp
        n2 = -1
      else
        nr2 = log(real(sizes(2) - 1, dp))/log(2.0_dp)
        n2 = nint(nr2)
      endif
      if (sizes(3) .eq. 1) then
        nr3 = -1.0_dp
        n3 = -1
      else
        nr3 = log(real(sizes(3) - 1, dp))/log(2.0_dp)
        n3 = nint(nr3)
      endif

      !Check if the mesh is suitable for extrapolation.
      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp) .or. &
          (abs(real(n2, dp) - nr2) .ge. 1.0E-6_dp) .or. &
          (abs(real(n3, dp) - nr3) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))* &
                             (int_bounds(4) - int_bounds(3))* &
                             (int_bounds(6) - int_bounds(5))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        allocate (e3(sizes(1), sizes(2), sizes(3)))
        !i) Pass from memory to arbitrary array,
        call expand_carray3(array(:), e3(:, :, :), ep)
        !ii) integrate as if the data was scattered from [0, 1]x[0, 1]x[0, 1],
        result = nscalar3_integrate_complex(e3(:, :, :), n1, n2, n3)
        deallocate (e3)
        !iii) multiply by integral bounds so the integral is properly normalized.
        result = result*((int_bounds(2) - int_bounds(1))* &
                         (int_bounds(4) - int_bounds(3))* &
                         (int_bounds(6) - int_bounds(5)))
        info = 1 !Sucess.

      endif

    elseif (size(sizes) .eq. 2) then

      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(4) - int_bounds(3)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1
      !  return

      !endif

      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif
      if (sizes(2) .eq. 1) then
        nr2 = -1.0_dp
        n2 = -1
      else
        nr2 = log(real(sizes(2) - 1, dp))/log(2.0_dp)
        n2 = nint(nr2)
      endif

      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp) .or. &
          (abs(real(n2, dp) - nr2) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))* &
                             (int_bounds(4) - int_bounds(3))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        allocate (e2(sizes(1), sizes(2)))
        call expand_carray2(array(:), e2(:, :), ep)
        result = nscalar2_integrate_complex(e2(:, :), n1, n2)
        deallocate (e2)
        result = result*((int_bounds(2) - int_bounds(1))* &
                         (int_bounds(4) - int_bounds(3)))
        info = 1

      endif

    elseif (size(sizes) .eq. 1) then

      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1
      !  return

      !endif

      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif

      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        result = nscalar1_integrate_complex(array(:), n1)
        result = result*((int_bounds(2) - int_bounds(1)))
        info = 1

      endif

    else

      info = -1
      return

    endif

  end subroutine scalar_integral_extrapolation_complex

  subroutine vector_integral_extrapolation_complex(array, sizes, int_bounds, result, info)

    complex(kind=dp), intent(in)  :: array(:, :), int_bounds(:)
    integer, intent(in)           :: sizes(:)

    complex(kind=dp), intent(out) :: result(size(array(1, :)))
    integer, intent(out)          :: info

    integer                       :: i

    result = cmplx(0.0_dp, 0.0_dp, dp)

    if (2*size(sizes) .ne. size(int_bounds)) then
      info = -1 !Error.
      return
    endif

    if (product(sizes) .ne. size(array(:, 1))) then
      info = -1 !Error.
      return
    endif

    do i = 1, size(array(1, :))
      call scalar_integral_extrapolation_complex(array(:, i), sizes, int_bounds, result(i), info)
    enddo

  end subroutine vector_integral_extrapolation_complex

  subroutine scalar_integral_extrapolation_complex_array_real_bounds(array, sizes, int_bounds, result, info)

    complex(kind=dp), intent(in)  :: array(:)
    real(kind=dp), intent(in)  :: int_bounds(:)
    integer, intent(in)  :: sizes(:)

    complex(kind=dp), intent(out) :: result
    integer, intent(out) :: info

    integer                       :: n1, n2, n3, ep
    real(kind=dp)                 :: nr1, nr2, nr3
    complex(kind=dp), allocatable :: e3(:, :, :), &
                                     e2(:, :)

    result = cmplx(0.0_dp, 0.0_dp, dp)

    if (2*size(sizes) .ne. size(int_bounds)) then
      info = -1 !Error.
      return
    endif

    if (product(sizes) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    if (size(sizes) .eq. 3) then

      !Approximate to 0.
      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(4) - int_bounds(3)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(6) - int_bounds(5)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1 !Sucess.
      !  return

      !endif

      !Calculate n_i = log_2(N_i - 1)
      !Exceptional case: N_i = 1:
      !Assign an otherwise unachievable negative number and use it
      !as a flag for nscalar1_integrate_complex. That dimension will
      !not be integrated over.
      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif
      if (sizes(2) .eq. 1) then
        nr2 = -1.0_dp
        n2 = -1
      else
        nr2 = log(real(sizes(2) - 1, dp))/log(2.0_dp)
        n2 = nint(nr2)
      endif
      if (sizes(3) .eq. 1) then
        nr3 = -1.0_dp
        n3 = -1
      else
        nr3 = log(real(sizes(3) - 1, dp))/log(2.0_dp)
        n3 = nint(nr3)
      endif

      !Check if the mesh is suitable for extrapolation.
      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp) .or. &
          (abs(real(n2, dp) - nr2) .ge. 1.0E-6_dp) .or. &
          (abs(real(n3, dp) - nr3) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))* &
                             (int_bounds(4) - int_bounds(3))* &
                             (int_bounds(6) - int_bounds(5))) &
                 /size(array)
        !Not suitable, returning the rectangle method approximation.
        return

      else

        allocate (e3(sizes(1), sizes(2), sizes(3)))
        !i) Pass from memory to arbitrary array,
        call expand_carray3(array(:), e3(:, :, :), ep)
        !ii) integrate as if the data was scattered from [0, 1]x[0, 1]x[0, 1],
        result = nscalar3_integrate_complex(e3(:, :, :), n1, n2, n3)
        deallocate (e3)
        !iii) multiply by integral bounds so the integral is properly normalized.
        result = result*((int_bounds(2) - int_bounds(1))* &
                         (int_bounds(4) - int_bounds(3))* &
                         (int_bounds(6) - int_bounds(5)))
        info = 1 !Sucess.

      endif

    elseif (size(sizes) .eq. 2) then

      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp) .or. &
      !    (abs(int_bounds(4) - int_bounds(3)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1
      !  return

      !endif

      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif
      if (sizes(2) .eq. 1) then
        nr2 = -1.0_dp
        n2 = -1
      else
        nr2 = log(real(sizes(2) - 1, dp))/log(2.0_dp)
        n2 = nint(nr2)
      endif

      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp) .or. &
          (abs(real(n2, dp) - nr2) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))* &
                             (int_bounds(4) - int_bounds(3))) &
                 /size(array)
        return

      else

        allocate (e2(sizes(1), sizes(2)))
        call expand_carray2(array(:), e2(:, :), ep)
        result = nscalar2_integrate_complex(e2(:, :), n1, n2)
        deallocate (e2)
        result = result*((int_bounds(2) - int_bounds(1))* &
                         (int_bounds(4) - int_bounds(3)))
        info = 1

      endif

    elseif (size(sizes) .eq. 1) then

      !if ((abs(int_bounds(2) - int_bounds(1)) .le. 1.0E-6_dp)) then

      !  result = cmplx(0.0_dp, 0.0_dp, dp)
      !  info = 1
      !  return

      !endif

      if (sizes(1) .eq. 1) then
        nr1 = -1.0_dp
        n1 = -1
      else
        nr1 = log(real(sizes(1) - 1, dp))/log(2.0_dp)
        n1 = nint(nr1)
      endif

      if ((abs(real(n1, dp) - nr1) .ge. 1.0E-6_dp)) then

        info = 0
        result = sum(array)*((int_bounds(2) - int_bounds(1))) &
                 /size(array)
        return

      else

        result = nscalar1_integrate_complex(array(:), n1)
        result = result*((int_bounds(2) - int_bounds(1)))
        info = 1

      endif

    else

      info = -1
      return

    endif

  end subroutine scalar_integral_extrapolation_complex_array_real_bounds

  subroutine vector_integral_extrapolation_complex_array_real_bounds(array, sizes, int_bounds, result, info)

    complex(kind=dp), intent(in)  :: array(:, :)
    real(kind=dp), intent(in)     :: int_bounds(:)
    integer, intent(in)           :: sizes(:)

    complex(kind=dp), intent(out) :: result(size(array(1, :)))
    integer, intent(out)          :: info

    integer                       :: i

    result = cmplx(0.0_dp, 0.0_dp, dp)

    if (2*size(sizes) .ne. size(int_bounds)) then
      info = -1 !Error.
      return
    endif

    if (product(sizes) .ne. size(array(:, 1))) then
      info = -1 !Error.
      return
    endif

    do i = 1, size(array(1, :))
      call scalar_integral_extrapolation_complex_array_real_bounds(array(:, i), sizes, int_bounds, result(i), info)
    enddo

  end subroutine vector_integral_extrapolation_complex_array_real_bounds

  function nscalar3_integrate_real(array, n1, n2, n3) result(u)

    !Utility to integrate a real f(x) (R^3->R) function over the
    ![0, 1]x[0, 1]x[0, 1] real interval using the extrapolation method.
    !The data is stored in the input array
    !such that array(i, j, k) = f(x_i, y_j, z_k), x_i = (i-1)/(2^n1), i \in [1, 2^n1 + 1],
    !y_i = (i-1)/(2^n2), i \in [1, 2^n2 + 1], z_i = (i-1)/(2^n2), i \in [1, 2^n3 + 1].
    !size(array) must be expressible as (2^n1 + 1)*(2^n2 + 1)*(2^n3 + 1), as such, n1, n2, n3 must be provided
    !on input.

    real(kind=dp), intent(in) :: array(:, :, :)
    integer, intent(in)       :: n1, n2, n3

    real(kind=dp)             :: u
    real(kind=dp)             :: aux(size(array(:, 1, 1)), size(array(1, :, 1)))
    integer                   :: i, j

    do i = 1, size(array(:, 1, 1))
      do j = 1, size(array(1, :, 1))
        aux(i, j) = nscalar1_integrate_real(array(i, j, :), n3)
      enddo
    enddo

    u = nscalar2_integrate_real(aux, n1, n2)

  end function nscalar3_integrate_real

  function nscalar2_integrate_real(array, n1, n2) result(u)

    !Utility to integrate a real f(x) (R^2->R) function over the
    ![0, 1]x[0, 1] real interval using the extrapolation method.
    !The data is stored in the input array
    !such that array(i, j) = f(x_i, y_j), x_i = (i-1)/(2^n1), i \in [1, 2^n1 + 1],
    !y_i = (i-1)/(2^n2), i \in [1, 2^n2 + 1].
    !size(array) must be expressible as (2^n1 + 1)*(2^n2 + 1), as such, n1, n2 must be provided
    !on input.

    real(kind=dp), intent(in) :: array(:, :)
    integer, intent(in)       :: n1, n2

    real(kind=dp)             :: u
    real(kind=dp)             :: aux(size(array(:, 1)))
    integer                   :: i

    do i = 1, size(array(:, 1))
      aux(i) = nscalar1_integrate_real(array(i, :), n2)
    enddo

    u = nscalar1_integrate_real(aux, n1)

  end function nscalar2_integrate_real

  function nscalar1_integrate_real(array, n) result(u)

    !Utility to integrate a real f(x) (R->R) function over the
    ![0, 1] real interval using the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, as such, n must be provided
    !on input.

    real(kind=dp), intent(in) :: array(:)
    integer, intent(in)       :: n

    real(kind=dp)             :: u
    real(kind=dp)             :: ord_array(size(array))

    if (n .eq. 0) then
      !Special case size(array) = 2, use simple trapezium method.
      u = sum(array)*0.5_dp
    else if (n .eq. -1) then
      !Exception size(array) = 1, do not integrate.
      u = array(1)
    else
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = organize_data(array, n)
      u = extrapolation(ord_array, n)
    endif

  contains

    function extrapolation(array, n) result(u)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, u, is the last extrapolated result.

      real(kind=dp), intent(in) :: array(:)
      integer, intent(in) :: n

      real(kind=dp)             :: u
      real(kind=dp)             :: workarray(n)
      integer                   :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trapezium_method(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_dp)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_dp + (4.0_dp)**(i - 1))
        enddo
      enddo

      u = workarray(1)

    end function extrapolation

    function trapezium_method(array) result(u)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally spaced
      !and practically ordered by organize_data.

      real(kind=dp), intent(in) :: array(:)
      real(kind=dp)             :: u

      integer                   :: i

      u = 0.5_dp*(array(1) + array(2))
      do i = 3, size(array)
        u = u + array(i)
      enddo
      u = u/real(size(array) - 1, dp)

    end function trapezium_method

    function organize_data(array, n) result(org_array)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates the quantity of data points
      !by: size(array) = 1 + 2^n.

      real(kind=dp), intent(in) :: array(:)
      integer, intent(in) :: n
      real(kind=dp)             :: org_array(size(array))

      integer                   :: i, j, k, l

      !Manual assignation to endpoints.
      org_array(1) = array(1)
      org_array(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_array(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function organize_data

  end function nscalar1_integrate_real

  function nscalar3_integrate_complex(array, n1, n2, n3) result(u)

    !Utility to integrate a complex f(x) (R^3->R) function over the
    ![0, 1]x[0, 1]x[0, 1] real interval using the extrapolation method.
    !The data is stored in the input array
    !such that array(i, j, k) = f(x_i, y_j, z_k), x_i = (i-1)/(2^n1), i \in [1, 2^n1 + 1],
    !y_i = (i-1)/(2^n2), i \in [1, 2^n2 + 1], z_i = (i-1)/(2^n2), i \in [1, 2^n3 + 1].
    !size(array) must be expressible as (2^n1 + 1)*(2^n2 + 1)*(2^n3 + 1), as such, n1, n2, n3 must be provided
    !on input.

    complex(kind=dp), intent(in) :: array(:, :, :)
    integer, intent(in)          :: n1, n2, n3

    complex(kind=dp)             :: u
    complex(kind=dp)             :: aux(size(array(:, 1, 1)), size(array(1, :, 1)))
    integer                      :: i, j

    do i = 1, size(array(:, 1, 1))
      do j = 1, size(array(1, :, 1))
        aux(i, j) = nscalar1_integrate_complex(array(i, j, :), n3)
      enddo
    enddo

    u = nscalar2_integrate_complex(aux, n1, n2)

  end function nscalar3_integrate_complex

  function nscalar2_integrate_complex(array, n1, n2) result(u)

    !Utility to integrate a complex f(x) (R^2->C) function over the
    ![0, 1]x[0, 1] real interval using the extrapolation method.
    !The data is stored in the input array
    !such that array(i, j) = f(x_i, y_j), x_i = (i-1)/(2^n1), i \in [1, 2^n1 + 1],
    !y_i = (i-1)/(2^n2), i \in [1, 2^n2 + 1].
    !size(array) must be expressible as (2^n1 + 1)*(2^n2 + 1), as such, n1, n2 must be provided
    !on input.

    complex(kind=dp), intent(in) :: array(:, :)
    integer, intent(in)          :: n1, n2

    complex(kind=dp)             :: u
    complex(kind=dp)             :: aux(size(array(:, 1)))
    integer                      :: i

    do i = 1, size(array(:, 1))
      aux(i) = nscalar1_integrate_complex(array(i, :), n2)
    enddo

    u = nscalar1_integrate_complex(aux, n1)

  end function nscalar2_integrate_complex

  function nscalar1_integrate_complex(array, n) result(u)

    !Utility to integrate a complex f(x) (R->C) function over the
    ![0, 1] real interval using the extrapolation method.
    !The data is stored in the input array
    !such that array(i) = f(x_i), x_i = (i-1)/(2^n), i \in [1, 2^n + 1].
    !size(array) must be expressible as 2^n + 1, as such, n must be provided
    !on input.

    complex(kind=dp), intent(in) :: array(:)
    integer, intent(in)          :: n

    complex(kind=dp)             :: u
    complex(kind=dp)             :: ord_array(size(array))

    if (n .eq. 0) then
      !Special case size(array) = 2, use simple trapezium method.
      u = sum(array)*0.5_dp
    else if (n .eq. -1) then
      !Exception size(array) = 1, do not integrate.
      u = array(1)
    else
      !size(array) = 3, 5, 9, ... proceed normally.
      ord_array = organize_data(array, n)
      u = extrapolation(ord_array, n)
    endif

  contains

    function extrapolation(array, n) result(u)

      !Implementation of the extrapolation scheme.
      !Calculates trapezium rule approximation for the data array and different
      !values of the step. Practical order is assumed in the data array.
      !The integer n is the "order" of the data. Relates the quantity of data points
      !by: size(array) = 1 + 2^n.
      !The different I_i(1/2^j) are then calculated using extrapolation.
      !The output, u, is the last extrapolated result.

      complex(kind=dp), intent(in) :: array(:)
      integer, intent(in) :: n

      complex(kind=dp)             :: u
      complex(kind=dp)             :: workarray(n)
      integer                      :: i, j

      do i = 1, n !Trapezium rule. Implemented using practical order.
        workarray(i) = trapezium_method(array(1:2**(i) + 1))
      enddo

      do i = 2, n !Extrapolation scheme. Calculates I_i(1/2^j).
        do j = 1, n - i + 1
          workarray(j) = ((4.0_dp)**(i - 1))*workarray(j + 1) - workarray(j)
          workarray(j) = workarray(j)/(-1.0_dp + (4.0_dp)**(i - 1))
        enddo
      enddo

      u = workarray(1)

    end function extrapolation

    function trapezium_method(array) result(u)

      !Implementation of the composite trapezium method.
      !Assigns a weight of 0.5 to the endpoints and a weight of 1
      !to internal points. Later divides by the number of steps.
      !Data from array is assumed to be in the [0, 1] interval equally spaced
      !and practically ordered by organize_data.

      complex(kind=dp), intent(in) :: array(:)
      complex(kind=dp)             :: u

      integer                      :: i

      u = 0.5_dp*(array(1) + array(2))
      do i = 3, size(array)
        u = u + array(i)
      enddo
      u = u/real(size(array) - 1, dp)

    end function trapezium_method

    function organize_data(array, n) result(org_array)

      !Implementation of the practical ordering map.
      !Passes data array form being canonically ordered to practically ordered.
      !The integer n is the "order" of the data. Relates the quantity of data points
      !by: size(array) = 1 + 2^n.

      complex(kind=dp), intent(in) :: array(:)
      integer, intent(in) :: n
      complex(kind=dp)             :: org_array(size(array))

      integer                      :: i, j, k, l

      !Manual assignation to endpoints.
      org_array(1) = array(1)
      org_array(2) = array(size(array))

      !Intermediate points.
      l = 2
      do i = 2, 1 + n
        j = 2**(n - i + 1) + 1
        do k = 1, 2**(i - 2)
          l = l + 1
          org_array(l) = array(j + 2*(j - 1)*(k - 1))
        enddo
      enddo

    end function organize_data

  end function nscalar1_integrate_complex

  subroutine shrink_rarray1(array, shrink, info)

    !Pass from dim = 1 arbitrary array to memory layout.

    real(kind=dp), intent(in)  :: array(:)

    real(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)       :: info

    integer                    :: a1, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      shrink(l) = array(a1)
      l = l + 1
    enddo

    info = 1 !Sucess.

  end subroutine shrink_rarray1

  subroutine shrink_rarray2(array, shrink, info)

    !Pass from dim = 2 arbitrary array to memory layout.

    real(kind=dp), intent(in)  :: array(:, :)

    real(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)       :: info

    integer                    :: a1, a2, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      do a2 = lbound(array, 2), ubound(array, 2)
        shrink(l) = array(a1, a2)
        l = l + 1
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine shrink_rarray2

  subroutine shrink_rarray3(array, shrink, info)

    !Pass from dim = 3 arbitrary array to memory layout.

    real(kind=dp), intent(in)  :: array(:, :, :)

    real(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)       :: info

    integer                    :: a1, a2, a3, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      do a2 = lbound(array, 2), ubound(array, 2)
        do a3 = lbound(array, 3), ubound(array, 3)
          shrink(l) = array(a1, a2, a3)
          l = l + 1
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine shrink_rarray3

  subroutine shrink_rarray4(array, shrink, info)

    !Pass from dim = 4 arbitrary array to memory layout.

    real(kind=dp), intent(in)  :: array(:, :, :, :)

    real(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)       :: info

    integer                    :: a1, a2, a3, a4, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      do a2 = lbound(array, 2), ubound(array, 2)
        do a3 = lbound(array, 3), ubound(array, 3)
          do a4 = lbound(array, 4), ubound(array, 4)
            shrink(l) = array(a1, a2, a3, a4)
            l = l + 1
          enddo
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine shrink_rarray4

  subroutine shrink_carray1(array, shrink, info)

    !Pass from dim = 1 arbitrary array to memory layout.

    complex(kind=dp), intent(in)  :: array(:)

    complex(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)          :: info

    integer                       :: a1, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      shrink(l) = array(a1)
      l = l + 1
    enddo

    info = 1 !Sucess.

  end subroutine shrink_carray1

  subroutine shrink_carray2(array, shrink, info)

    !Pass from dim = 2 arbitrary array to memory layout.

    complex(kind=dp), intent(in)  :: array(:, :)

    complex(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)          :: info

    integer                       :: a1, a2, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      do a2 = lbound(array, 2), ubound(array, 2)
        shrink(l) = array(a1, a2)
        l = l + 1
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine shrink_carray2

  subroutine shrink_carray3(array, shrink, info)

    !Pass from dim = 3 arbitrary array to memory layout.

    complex(kind=dp), intent(in)  :: array(:, :, :)

    complex(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)          :: info

    integer                       :: a1, a2, a3, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      do a2 = lbound(array, 2), ubound(array, 2)
        do a3 = lbound(array, 3), ubound(array, 3)
          shrink(l) = array(a1, a2, a3)
          l = l + 1
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine shrink_carray3

  subroutine shrink_carray4(array, shrink, info)

    !Pass from dim = 4 arbitrary array to memory layout.

    complex(kind=dp), intent(in)  :: array(:, :, :, :)

    complex(kind=dp), intent(out) :: shrink(:)
    integer, intent(out)          :: info

    integer                       :: a1, a2, a3, a4, l

    shrink = 0.0_dp
    info = 0

    if (size(shrink) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(array, 1), ubound(array, 1)
      do a2 = lbound(array, 2), ubound(array, 2)
        do a3 = lbound(array, 3), ubound(array, 3)
          do a4 = lbound(array, 4), ubound(array, 4)
            shrink(l) = array(a1, a2, a3, a4)
            l = l + 1
          enddo
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine shrink_carray4

  subroutine expand_rarray1(array, expand, info)

    !Pass from memory layout to dim = 1 arbitrary array.

    real(kind=dp), intent(in)  :: array(:)

    real(kind=dp), intent(out) :: expand(:)
    integer, intent(out)       :: info

    integer                    :: a1, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      expand(a1) = array(l)
      l = l + 1
    enddo

    info = 1 !Sucess.

  end subroutine expand_rarray1

  subroutine expand_rarray2(array, expand, info)

    !Pass from memory layout to dim = 2 arbitrary array.

    real(kind=dp), intent(in)  :: array(:)

    real(kind=dp), intent(out) :: expand(:, :)
    integer, intent(out)       :: info

    integer                    :: a1, a2, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      do a2 = lbound(expand, 2), ubound(expand, 2)
        expand(a1, a2) = array(l)
        l = l + 1
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine expand_rarray2

  subroutine expand_rarray3(array, expand, info)

    !Pass from memory layout to dim = 3 arbitrary array.

    real(kind=dp), intent(in)  :: array(:)

    real(kind=dp), intent(out) :: expand(:, :, :)
    integer, intent(out)       :: info

    integer                    :: a1, a2, a3, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      do a2 = lbound(expand, 2), ubound(expand, 2)
        do a3 = lbound(expand, 3), ubound(expand, 3)
          expand(a1, a2, a3) = array(l)
          l = l + 1
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine expand_rarray3

  subroutine expand_rarray4(array, expand, info)

    !Pass from memory layout to dim = 4 arbitrary array.

    real(kind=dp), intent(in)  :: array(:)

    real(kind=dp), intent(out) :: expand(:, :, :, :)
    integer, intent(out)       :: info

    integer                    :: a1, a2, a3, a4, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      do a2 = lbound(expand, 2), ubound(expand, 2)
        do a3 = lbound(expand, 3), ubound(expand, 3)
          do a4 = lbound(expand, 4), ubound(expand, 4)
            expand(a1, a2, a3, a4) = array(l)
            l = l + 1
          enddo
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine expand_rarray4

  subroutine expand_carray1(array, expand, info)

    !Pass from memory layout to dim = 1 arbitrary array.

    complex(kind=dp), intent(in)  :: array(:)

    complex(kind=dp), intent(out) :: expand(:)
    integer, intent(out)          :: info

    integer                       :: a1, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      expand(a1) = array(l)
      l = l + 1
    enddo

    info = 1 !Sucess.

  end subroutine expand_carray1

  subroutine expand_carray2(array, expand, info)

    !Pass from memory layout to dim = 2 arbitrary array.

    complex(kind=dp), intent(in)  ::  array(:)

    complex(kind=dp), intent(out) :: expand(:, :)
    integer, intent(out)          :: info

    integer                       :: a1, a2, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      do a2 = lbound(expand, 2), ubound(expand, 2)
        expand(a1, a2) = array(l)
        l = l + 1
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine expand_carray2

  subroutine expand_carray3(array, expand, info)

    !Pass from memory layout to dim = 3 arbitrary array.

    complex(kind=dp), intent(in)  :: array(:)

    complex(kind=dp), intent(out) :: expand(:, :, :)
    integer, intent(out)          :: info

    integer                       :: a1, a2, a3, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      do a2 = lbound(expand, 2), ubound(expand, 2)
        do a3 = lbound(expand, 3), ubound(expand, 3)
          expand(a1, a2, a3) = array(l)
          l = l + 1
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine expand_carray3

  subroutine expand_carray4(array, expand, info)

    !Pass from memory layout to dim = 4 arbitrary array.

    complex(kind=dp), intent(in)  :: array(:)

    complex(kind=dp), intent(out) :: expand(:, :, :, :)
    integer, intent(out)          :: info

    integer                       :: a1, a2, a3, a4, l

    expand = 0.0_dp
    info = 0

    if (size(expand) .ne. size(array)) then
      info = -1 !Error.
      return
    endif

    l = 1
    do a1 = lbound(expand, 1), ubound(expand, 1)
      do a2 = lbound(expand, 2), ubound(expand, 2)
        do a3 = lbound(expand, 3), ubound(expand, 3)
          do a4 = lbound(expand, 4), ubound(expand, 4)
            expand(a1, a2, a3, a4) = array(l)
            l = l + 1
          enddo
        enddo
      enddo
    enddo

    info = 1 !Sucess.

  end subroutine expand_carray4

end module extrapolation_integration
