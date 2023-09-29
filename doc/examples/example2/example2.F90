program example02

  USE OMP_LIB
  USE MPI_F08

  use SsTC

  implicit none

  integer, parameter :: dp = 8

  integer :: ierror

  type(SsTC_sys) :: GaAs

  call MPI_INIT(ierror)

  call SsTC_init()

  GaAs = SsTC_sys_constructor("GaAs", "./", efermi = 7.7414_dp)

  block
  
    type(optical_BZ_integral_task) :: jerk

    call jerk_current_constructor(optical_task = jerk, method = "rectangle", samples = (/100, 100, 100/), &
                                  omegastart = 0.0_dp, omegaend = 10.0_dp, omegasteps = 100)

    call SsTC_sample_and_integrate_BZ_integral_task(task = jerk, &
                                                    system = GaAs)

    call SsTC_print_BZ_integral_task(task = jerk, &
                                     system = GaAs)

  end block

  call MPI_FINALIZE(ierror)

end program example02