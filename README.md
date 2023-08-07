# (S)olid (s)tate (T)ask (C)onstructor

## A high perfomance library to create integration and sampling tasks in the BZ of a system for given k-dependent calculators.

# Prerequisites

## Make software.

## Fortran compiler:

Fortran 2018 (ISO/IEC 1539:2018) complying compiler.

We recommend ifort (IFORT), at least, `--version` 2021.5.0 20211109.

-pthread compile flag.

## Libraries:

Intel's Math Kernel Library (MLK).

OpenMP Library (OMP).

# Installation:

1. Clone this repository on a destination of your choice.

        bash:/path/of/your/choice$ git clone https://github.com/irukoa/SsTC.git

2. Change directory to SsTC.

        bash:/path/of/your/choice$ cd SsTC/

3. Run make. This will create the (static) library file `./bin/libSsTC.a` and the module file `./bin/sstc.mod`.

        bash:/path/of/your/choice/SsTC$ make

4. To uninstall run `make uninstall` in the installation directory.

# Linking to your application:

1. Include the line `use SsTC` in your application preamble.

2. Before creating any SsTC task, sampling/integrating in the BZ or printing to files, make sure to have the line `call SsTC_init()` in your application.

Note 1: By default SsTC uses double precision `kind=8` numbers except for integers, which are the default `kind=4`.

Note 2: It is recommended that each task is defined within a [BLOCK](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/block.html)
 construct to help in derived-type finalization and thus prevent memory leaks.

For example, an application calculating the joint density of states (JDOS) of the system GaAs should look like:

       bash:/path/to/application/$ cat my_jdos_application.F90
<!-- tsk -->
       program my_jdos_application

         use SsTC

         integer, parameter :: dp = 8

         type(SsTC_sys) :: GaAs

         call SsTC_init()

         GaAs = SsTC_sys_constructor("GaAs", "./", efermi = 7.7414_dp)

         JDOS: block

           type(SsTC_optical_BZ_integral_task) :: jdostsk

           call SsTC_default_jdos_constructor(optical_task = jdostsk, &
                                              method = "extrapolation", samples = (/65, 65, 65/), &
                                              omegastart = 0.0_dp, omegaend = 10.0_dp, omegasteps = 100)

           call SsTC_sample_and_integrate_BZ_integral_task(task = jdostsk, &
                                                           system = GaAs)

           call SsTC_print_BZ_integral_task(task = jdostsk, &
                                            system = GaAs)

         end block JDOS

       end program my_jdos_application

3. To link SsTC to your program, the compilation command should have the form:

       bash:/path/to/application/$ $(F90) $(F90FLAGS) my_jdos_application.F90 -I/path/of/your/choice/SsTC/bin /path/of/your/choice/SsTC/bin/libSsTC.a -o "my_jdos_application.x"

   We recommend `$(F90) = ifort`, and `$(F90FLAGS)` should include, at least, `-qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -pthread` or other compiler-specific analogous flags.

# Usage

See user's manual in the documentation.
