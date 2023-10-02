# Solid state Task Constructor - SsTC

## A high-perfomance computing oriented library to create integration and sampling tasks in the BZ of a crystal for k-dependent functions.

# Prerequisites

### Fortran compiler:

[Intel Fortran oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) compilers `mpiifort` (recommended) or `mpiifx`.

### Make software.

### Python 3 (> v3.9):

We recommend Python 3.9.7 :: Intel Corporation.

### Libraries:

Intel oneAPI Math Kernel Libraries (MKL), OpenMP library and MPI library.

Python3's 're' and 'glob' libraries.

# Installation:

1. Clone this repository with the flag `--recurse-submodules` on a destination of your choice.

       bash:/path/of/your/choice$ git clone --recurse-submodules https://github.com/irukoa/SsTC.git

2. Change directory to SsTC.

       bash:/path/of/your/choice$ cd SsTC/

3. Run make. This will create the (static) library file `./bin/libSsTC.a` and the module file `./bin/sstc.mod`.

       bash:/path/of/your/choice/SsTC$ make

4. To uninstall run `make uninstall` in the installation directory.

# Linking to your application:

1. Include the line `use SsTC` in your application preamble.

2. Before creating any SsTC task, sampling/integrating in the BZ or printing to files, make sure to have the MPI enviroment initialized and the line `call SsTC_init()` in your application.

Note 1: SsTC uses double precision numbers for real and complex kinds.

Note 2: It is recommended that each task is defined within a [BLOCK](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/block.html)
 construct to help in derived-type finalization and thus prevent memory leaks.

For example, an application calculating the jerk current of the system GaAs, as in Example 2 os the User's guide, should look like:

    bash:/path/to/application/$ cat my_jerk_application.F90
<!-- tsk -->
    program my_jerk_application

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

    end program my_jerk_application

3. To link SsTC to your program, the compilation command should have the form:

       bash:/path/to/application/$ $(F90) $(F90FLAGS) my_jerk_application.F90 -I/path/of/your/choice/SsTC/bin /path/of/your/choice/SsTC/bin/libSsTC.a -o "my_jdos_application.x"

   Where `$(F90) = mpiifort/mpiifx`, and `$(F90FLAGS)` should include, at least, `-qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -pthread`.

# Usage

See user's manual in the documentation folder.
