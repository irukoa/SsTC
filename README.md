[![Language](https://img.shields.io/badge/-Fortran-734f96?logo=fortran&logoColor=white)](https://github.com/topics/fortran)
[![DOI](https://zenodo.org/badge/659820914.svg)](https://zenodo.org/badge/latestdoi/659820914)
[![codecov](https://codecov.io/gh/irukoa/SsTC/graph/badge.svg?token=OJEVB2LAG9)](https://codecov.io/gh/irukoa/SsTC)
[![Testing suite.](https://github.com/irukoa/SsTC/actions/workflows/CI.yml/badge.svg)](https://github.com/irukoa/SsTC/actions/workflows/CI.yml)

# Solid state Task Constructor - SsTC

A library to automate integration and sampling tasks of quantities in the Brillouin zone of a crystal. Oriented to high-perfomance computing.

This is a work in progress!

# Prerequisites

### Fortran compiler:

We strongly recommend the compilers contained in the [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html) toolkit `ifort` or `ifx`. The [GNU Fortran](https://gcc.gnu.org/wiki/GFortran) compiler `gfortran` is also an option. As for this release, intense testing has been conducted only for these compilers. Future plans involve testing and creating builds for more compilers, with special emphasis on the Nvidia HPC compiler.

### Python 3

### Libraries:

-Intel oneAPI Math Kernel Libraries (MKL), contained in the Intel oneAPI toolkit, or any other implementation of BLAS and LAPACK libraries.

-Any implementation of OpenMP and MPI libraries, we recommend those of the Intel oneAPI toolkit.

Python3's 're' and 'glob' libraries.

# Installation:

1. Clone this repository with the flag `--recurse-submodules` on a destination of your choice.

       bash:/path/of/your/choice$ git clone --recurse-submodules https://github.com/irukoa/SsTC.git

2. Change directory to SsTC.

       bash:/path/of/your/choice$ cd SsTC/

3. Run make. This will create the (static) library file `./bin/libSsTC.a` and the module file `./bin/sstc.mod`.

       bash:/path/of/your/choice/SsTC$ make

   The default build is done with the `ifort` compiler. The `config` folder contains several Makefiles for different compilers. The user can copy the Makefile of his/hers choice into the root SsTC directoty and run `make`.

5. To uninstall run `make uninstall` in the installation directory.

# Working principle

The idea behind SsTC is to automate the most common tasks involved in solid state physics. These include sampling (and sometimes integration) of functions with the form $C^{\alpha}(k;\beta)$, for some or all $k$ in the Brillouin zone (BZ) for a given crystalline system. The indices $\alpha$ represent "integer"-like variables while $\beta$ represents "continuous" (real) variables.

Sampling/integration tasks come in different flavours and are specified as Fortran derived types in the implementation level. These include:

1. type(SsTC_BZ_integral_task). Holds the task to integrate $C^{\alpha}(k;\beta)$ in the BZ.
2. type(SsTC_kpath_task). Holds the task to sample $C^{\alpha}(k;\beta)$ along a predefined path in the BZ.
3. type(SsTC_kslice_task). Holds the task to sample $C^{\alpha}(k;\beta)$ along a predefined slice in the BZ.
4. type(SsTC_sampling_task). Holds the task to sample $C^{\alpha}(k;\beta)$ in a regular grid, or along a predefined set of points in the BZ.

To create tasks one needs to employ the task constructors `SsTC_BZ_integral_task_constructor`, `SsTC_kpath_constructor`, `SsTC_kslice_task_constructor` or `SsTC_sampling_task_constructor`, respectively. While the constructors take different arguments related to the particular task being created (a k-path descriptor in the case of a kpath task... see the User's guide), all four take as inputs some common descriptors of $C^{\alpha}(k;\beta)$. These are:

1. The number of $\alpha$ indices `N_int_ind`. An integer greater than 0.
2. The number of possible values of each of the $\alpha_i$ indices. An integer array with elements greater than 0 and of size `N_int_ind`.
3. The number of $\beta$ indices `N_ext_vars`. An integer greater than 0.
4. The starting and ending value of each of the $\beta_i$ indices and the number of steps between these two values. The first two are double precision arrays of size `N_ext_vars` and the last is an integer array with elements greater than 0 and of size `N_ext_vars`.
5. A procedure pointer `g_calculator` to the user's implementation of $C^{\alpha}(k;\beta)$. This must be a function with interface,

       function SsTC_global_calculator(task, system, k, error) result(u)
         class(SsTC_global_k_data), intent(in) :: task
         type(SsTC_sys), intent(in)            :: system
         real(kind=dp), intent(in)             :: k(3)
         logical, intent(inout)                :: error

         complex(kind=dp) :: u(product(task%integer_indices), product(task%continuous_indices))
       end function SsTC_global_calculator

   so the output is an array storing $C^{\alpha}(k;\beta)$ for each combination of the values $\alpha_i$ and $\beta_i$. It is suggested to iterate through all these possible combinations by using the SsTC jagged array indexing functions `SsTC_integer_memory_element_to_array_element` and `SsTC_continuous_memory_element_to_array_element` as

       do i = 1, product(task%integer_indices)
         i_arr = SsTC_integer_memory_element_to_array_element(task, i)
         !i_arr(j) contains the particular value $\alpha_j$ has in this iteration.
         do r = 1, product(task%continuous_indices)
           r_arr = SsTC_continuous_memory_element_to_array_element(task, r)
           !r_arr(j) contains the particular iteration of the variable $\beta_j$,
           !with value task%ext_var_data(j)%data(r_arr(j)).
           u(i, r) = ...
         enddo
       enddo

The user can also pass as an input a cristalline system, which corresponds to the tight-binding representation of said system. This is a file given in the [Wannier90](https://wannier.org/) format, containing the Hamiltonian matrix elements and the position operator's matrix elements in the Wannier basis. A system is a derived type `type(SsTC_sys)`, which can be loaded from a file by using the `SsTC_sys_constructor`. The user is then encouraged to employ the Wannier interpolation routines contained in the `local_k_quantities.F90` module to compute the Hamiltonian and Berri connection matrices and its derivatives in the Wannier gauge.

Once created a task, this can be sampled/integrated by using the task specific routines `SsTC_sample_and_integrate_BZ_integral_task`, `SsTC_kpath_sampler`, `SsTC_sample_kslice_task` and `SsTC_sample_sampling_task`. Which take as inputs a task of the corresponding type and a system of type `type(SsTC_sys)`. This will write to the particular result variable contained in the type (see User's guide). Lastly, the results can be plotted to files by using the routines `SsTC_print_BZ_integral_task`, `SsTC_print_kpath`, `SsTC_print_kslice`, and `SsTC_print_sampling`. These will write a file for each possible combination of $\alpha_i$ integer indices and each of the columns will correspond to a particular $\beta_i$ index and/or a $k$ value (in the case of sampling).

## Parallelization and concerns

The sampling/integration routines are hybrid MPI+OpenMP routines. The working principle is that for every batch of $k$ points these are distributed first between MPI ranks and then between OpenMP threads. The idea behind this model is that SsTC can run in multinode clusters with the number of MPI ranks being the number of nodes and the number of OpenMP threads the number of threads in each node. The user employing SsTC in his/her program needs to make sure that MPI has been initialized at the start of the application followed by the routine `SsTC_init`, before any SsTC tasks/routines have been employed. For personal computers, we recommend running applications using in the MPI singleton case (a single MPI rank).

Lastly, the user is free to employ SsTC routines inside OpenMP parallel regions or other parallel processes in his/hers code (such as parallelization of some or all external variables $\beta_i$). In this case, it is recommended to turn off the OpenMP parallelization done by SsTC by using the OpenMP's `omp_set_max_active_levels` routine (to avoid unecessary overhead) and to check that there is no interference between the user's intended use of his/hers MPI directives and those of SsTC.

# Linking to your application:

1. Include the line `use SsTC` in your application preamble.

2. Before creating any SsTC task, sampling/integrating in the BZ or writing to files, make sure to have the MPI enviroment initialized and the line `call SsTC_init()` in your application after the MPI initialization call.

3. To link SsTC to your program, the compilation command should have the form:

       bash:/path/to/application/$ $(F90) $(F90FLAGS) my_application.F90 -I/path/of/your/choice/SsTC/bin /path/of/your/choice/SsTC/bin/libSsTC.a -o "my_application.x"

   Where `$(F90) = mpiifort/mpiifx/mpif90`, and `$(F90FLAGS)` should include, at least, `-qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -pthread` or the compiler-specific alternatives.

Note 1: SsTC uses double precision numbers for real and complex kinds.

Note 2: It is recommended that each task is defined within a [BLOCK](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2023-0/block.html)
 construct to help in derived-type finalization.

As an example, consider Example 1 in the User's guide,

    program example01

      USE OMP_LIB
      USE MPI_F08

      use SsTC

      implicit none

      integer, parameter :: dp = 8

      integer :: ierror

      type(SsTC_sys) :: dummy

      call MPI_INIT(ierror)

      call SsTC_init()

      dummy = SsTC_sys_constructor("dummy", "./", efermi=0.0_dp)

      block

        type(SsTC_BZ_integral_task) :: test_integral

        call SsTC_BZ_integral_task_constructor(task=test_integral, name="example01-test", &
                                               g_calculator=test_calculator, &
                                               method="rectangle", samples=(/33, 1, 1/), &
                                               N_int_ind=2, int_ind_range=(/3, 3/), &
                                               N_ext_vars=1, ext_vars_start=(/1.0_dp/), &
                                               ext_vars_end=(/10.0_dp/), ext_vars_steps=(/10/))

        call SsTC_sample_and_integrate_BZ_integral_task(task=test_integral, &
                                                        system=dummy)

        call SsTC_print_BZ_integral_task(task=test_integral, &
                                         system=dummy)

      end block

      call MPI_FINALIZE(ierror)

    contains

      function test_calculator(task, system, k, error) result(u)

        class(SsTC_global_k_data), intent(in) :: task
        type(SsTC_sys), intent(in)            :: system
        real(kind=dp), intent(in)             :: k(3)
        logical, intent(inout)                :: error

        complex(kind=dp) :: u(product(task%integer_indices), product(task%continuous_indices))

        integer :: i, r, & !Indices for memory layout.
                   i_arr(2), r_arr(1) !Indices for array layout.

        u = cmplx(0.0, 0.0, dp)

        do i = 1, product(task%integer_indices)
          i_arr = SsTC_integer_memory_element_to_array_element(task, i)
          do r = 1, product(task%continuous_indices)
            r_arr = SsTC_continuous_memory_element_to_array_element(task, r)
	        !C^{\alpha}(k;\beta) = (\alpha_1 + \alpha_2)exp(k_x*\beta_1).
            u(i, r) = real(i_arr(1) + i_arr(2), dp)* &
                      cmplx(exp(k(1)*task%ext_var_data(1)%data(r_arr(1))))
          enddo
        enddo

      end function test_calculator

    end program example01

As another example, an application calculating the jerk current of the system GaAs, as in Example 2 of the User's guide, should look like:

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
