#!/bin/bash
cp Makefile ../Makefile      #Load testing Makefile and build.
(cd ../ && make)
#Compile.
mpif90 testing_driver.F90 test-drive/src/testdrive.F90 suites/suite_l_k_quantities.F90 \
                                                       -fprofile-arcs -ftest-coverage \
                                                       -g -fbacktrace --warn-extra --warn-unused --warn-all \
                                                       -Wno-maybe-uninitialized --check=bounds -fstack-check -fimplicit-none \
                                                       -O0 -fopenmp \
                                                       -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread \
                                                       -I../bin/ ../bin/libSsTC.a \
                                                       -o "tests.x"
