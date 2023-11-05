mv ../Makefile temp_Makefile #Save previous Makefile.
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
#Serial test.
mpirun -np 1 ./tests.x
#Retrieve coverage.
(cd ../ && gcov -pb ./src/obj/*.o)
#gcov -pb ../src/obj/*.o
lcov --gcov-tool gcov --capture --directory . --output-file coverage.info
genhtml --output-directory html coverage.info
#Restore previous Makefile and rebuild.
mv temp_Makefile ../Makefile
(cd ../ && make)
