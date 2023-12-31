default: main

F90 = mpif90 -L/usr/lib -llapack -L/usr/lib -lblas
F90FLAGS = -fPIE -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -cpp -std=f2018
F90FLAGS += -g -fbacktrace --warn-extra --warn-unused --warn-all -Wno-maybe-uninitialized --check=bounds -fstack-check -fimplicit-none #Debug flags.
F90FLAGS += -fprofile-arcs -ftest-coverage
F90FLAGS += -O0 #Optimization flags.

PY = python3

SRC = ./src
OBJ = ./src/obj
BIN = ./bin
CALC = ./src/calculators

#Base deps of SsTC.o, these get updated if mods are detected.
DEPS = utility.o extrapolation_integration.o data_structures.o kpath.o kslice.o sampler.o integrator.o local_k_quantities.o comms.o

#Version
MAYOR = 0
MINOR = 3
REVISION = 1
F90FLAGS += -D_VERSION="'$(MAYOR).$(MINOR).$(REVISION)'"

include ./src/calculators/Makefile #Checks for mods.

utility.o: $(SRC)/utility.F90
					 $(F90) $(F90FLAGS) -c $(SRC)/utility.F90 -o "$(OBJ)/utility.o"

comms.o: $(SRC)/comms.F90
				 $(F90) $(F90FLAGS) -c $(SRC)/comms.F90 -o "$(OBJ)/comms.o"

extrapolation_integration.o: $(SRC)/F90-Extrapolation-Integration/integration.F90
														 $(F90) $(F90FLAGS) -c $(SRC)/F90-Extrapolation-Integration/integration.F90 -o "$(OBJ)/extrapolation_integration.o"


data_structures.o: $(SRC)/data_structures.F90 utility.o comms.o
									 $(F90) $(F90FLAGS) -c $(SRC)/data_structures.F90 -o "$(OBJ)/data_structures.o"


kpath.o: $(SRC)/kpath.F90 utility.o comms.o data_structures.o
				 $(F90) $(F90FLAGS) -c $(SRC)/kpath.F90 -o "$(OBJ)/kpath.o"

kslice.o: $(SRC)/kslice.F90 utility.o comms.o extrapolation_integration.o data_structures.o
				  $(F90) $(F90FLAGS) -c $(SRC)/kslice.F90 -o "$(OBJ)/kslice.o"

sampler.o: $(SRC)/sampler.F90 utility.o comms.o extrapolation_integration.o data_structures.o
				  $(F90) $(F90FLAGS) -c $(SRC)/sampler.F90 -o "$(OBJ)/sampler.o"

integrator.o : $(SRC)/integrator.F90 utility.o comms.o extrapolation_integration.o data_structures.o
							 $(F90) $(F90FLAGS) -c $(SRC)/integrator.F90 -o "$(OBJ)/integrator.o"


local_k_quantities.o: $(SRC)/local_k_quantities.F90 utility.o comms.o data_structures.o
											$(F90) $(F90FLAGS) -c $(SRC)/local_k_quantities.F90 -o "$(OBJ)/local_k_quantities.o"


SsTC.o: $(SRC)/SsTC.F90 $(DEPS)
				$(PY) $(SRC)/mod_setup.py
				$(F90) $(F90FLAGS) -c $(SRC)/SsTC_mod.F90 -o "$(OBJ)/SsTC.o"

main: SsTC.o
			ar cr "$(BIN)/libSsTC.a" $(OBJ)/*.o
			mv *.mod $(BIN)
			rm -f $(SRC)/*.mod $(CALC)/*.mod

.PHONY: uninstall
uninstall:
	rm -rf $(OBJ)/*.o *.mod $(BIN)/*
