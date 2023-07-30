default: main

F90 = ifort
F90FLAGS = -fPIE -g -warn all -check bounds -O2 -qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread

SRC = ./src
OBJ = ./src/obj
BIN = ./bin

CALC = ./src/calculators

utility.o: $(SRC)/utility.F90
					 $(F90) $(F90FLAGS) -c $(SRC)/utility.F90 -o "$(OBJ)/utility.o"

data_structures.o: $(SRC)/data_structures.F90 utility.o
									 $(F90) $(F90FLAGS) -c $(SRC)/data_structures.F90 -o "$(OBJ)/data_structures.o"

extrapolation_integration.o: $(SRC)/extrapolation_integration.F90 utility.o
					 									 $(F90) $(F90FLAGS) -c $(SRC)/extrapolation_integration.F90 -o "$(OBJ)/extrapolation_integration.o"


local_k_quantities.o: $(SRC)/local_k_quantities.F90 utility.o data_structures.o
											$(F90) $(F90FLAGS) -c $(SRC)/local_k_quantities.F90 -o "$(OBJ)/local_k_quantities.o"

kpath.o: $(SRC)/kpath.F90 utility.o data_structures.o
				 $(F90) $(F90FLAGS) -c $(SRC)/kpath.F90 -o "$(OBJ)/kpath.o"

kslice.o: $(SRC)/kslice.F90 utility.o data_structures.o
				  $(F90) $(F90FLAGS) -c $(SRC)/kslice.F90 -o "$(OBJ)/kslice.o"

integrator.o : $(SRC)/integrator.F90 utility.o extrapolation_integration.o data_structures.o
							 $(F90) $(F90FLAGS) -c $(SRC)/integrator.F90 -o "$(OBJ)/integrator.o"

calculators_general.o: $(CALC)/calculators_general.F90 utility.o data_structures.o local_k_quantities.o kpath.o
									 $(F90) $(F90FLAGS) -c $(CALC)/calculators_general.F90 -o "$(OBJ)/calculators_general.o"

calculators_floquet.o: $(CALC)/calculators_floquet.F90 utility.o data_structures.o local_k_quantities.o kpath.o
									 $(F90) $(F90FLAGS) -c $(CALC)/calculators_floquet.F90 -o "$(OBJ)/calculators_floquet.o"

calculators_optical.o: $(CALC)/calculators_optical.F90 utility.o integrator.o data_structures.o local_k_quantities.o
									 $(F90) $(F90FLAGS) -c $(CALC)/calculators_optical.F90 -o "$(OBJ)/calculators_optical.o"


SsTC.o: $(SRC)/SsTC.F90 utility.o data_structures.o extrapolation_integration.o local_k_quantities.o kpath.o kslice.o integrator.o calculators_general.o calculators_floquet.o calculators_optical.o
			$(F90) $(F90FLAGS) -c $(SRC)/SsTC.F90 -o "$(OBJ)/SsTC.o"

main: SsTC.o
			ar cr "$(BIN)/libSsTC.a" $(OBJ)/*.o
			mv *.mod $(BIN)

.PHONY: uninstall
uninstall:
	rm -rf $(OBJ)/*.o *.mod $(BIN)/*.x
