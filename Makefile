default: main

F90 = ifort
F90FLAGS = -g -warn all -check bounds -O2 -qopenmp -llapack -lblas -lfftw3

SRC = ./src
OBJ = ./src/obj
BIN = ./bin

utility.o: $(SRC)/utility.F90
					 $(F90) $(F90FLAGS) -c $(SRC)/utility.F90 -o "$(OBJ)/utility.o"

data_structures.o: $(SRC)/data_structures.F90
									 $(F90) $(F90FLAGS) -c $(SRC)/data_structures.F90 -o "$(OBJ)/data_structures.o"

extrapolation_integration.o: $(SRC)/extrapolation_integration.F90
					 									 $(F90) $(F90FLAGS) -c $(SRC)/extrapolation_integration.F90 -o "$(OBJ)/extrapolation_integration.o"

integrator.o : $(SRC)/integrator.F90 extrapolation_integration.o data_structures.o
							 $(F90) $(F90FLAGS) -c $(SRC)/integrator.F90 -o "$(OBJ)/integrator.o"

main: $(SRC)/main.F90 utility.o extrapolation_integration.o integrator.o data_structures.o
			$(F90) $(F90FLAGS) $(SRC)/*.F90 -o "$(BIN)/tb.x"
			rm *.mod

.PHONY: clean
clean: 
	rm -rf $(OBJ)/*.o *.mod

.PHONY: uninstall
uninstall: 
	rm -rf $(OBJ)/*.o *.mod $(BIN)/*.x