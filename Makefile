default: main

F90 = ifort
F90FLAGS = -g -warn all -O2 -qopenmp -llapack -lblas -lfftw3

SRC = ./src
OBJ = ./src/obj
BIN = ./bin

utility.o: $(SRC)/utility.F90
					 $(F90) $(F90FLAGS) -c $(SRC)/utility.F90 -o "$(OBJ)/utility.o"

calculator.o: $(SRC)/calculator.F90 utility.o
							$(F90) $(F90FLAGS) -c $(SRC)/calculator.F90 $(OBJ)/utility.o -o "$(OBJ)/calculator.o"

extrapolation_integration.o: $(SRC)/extrapolation_integration.F90 utility.o
					 $(F90) $(F90FLAGS) -c $(SRC)/extrapolation_integration.F90 $(OBJ)/utility.o -o "$(OBJ)/extrapolation_integration.o"

integrator.o : $(SRC)/integrator.F90 utility.o extrapolation_integration.o calculator.o
							 $(F90) $(F90FLAGS) -c $(SRC)/integrator.F90 $(OBJ)/utility.o $(OBJ)/extrapolation_integration.o $(OBJ)/calculator.o -o "$(OBJ)/integrator.o"

main: $(SRC)/main.F90 utility.o extrapolation_integration.o calculator.o integrator.o
			$(F90) $(F90FLAGS) $(SRC)/main.F90 $(OBJ)/utility.o $(OBJ)/extrapolation_integration.o $(OBJ)/calculator.o $(OBJ)/integrator.o -o "$(BIN)/tb.x"
			rm *.mod

.PHONY: clean
clean: 
	rm -rf $(OBJ)/*.o *.mod

.PHONY: uninstall
uninstall: 
	rm -rf $(OBJ)/*.o *.mod $(BIN)/*.x