FC = gfortran
FFLAGS = -c -O2 -g -C

all: hydro1d

FSOURCE = main.f90 grid.f90 datatypes.f90 params.f90 \
          init.f90 variables.f90 eos.f90 output.f90 \
          bcs.f90 dt.f90 interface_states.f90 riemann.f90 update.f90

# dependencies
bcs.o : bcs.f90 params.o variables.o grid.o datatypes.o 
datatypes.o : datatypes.f90 
dt.o : dt.f90 variables.o eos.o params.o grid.o datatypes.o 
eos.o : eos.f90 params.o datatypes.o 
grid.o : grid.f90 datatypes.o 
init.o : init.f90 variables.o eos.o grid.o params.o datatypes.o 
interface_states.o : interface_states.f90 variables.o grid.o datatypes.o params.o
main.o : main.f90 update.o riemann.o interface_states.o dt.o bcs.o output.o variables.o init.o params.o grid.o datatypes.o 
output.o : output.f90 eos.o variables.o params.o datatypes.o grid.o 
params.o : params.f90 datatypes.o 
riemann.o : riemann.f90 variables.o grid.o datatypes.o params.o eos.o
update.o : update.f90 variables.o grid.o datatypes.o 
variables.o : variables.f90


%.o: %.f90
	$(FC) $(FFLAGS) $< 

OBJECTS = $(FSOURCE:.f90=.o)

hydro1d: $(OBJECTS)
	gfortran -o hydro1d $(OBJECTS)

.PHONY: clean

clean:
	rm -f *.o *.mod
