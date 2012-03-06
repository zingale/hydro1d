FC = gfortran
FFLAGS = -c -O2 -g -C

all: hydro1d

FSOURCE = main.f90 grid.f90 datatypes.f90 params.f90 \
          init.f90 variables.f90 eos.f90 output.f90 \
          bcs.f90 dt.f90 states.f90 riemann.f90 update.f90

# dependencies
bcs.o       : params.o variables.o grid.o datatypes.o 
datatypes.o : 
dt.o        : variables.o eos.o params.o grid.o datatypes.o 
eos.o       : params.o datatypes.o 
grid.o      : datatypes.o 
init.o      : variables.o eos.o grid.o params.o datatypes.o 
states.o    : variables.o grid.o datatypes.o params.o
output.o    : eos.o variables.o params.o datatypes.o grid.o 
params.o    : datatypes.o 
riemann.o   : variables.o grid.o datatypes.o params.o eos.o
update.o    : variables.o grid.o datatypes.o 
variables.o : 

main.o      : update.o riemann.o states.o dt.o bcs.o output.o \
              variables.o init.o params.o grid.o datatypes.o 


%.o: %.f90
	$(FC) $(FFLAGS) $< 

OBJECTS = $(FSOURCE:.f90=.o)

hydro1d: $(OBJECTS)
	gfortran -o hydro1d $(OBJECTS)

.PHONY: clean

clean:
	rm -f *.o *.mod
