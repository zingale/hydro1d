FC = gfortran
FFLAGS = -c -O2

all: hydro1d

FSOURCE = main.f90 grid.f90 datatypes.f90 params.f90 \
          init.f90 variables.f90 eos.f90 output.f90 \
          bcs.f90 dt.f90

# dependencies
grid.o: datatypes.o
main.o: datatypes.o grid.o params.o variables.o init.o output.o bcs.o dt.o
params.o: datatypes.o
init.o: datatypes.o params.o grid.o eos.o variables.o
eos.o: datatypes.o params.o
output.o: datatypes.o grid.o params.o variables.o eos.o
bcs.o: datatypes.o grid.o variables.o params.o
dt.o: datatypes.o grid.o params.o eos.o variables.o

%.o: %.f90
	$(FC) $(FFLAGS) $< 

OBJECTS = $(FSOURCE:.f90=.o)

hydro1d: $(OBJECTS)
	gfortran -o hydro1d $(OBJECTS)

.PHONY: clean

clean:
	rm -f *.o *.mod
