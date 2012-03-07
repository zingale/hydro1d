all: hydro1d

FSOURCE += main.f90 grid.f90 datatypes.f90 params.f90 \
           init.f90 variables.f90 eos.f90 output.f90 \
           bcs.f90 dt.f90 godunov.f90 plm.f90 ppm.f90 \
	   riemann.f90 update.f90 probparams.f90


vpath %.f90 . ..


odir := _build

# dependencies
include ../Ghydro.dep


FFLAGS := -c -O2 -g -C
FFLAGS += -J $(odir) -I $(odir)


$(odir)/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	$(FC) $(FFLAGS) -o $@ $< 


OBJECTS = $(addprefix $(odir)/, $(FSOURCE:.f90=.o))

hydro1d: $(OBJECTS)
	gfortran -o hydro1d $(OBJECTS)

.PHONY: clean

clean:
	rm -f $(odir)/*.o $(odir)/*.mod

realclean: clean
	rmdir $(odir)
	rm -f hydro1d
