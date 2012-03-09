# basic rules for building the source.  Problems are built in
# sub-directories.  Some bits taken from BoxLib (CCSE/LBNL).

ALL: hydro1d

vpath %.f90 . ..


odir := _build

# source files
include ../Ghydro.src

# dependencies
include ../Ghydro.dep

ifeq ($(FC),gfortran)
  FFLAGS := -c -O2 -g -C
  FFLAGS += -J $(odir) -I $(odir)
else
  $(error ERROR: compiler $(FC) invalid)
endif


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
	@if [ -d $(odir) ]; then rmdir $(odir); echo "removing $(odir)"; fi
	rm -f hydro1d
