# basic rules for building the source.  Problems are built in
# sub-directories.  Some make tricks taken from BoxLib (CCSE/LBNL).

ALL: hydro1d

# vpath is the list of directories to search for source files.  
vpath %.f90 . ..


# odir is the directory where we put the .o and .mod files as we build
odir := _build


# source files and dependencies
SRC_DIRS := . ../

FSOURCE_DIR := $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.f90))
FSOURCE := $(sort $(notdir $(FSOURCE_DIR)))

$(odir)/deps: $(FSOURCE)                                                                                
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi                                                        
	../util/dep.py --prefix $(odir)/  --search_path "$(SRC_DIRS)" $(FSOURCE) > $(odir)/deps

include $(odir)/deps    


# set the compiler flags for those compilers we know about
ifeq ($(FC),gfortran)
  FFLAGS := -c -O2 -g -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid -finit-real=nan
  FFLAGS += -J $(odir) -I $(odir)

else ifeq ($(FC),mycomp)
  $(info mycomp stuff goes here)

else
  $(error ERROR: compiler $(FC) invalid)
endif


# default rule for building the object files
$(odir)/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	$(FC) $(FFLAGS) -o $@ $< 


# create the list of dependencies for the final build (all the .o files)
OBJECTS = $(addprefix $(odir)/, $(FSOURCE:.f90=.o))


# default target -- this is the executable
hydro1d: $(OBJECTS)
	@echo " "
	@echo "Linking..."
	gfortran -o hydro1d $(OBJECTS)


.PHONY: clean


# targets for cleaning up
clean:
	rm -f $(odir)/*.o $(odir)/*.mod $(odir)/deps

realclean: clean
	@if [ -d $(odir) ]; then rmdir $(odir); echo "removing $(odir)"; fi
	rm -f hydro1d

print-%: ; @echo $* is $($*)                                                                                   
