# hydro1d 

*A Fortran 95 1-d hydrodynamics code*

http://zingale.github.io/hydro1d/


Important: this project is no longer developed.  It has largely been
replaced by a python implementation of PPM: https://zingale.github.io/ppmpy/


## About

`hydro1d` is a simple one-dimensional finite-volume Eulerian
hydrodynamics code that implements piecewise constant, piecewise
linear, and piecewise parabolic (PPM) reconstruction, and supports
both Cartesian and spherical geometries.

It is written in modern Fortran (a F2003+ compiler is needed).

At the moment, a constant-gamma equation of state is assumed.


## Problems

Individual problems are built in subdirectories.  For example,
to build and run the Sod shock tube problem, do:

* `cd hydro1d/sod`

* `make`

* `./hydro1d inputs-sod-xp`

As the code is built, object files and modules will be output into the
`_build` subdirectory.  `make realclean` will clean up the objects.
Things are setup for gfortran by default -- you will need to edit the
`Ghydro.mak` with different options for different compilers.  Some bits
of Fortran 2003 and 2008 are used, so an up-to-date compiler is
needed.


## Runtime Parameters

A number of runtime options can be specified -- look in `params.f90`
for the available runtime parameters.  These are set in a namelist in
the inputs file.  Problems can specify their own runtime parameters
(with a `probparams.f90` in the problem directory).  These have a
separate namelist in the same inputs file.  See the sod problem for
examples.


## Development

 * The PPM implementation should be synced up to what was done in 
   Castro's PPM -- in particular, support for using the reference
   state in the eigenvectors

 * The boundary filling stuff should be generalized so we can use the
   same logic for gravity as for the main conserved state.

 * The monopole gravity is not quite right for the conservative update.
   We need to time-center the gravitational acceleration.
   
 * We do not trace under the geometry or gravity source terms in the
   interface state construction
   
