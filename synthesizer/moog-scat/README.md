MOOG_SCAT
==========

MOOG_SCAT .

Authors
-------

- **Jennifer Sobeck** (CFHT/UW)
- **Chris Sneden** (UT Austin)

Installation
------------

MOOG_SCAT is a Fortran 77 (F77) code.  To install the latest version of the code, employ one of the OS-specific Makefiles (Makefile_SCAT.macos; Makefile_SCAT.rh) and 
use the following command:

    make -f Makefile*

Accordingly, a Fortran compiler is necessary (succesful MOOG_SCAT compilation has been achieved with a ``gfortran`` compiler).  In addition, MOOG_SCAT relies 
completely upon [SuperMongo](https://www.astro.princeton.edu/~rhl/sm/) for plotting/visualziation (latest release 2.4.1).	
    
Documentation
-------------

Refer to the [Documentation](https://moog-scat.readthedocs.io/en/latest/) for information on how
to install and use MOOG_SCAT.
