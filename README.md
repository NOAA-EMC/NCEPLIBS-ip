![Status](https://github.com/NOAA-EMC/NCEPLIBS-ip/workflows/developer/badge.svg)

# Interpolation Library

The NCEP general interpolation library contains Fortran 90 subprograms to be
used for interpolating between nearly all grids used at NCEP. The library is
particularly efficient when interpolating many fields at one time. It also
contains routines for spectral transforms and other processing, including those
previously contained in the NCEPLIBS-sp library.

This is part of the [NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS)
project.

There are currently six interpolation methods available in the library:
- bilinear
- bicubic
- neighbor
- budget
- spectral
- neighbor-budget

For full documentation see https://noaa-emc.github.io/NCEPLIBS-ip/.

To submit bug reports, feature requests, or other code-related issues including installation and usage questions, please create a [GitHub issue](https://github.com/NOAA-EMC/NCEPLIBS-ip/issues). For general NCEPLIBS inquiries, contact [Ed Hartnett](mailto:edward.hartnett@noaa.gov) (secondary point of contact [Alex Richert](mailto:alexander.richert@noaa.gov)).

### Authors

* NCEP/EMC Developers

Code Manager: [Alex Richert](mailto:alexander.richert@noaa.gov)

### Prerequisites

This package requires a BLAS/LAPACK library to provide several LU decomposition-related
routines, and requires CMake (version 3.15+) to build. In spack-stack, OpenBLAS will
be used as the provider of BLAS/LAPACK routines for all compilers.

### Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install /path/to/NCEPLIBS-ip
make -j2
make test # or 'ctest --verbose'; use 'ctest -LE SLOW_TEST' to skip long-running tests
make install
```

### Usage

Most routines and any public interfaces required can be accessed by adding `use
ip_mod` to your Fortran code. Most spectral transform and processing subroutines
can be accessed by calling them in your code (no `use` statement) and linking
to the ip library at build time.

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an "as is" basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their endorsement,
recommendation or favoring by the Department of Commerce. The Department of
Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used
in any manner to imply endorsement of any commercial product or activity by DOC
or the United States Government.
