![Status](https://github.com/NOAA-EMC/NCEPLIBS-sp/workflows/Build%20and%20Test/badge.svg)

# Interpolation Library

The NCEP general interpolation library contains Fortran 90
subprograms to be used for interpolating between nearly all grids used
at NCEP. The library is particularly efficient when interpolating many
fields at one time. The library has been extensively tested with AIX
and Intel Fortran compilers.

This is part of the [NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS)
project.

There are currently six interpolation methods available in the
library:
- bilinear
- bicubic
- neighbor
- budget
- spectral
- neighbor-budget

### Authors

* NCEP/EMC Developers

Code Manager: George Gayno

### Prerequisites

This package requires the [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp) library.

### Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install /path/to/NCEPLIBS-ip2
make -j2
make install
```

### Testing

```
cmake -DENABLE_TESTS=ON -DCMAKE_PREFIX_PATH="/path/to/pfunit;/path/to/NCEPLIBS" /path/to/NCEPLIBS-ip
make
make test
```

### Usage

To use the ip library add `use ip_mod` to your Fortran code. It contains all the necessary public interfaces.

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.
