# Interpolation Library

The NCEP general interpolation library (iplib) contains Fortran 90 subprograms
to be used for interpolating between nearly all (rectilinear) grids used at NCEP.
For more detailed documentation see https://noaa-emc.github.io/NCEPLIBS-ip/.

### Authors

* NCEP/EMC Developers

Code Manager: George Gayno

### Prerequisites

This package requires the [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp) library.

### Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install /path/to/NCEPLIBS-ip
make -j2
make install
```

### Testing

Testing requires [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit).

```
cmake -DENABLE_TESTS=ON -DCMAKE_PREFIX_PATH="/path/to/pfunit;/path/to/NCEPLIBS" /path/to/NCEPLIBS-ip
make -j2
make test
```

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
