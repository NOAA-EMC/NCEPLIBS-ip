# Interpolation Library

The NCEP general interpolation library (iplib) contains Fortran 90 subprograms
to be used for interpolating between nearly all (rectilinear) grids used at NCEP.
For more detailed documentation see [README.iplib](README.iplib).

### Prerequisites

Compilers: GNU | Intel | Clang | AppleClang | PGI


### Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install /path/to/NCEPLIBS-ip
make -j2
make install
```


### Version
3.0.0


### Authors

* **[NCEP/EMC](NCEP.List.EMC.nceplibs.Developers@noaa.gov)**
