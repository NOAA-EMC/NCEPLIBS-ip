## Description

Unit test for the NCEP IPOLATES library (IPLIB).

This test invokes IPOLATES to interpolate scalar and
vector data to several grids of various map projections
using all IPOLATES interpolation options.  The output 
is compared to baseline data and a summary of differences is
sent to standard output.


## How the Test Works

There are two separate programs to test the scalar and vector interpolation.
The scalar program calls routine 'ipolates' to interpolate a global field
of snow albedo to several grids of various map projections and
using all IPOLATES interpolation options.  Likewise, the vector program
calls routine 'ipolatev' to interpolate a global field of 500mb u/w wind.
The specific tests are:

  - grid 3 (global one-deg lat/lon) using bilinear (IP option "0")
  - grid 8 (mercator) using bicubic (IP option "1")
  - grid 127 (gaussian lat/lon) using neighbor (IP option "2")
  - grid 203 (rotated lat/lon "E") using budget (IP option "3")
  - grid 205 (rotated lat/lon "B") using spectral (IP option "4")
  - grid 212 (polar stereographic) using neighbor-budget (IP option "6")
  - grid 218 (lambert conformal) using bilinear (IP option "0")

Some grid numbers refer to an NCEP grib 1 standard grid.  Others
refer to the grib 1 data representation type.

The output from each call to ipolates and ipolatev is compared to its
baseline counterpart in the ./baseline data directory.  This baseline 
data was created using the Intel Fortran compiler.  Differences from
the baseline data are computed and sent to standard output. For this
unit test to pass the differences should be "small".


## Data

### Input Data
Contains the input scalar and vector data.  The 
data are in binary little endian format:

  - data/scalar/global_snoalb.bin    (global snow albedo)
  - data/vector/global_uv_wind.bin   (global 500mb u/v wind)

### Baseline Data
Contains the baseline set of scalar and vector data 
interpolated to each target grid. Binary, little endian
format.  Output from each call to routine ipolates and
ipolatev is compared to its baseline counterpart.  There
are sub-directories for the scalar and vector data.
And under these sub-directories, there are sub-directories
for the single (4_byte_bin) and mixed/double (8_byte_bin)
versions of the IPOLATES library.  The file names
contain the grid number and IPOLATES interpolation
option (defined below); i.e.:

- grid${gridnum}.opt${ip_option}.bin_4/8

Grads control files to view each file are located in the
grads subdirectory.  These files do not account
for map projection.  Rather, they display the data
as a 2-D field with no interpolation.  So when
viewing, do a "set mproj off" during your grads session.
                  

## How to Run

From the build directory simply run `make test`

Alternatively, to run all tests and show stdout run `ctest --verbose`

The individual test executables can be run from `<build_directory>/tests`. `test_scalar` and `test_vector` take two command-line arguments for the grid type and interpolation scheme.


