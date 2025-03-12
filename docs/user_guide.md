## Introduction

The NCEP general interpolation library (NCEPLIBS-ip) contains Fortran 90
subprograms to be used for interpolating between nearly all grids used at NCEP.
The library is particularly efficient when interpolating many fields at one
time. It also contains functionality for interpolating, transforming, and
otherwise manipulating spectral data (these functions were formerly contained in
the NCEPLIBS-sp library).

NCEPLIBS-ip supports compilation with the GNU Compiler Collection (gfortran),
Intel Classic (ifort), and Intel OneAPI (ifx) compilers. In the case of Intel
OneAPI (IntelLLVM), it is recommended to use at least version 2023.2.1 to avoid
any number of compiler issues.

\note Some routines may behave poorly or unpredictably when using 4-byte reals
(libip_4). For instance, there is an ATAN2 function used for polar stereo grids
where for certain regions of certain grids, floating point differences between
4-byte output values (~1e-7) can be amplified into sizable differences in output
field values. Some applications may therefore benefit from the use of 8-byte
reals (libip_d or libip_8).

NCEPLIBS-ip uses several BLAS/LAPACK routines in the splat() subroutine, and
therefore requires an external BLAS/LAPACK provider. In practice, this should
generally be OpenBLAS, which is the [spack-stack](https://github.com/JCSDA/spack-stack)
BLAS/LAPACK provider.

\note When running unit tests, excluding the CTest label SLOW_TEST may be used
to avoid slow-running tests (15 seconds to 1 or 2 minutes depending on
optimization settings). The NO_INPUT_DATA label may be used to run only tests
which do not require files obtained remotely when FTP_TEST_FILES=ON.

## Interpolation

### Interpolation Methods

There are currently six interpolation methods available in the library:
- bilinear
- bicubic
- neighbor
- budget
- spectral
- neighbor-budget

Some of the methods have interpolation sub-options. A few methods have
restrictions on the type of input or output grids.

Several methods can perform interpolation on fields with bitmaps (i.e. some
points on the input grid may be undefined). In this case, the bitmap is
interpolated to the output grid. Only valid input points are used to interpolate
to valid output points. An output bitmap will also be created to locate invalid
data where the output grid extends outside the domain of the input grid.

The driver routines for interpolating scalars and vectors may be found in
ipolates_mod. The interpolation method is chosen via the first argument of these
routines (variable IP). Sub-options are set via the IPOPT array.

#### Bilinear Interpolation Method

Bilinear interpolation is chosen by setting IP=0.

This method has two sub-options:

1. The percent of valid input data required to make output data (the default is
   50%).

2. If valid input data is not found near an a spiral search may be performed.
   The spiral search is only an option for scalar data.

The bilinear method has no restrictions and can interpolate with bitmaps.

#### Bicubic Interpolation Method

Bicubic interpolation is chosen by setting IP=1.

This method has two sub-options:

1. A monotonic constraint option for straight bicubic or for constraining the
   output value to be within the range of the four surrounding input values.
2. The percent of valid input data required to make output data, which defaults
   to 50%.

The bicubic method cannot interpolate data with bitmaps.

#### Neighbor Interpolation Method

Neighbor interpolation is chosen by setting IP=2.

Neighbor interpolation means that the output value is set to the nearest input
value. It would be appropriate for interpolating integer fields such as
vegetation index.

This method has one sub-option: If valid input data is not found near an an
output point, a spiral search is optionally performed.

The neighbor method has no restrictions and can interpolate with bitmaps.

#### Budget Interpolation Method

Budget interpolation is chosen by setting IP=3.

Budget interpolation means a low-order interpolation method that quasi-conserves
area averages. It would be appropriate for interpolating budget fields such as
precipitation.

This method assumes that the field really represents box averages where each box
extends halfway to its neighboring grid point in each direction. The method
actually averages bilinearly interpolated values in a square array of points
distributed within each output grid box.

There are several sub-options:

1. The number of points in the radius of the square array may be set. The
   default is 2, meaning that 25 sample points will be averaged for each output
   value.
  
2. The respective averaging weights for the radius points are adjustable. The
   default is for all weights equal to 1, giving an unweighted average.
  
3. Optionally, one may assume the boxes stretch nearly all the way to each of
   the neighboring grid points and the weights are the adjoint of the bilinear
   interpolation weights.
   
4. The percent of valid input data required to make output data is adjustable.
   The default is 50%.

5. In cases where there is no or insufficient valid input data, a spiral search
   may be invoked to search for the nearest valid data. search square (scalar
   interpolation only).

This method can interpolate data with bitmaps.

#### Spectral Interpolation Method

The spectral interpolation scheme is chosen by setting IP=4.

This method has two sub-options:

1. set the spectral shape (triangular or rhomboidal)

2. set the spectral truncation.

The input grid must be a global cylindrical grid (either Gaussian or
equidistant). This method cannot interpolate data with bitmaps.

Unless the output grid is a global cylindrical grid, a polar stereographic grid
centered at the pole, or a Mercator grid, this method can be quite expensive.

#### Neighbor-Budget Interpolation Method

Neighbor-budget interpolation is chosen by setting IP=6.

This method computes weighted averages of neighbor points arranged in a square
box centered around each output grid point and stretching nearly halfway to each
of the neighboring grid points. The main difference with the budget
interpolation (IP=3) is neighbor vs bilinear interpolation of the square box of
points.

There are the following sub-options:

1. The number of points in the radius of the square array may be set. The
   default is 2, meaning that 25 sample points will be averaged for each output
   value.

2. The respective averaging weights for the radius points are adjustable. The
   default is for all weights equal to 1, giving an unweighted average.

3. The percent of valid input data required to make output data is adjustable.
   The default is 50%.

### Vectors and Scalars

The library can handle two-dimensional vector fields as well as scalar fields.
The input and output vectors are rotated if necessary so that they are either
resolved relative to their defined grid in the direction of increasing x and y
coordinates or resolved relative to eastward and northward directions on the
earth. The rotation is determined by the grid definitions.

Vectors are generally interpolated (by all methods except spectral
interpolation) by moving the relevant input vectors along a great circle to the
output point, keeping their orientations with respect to the great circle
constant, before independently interpolating the respective components. This
ensures that vector interpolation will be consistent over the whole globe
including the poles.

### Grids

The input and output grids are defined by their respective GRIB2 grid definition
template and template number as decoced by the NCEP G2 library. There are six
map projections recognized by the library:

Grid Template Number | Map projection
---------------------|---------------
00 | Equidistant cyclindrical
01 | Rotated equidistant cylindrical
10 | Mercator cyclindrical
20 | Polar stereographic azimuthal
30 | Lambert conformal conical
40 | Gaussian equidistant cyclindrical

If the output grid definition template number is negative, then the output data
may be just a set of station points. In this case, the user must pass the number
of points to be output along with their latitudes and longitudes.

For vector interpolation, the vector rotations parameters must also be passed.
On the other hand, for non-negative output data representation types, the number
of output grid points and their latitudes and longitudes (and the vector
rotation parameters for vector interpolation) are all returned by the
interpolation subprograms.

If an output equidistant cylindrical grid contains multiple pole points, then
the pole points are forced to be self-consistent. That is, scalar fields are
obliged to be constant at the pole and vector components are obliged to exhibit
a wavenumber one variation at the pole.

Generally, only regular grids can be interpolated in this library. However, the
thinned WAFS grids may be expanded to a regular grid (or vice versa) using
subprograms ipxwafs(), ipxwafs2(), or ipxwafs3(). Eta data (with Arakawa "E"
staggering) on the "H" or "V" grid may be expanded to a filled regular grid (or
vice versa) using subprogram ipxetas().

### Return Codes

The return code issued by an interpolation subprogram determines whether it ran
successfully or how it failed. Check nonzero return codes against the docblock
of the respective subprogram.

### Entry point list: interpolation

Scalar and vecotr field interpolation subprograms can be found in the relevant
module documentation:

Name | Function
---- |---------
ipolates_mod | Iredell's polate
bilinear_interp_mod | bilinear interpolation
bicubic_interp_mod | bicubic interpolation
neighbor_interp_mod | neighbor interpolation
budget_interp_mod | budget interpolation
spectral_interp_mod | spectral interpolation
neighbor_budget_interp_mod | neighbor-budget interpolation
polfixs() | make multiple pole scalar values consistent
movect() | move a vector along a great circle
polfixv()| make multiple pole vector values consistent

Grid description section decoders:

Name | Function
---- | --------
gdswzd() | grid description section (GDS) wizard
gdswzd_c() | C wrapper for calling gdswzd
gdswzd_equid_cylind() | GDS wizard for equidistant cyclindrical
gdswzd_mercator() | GDS wizard for mercator cyclindrical
gdswzd_lambert_conf() | GDS wizard for lambert conformal conical
gdswzd_gaussian() | GDS wizard for gaussian cyclindrical
gdswzd_polar_stereo() | GDS wizard for polar stereographic
gdswzd_rot_equid_cylind_egrid() | GDS wizard for rotated equidistant cyclindrical "e" stagger.
gdswzd_rot_equid_cylind() | GDS wizard for rotated equidistant cyclindrical non "e" stagger.
field_pos() | return field position for a given grid point
   
Transform subprograms for special irregular grids:

Name | Function
---- | --------
ipxwafs() |  expand or contract wafs grids
ipxwafs2() |  expand or contract wafs grids
ipxwafs3() |  expand or contract wafs grids

## Spectral Transformation & Processing

The library's spectral processing subroutines can handle both scalar and
two-dimensional vector fields. Each vector field will be represented in spectral
space appropriately by its respective spherical divergence and curl (vorticity),
thus avoiding the pole problems associated with representing components
separately.

Some of the functions performed by the library are spectral interpolations
between two grids, spectral truncations in place on a grid, and basic spectral
transforms between grid and wave space. Only global Gaussian or global
equidistant cylindrical grids are allowed for transforming into wave space.
There are no such restricitions on grids for transforming from wave space.
However, there are special fast entry points for transforming wave space to
polar stereographic and Mercator grids as well as the aforementioned cylindrical
grids.

The indexing of the cylindrical transform grids is totally general. The grids
may run north to south or south to north; they may run east to west or west to
east; they may start at any longitude as long as the prime meridian is on the
grid; they may be dimensioned in any order (e.g. (i,j,k), (k,j,i),
(i,k,nfield,j), etc.). Furthermore, the transform may be performed on only some
of the latitudes at one time as long as both hemisphere counterparts are
transformed at the same time (as in the global spectral model). The grid
indexing will default to the customary global indexing, i.e. north to south,
east to west, prime meridian as first longitude, and (i,j,k) order.

The wave space may be either triangular or rhomboidal in shape. Its internal
indexing is strictly "IBM order", i.e. zonal wavenumber is the slower index with
the real and imaginary components always paired together. The imaginary
components of all the zonally symmetric modes should always be zero, as should
the global mean of any divergence and vorticity fields. The stride between the
start of successive wave fields is general, defaulting to the computed length of
each field.



### Entry Point List: Spectral Interpolation & Transformation

Spectral interpolations or truncations between grid and grid

   Name        | Function
   ----        | --------
   sptrun()    | Spectrally truncate gridded scalar fields
   sptrunv()   | Spectrally truncate gridded vector fields
   sptrung()   | Spectrally interpolate scalars to stations
   sptrungv()  | Spectrally interpolate vectors to stations
   sptruns()   | Spectrally interpolate scalars to polar stereo
   sptrunsv()  | Spectrally interpolate vectors to polar stereo
   sptrunm()   | Spectrally interpolate scalars to Mercator
   sptrunmv()  | Spectrally interpolate vectors to Mercator

Spectral transforms between wave and grid

   Name        | Function
   ----        | ------------------------------------------------------------------
   sptran()    | Perform a scalar spherical transform
   sptranv()   | Perform a vector spherical transform
   sptrand()   | Perform a gradient spherical transform
   sptgpt()    | Transform spectral scalar to station points
   sptgptv()   | Transform spectral vector to station points
   sptgptd()   | Transform spectral to station point gradients
   sptgps()    | Transform spectral scalar to polar stereo
   sptgpsv()   | Transform spectral vector to polar stereo
   sptgpsd()   | Transform spectral to polar stereo gradients
   sptgpm()    | Transform spectral scalar to Mercator
   sptgpmv()   | Transform spectral vector to Mercator
   sptgpmd()   | Transform spectral to Mercator gradients

Spectral transform utilities

   Name        | Function
   ----        | ------------------------------------------------------------------
   spwget()    | Get wave-space constants
   splat()     | Compute latitude functions
   speps()     | Compute utility spectral fields
   splegend()  | Compute Legendre polynomials
   spanaly()   | Analyze spectral from Fourier
   spsynth()   | Synthesize Fourier from spectral
   spdz2uv()   | Compute winds from divergence and vorticity
   spuv2dz()   | Compute divergence and vorticity from winds
   spgradq()   | Compute gradient in spectral space
   splaplac()  | Compute Laplacian in spectral space


## Examples: Interpolation Routines

Example 1. Read a grib 2 file of scalar data on a global regular 1-deg lat/lon
grid and call ipolates to interpolate it to NCEP standard grid 218, a lambert
conformal grid. Uses the NCEP G2 library to degrib the data.

\code{fortran}
 program example_1

use ip_mod
 use grib_mod  ! ncep grib 2 library

 implicit none

 character(len=100)      :: input_file

 integer                 :: iunit, iret, lugi
 integer                 :: mi, mo, no
 integer, allocatable    :: ibi(:), ibo(:)
 integer                 :: ip, ipopt(20)
 integer                 :: j, jdisc, jpdtn, jgdtn, k, km
 integer                 :: jids(200), jgdt(200), jpdt(200)
 integer                 :: idim_input, jdim_input
 integer                 :: idim_output, jdim_output

 logical                 :: unpack
 logical*1, allocatable  :: input_bitmap(:,:), output_bitmap(:,:)

 real, allocatable       :: input_data(:,:)
 real, allocatable       :: output_rlat(:), output_rlon(:)
 real, allocatable       :: output_data(:,:)

 type(gribfield)         :: gfld_input

!---------------------------------------------------------------------------
! the output grid specs. this is ncep grid 218, a lambert conformal
! grid. the grid definition information is stored in section 3
! of a grib 2 message.
!---------------------------------------------------------------------------

 integer, parameter :: igdtnum218 = 30 ! grid definition template number.
                                       ! "30" is lambert conformal.
 integer, parameter :: igdtlen218 = 22 ! number of array elements needed
                                       ! for a lambert conf. grid definition
                                       ! template.
 integer     :: igdtmpl218(igdtlen218) ! the grid definition template.
                                       ! the entries are:
                                       ! 1 -shape of earth, oct 15
                                       ! 2 -scale factor, spherical earth, oct 16
                                       ! 3 -scaled value, spherical earth, octs 17-20
                                       ! 4 -scale factor, major axis of
                                       !    elliptical earth, oct 21
                                       ! 5 -scaled value of major axis of
                                       !    elliptical earth, octs 22-25
                                       ! 6 -scale factor, minor axis of
                                       !    elliptical earth, oct 26
                                       ! 7 -scaled value of minor axis of
                                       !    elliptical earth, octs 27-30
                                       ! 8 -number points along x-axis, octs 31-34
                                       ! 9 -number points along y-axis, octs 35-38
                                       ! 10-latitude of first point, octs 39-42
                                       ! 11-longitude of first point, octs 43-46
                                       ! 12-resolution and component flags, oct 47
                                       ! 13-latitude where grid lengths specified, 
                                       !    octs 48-51
                                       ! 14-longitude parallel to y-axis, octs 52-55
                                       ! 15-x-direction grid length, octs 56-59
                                       ! 16-y-direction grid length, octs 60-63
                                       ! 17-projection center flag, oct 64
                                       ! 18-scanning mode, oct 65
                                       ! 19-first tangent latitude from pole, octs 66-69
                                       ! 20-second tangent latitude from pole, octs 70-73
                                       ! 21-latitude of south pole, octs 74-77
                                       ! 22-longitude of south pole, octs 78-81

 integer, parameter :: missing=b'11111111111111111111111111111111'
 data igdtmpl218 / 6, 255, missing, 255, missing, 255, missing, 614, 428, &
                  12190000, 226541000, 56, 25000000, 265000000, &
                  12191000, 12191000, 0, 64, 25000000, 25000000, -90000000, 0/

!---------------------------------------------------------------------------
! open the grib 2 file containing data to be interpolated. for this
! example, there are two data records.
!---------------------------------------------------------------------------

 iunit=9
 input_file="${path}/input.data.grib2"
 call baopenr (iunit, input_file, iret)

!---------------------------------------------------------------------------
! prep for call to g2 library to degrib data. the data are on a regular
! lat/lon grid with i/j dimension of 360/181. 
!---------------------------------------------------------------------------

 idim_input = 360  ! the i/j dimensions of input grid
 jdim_input = 181
 mi         = idim_input * jdim_input   ! total number of pts, input grid

 jdisc   = -1         ! search for any discipline
 jpdtn   = -1         ! search for any product definition template number
 jgdtn   =  0         ! search for grid definition template number 0 - regular lat/lon grid
 jids    = -9999      ! array of values in identification section, set to wildcard
 jgdt    = -9999      ! array of values in grid definition template 3.m
 jgdt(8) = idim_input ! search for grid with i/j of 360/181
 jgdt(9) = jdim_input
 jpdt    = -9999      ! array of values in product definition template 4.n
 unpack  = .true.     ! unpack data
 lugi    = 0          ! no index file

 nullify(gfld_input%idsect)
 nullify(gfld_input%local)
 nullify(gfld_input%list_opt)
 nullify(gfld_input%igdtmpl)  ! holds the grid definition template information
 nullify(gfld_input%ipdtmpl)
 nullify(gfld_input%coord_list)
 nullify(gfld_input%idrtmpl)
 nullify(gfld_input%bmap)     ! holds the bitmap
 nullify(gfld_input%fld)      ! holds the data

!---------------------------------------------------------------------------
! degrib the data.  non-zero "iret" indicates a problem during degrib.
!---------------------------------------------------------------------------

 km = 2                  ! number of records to interpolate

 allocate(ibi(km))
 allocate(input_bitmap(mi,km))
 allocate(input_data(mi,km))

 do j = 0, (km-1)    ! number of records to skip

   call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
               unpack, k, gfld_input, iret)

   if (iret /= 0) stop

!---------------------------------------------------------------------------
! does input data have a bitmap?
!---------------------------------------------------------------------------

   if (gfld_input%ibmap==0) then  ! input data has bitmap
     ibi(k)            = 1        ! tell ipolates to use bitmap
     input_bitmap(:,k) = gfld_input%bmap
   else                           ! no bitmap, data everywhere
     ibi(k)            = 0        ! tell ipolates there is no bitmap
     input_bitmap(:,k) = .true.
   endif

   input_data(:,k) = gfld_input%fld  ! the input data field
 
 enddo

 call baclose (iunit, iret)

!---------------------------------------------------------------------------
! setup arguments for ipolates (scalar interpolation) call.
!---------------------------------------------------------------------------

 ip       = 0                         ! bilinear interpolation
 ipopt    = 0                         ! options for bilinear:
 ipopt(1) = 75                        ! set minimum mask to 75%

!---------------------------------------------------------------------------
! the i/j dimensions of the output grid.
!---------------------------------------------------------------------------

 idim_output = igdtmpl218(8)
 jdim_output = igdtmpl218(9)
 mo          = idim_output * jdim_output ! total number of output pts

!---------------------------------------------------------------------------
! will hold the latitude, longitude, data and bitmap on the output grid,
! which are computed in ipolates.
!---------------------------------------------------------------------------

 allocate (ibo(km))              ! bitmap flags on output grid
 allocate (output_rlat(mo))
 allocate (output_rlon(mo))
 allocate (output_data(mo,km))
 allocate (output_bitmap(mo,km))

!---------------------------------------------------------------------------
! call ipolates to interpolate scalar data.  non-zero "iret" indicates
! a problem.
!---------------------------------------------------------------------------

 call ipolates(ip, ipopt, gfld_input%igdtnum, gfld_input%igdtmpl, &
               gfld_input%igdtlen, igdtnum218, igdtmpl218, igdtlen218, &
               mi, mo, km, ibi, input_bitmap, input_data, no, output_rlat, &
               output_rlon, ibo, output_bitmap, output_data, iret)

 if (iret /= 0) stop

!---------------------------------------------------------------------------
! write interpolated data to file.  if ipolates computed a bitmap (ibo==1) 
! for the output grid, one may mask out this data with a flag value.
!---------------------------------------------------------------------------

 open (10, file="./output.bin", access='direct', recl=idim_output*jdim_output*4)

 do k = 1, km
   if(ibo(k)==1) where (.not. output_bitmap(:,k)) output_data(:,k) = -999.
   write(10, rec=k) output_data(:,k)
 enddo
 write(10, rec=km+1) output_rlat
 write(10, rec=km+2) output_rlon

 close(10)

 end program example_1
\endcode

Example 2.  Read a grib 2 file of u/v wind data on a global regular
            1-deg lat/lon grid and call ipolatev to interpolate
            it to four random station points.  Uses the NCEP
            G2 library to degrib the data.

\code{fortran}
 program example_2

 use grib_mod  ! ncep grib 2 library

 implicit none

 character(len=100)      :: input_file

 integer                 :: iunit, iret, lugi
 integer                 :: mi, mo, no
 integer                 :: ibi, ibo
 integer                 :: ip, ipopt(20)
 integer                 :: j, jdisc, jpdtn, jgdtn, k, km
 integer                 :: jids(200), jgdt(200), jpdt(200)
 integer                 :: idim_input, jdim_input

 logical                 :: unpack
 logical*1, allocatable  :: input_bitmap(:), output_bitmap(:)

 real, allocatable       :: input_u_data(:), input_v_data(:)
 real, allocatable       :: output_rlat(:), output_rlon(:)
 real, allocatable       :: output_crot(:), output_srot(:)
 real, allocatable       :: output_u_data(:), output_v_data(:)

 type(gribfield)         :: gfld_input

!---------------------------------------------------------------------------
! the output "grid" is a series of random station points.  in this case,
! set the grid definition template number of a negative number.
! the grid definition template array information is not used, so set
! to a flag value.
!---------------------------------------------------------------------------

 integer, parameter      :: igdtnumo = -1 
 integer, parameter      :: igdtleno =  1
 integer                 :: igdtmplo(igdtleno)

 data igdtmplo / -9999 /

!---------------------------------------------------------------------------
! open the grib 2 file containing data to be interpolated.  for this
! example, there is one record of u-wind and v-wind.
!---------------------------------------------------------------------------

 iunit=9
 input_file="./reg_tests/copygb2/data/uv_wind.grb2"
 call baopenr (iunit, input_file, iret)

!---------------------------------------------------------------------------
! prep for call to g2 library to degrib data. the data are on a regular
! lat/lon grid with i/j dimension of 360/181.  
!---------------------------------------------------------------------------

 idim_input = 360  ! the i/j dimensions of input grid
 jdim_input = 181
 mi         = idim_input * jdim_input   ! total number of pts, input grid

 jdisc   = -1         ! search for any discipline
 jpdtn   = -1         ! search for any product definition template number
 jgdtn   =  0         ! search for grid definition template number 0 - regular lat/lon grid
 jids    = -9999      ! array of values in identification section, set to wildcard
 jgdt    = -9999      ! array of values in grid definition template 3.m
 jgdt(8) = idim_input ! search for grid with i/j of 360/181
 jgdt(9) = jdim_input
 jpdt    = -9999      ! array of values in product definition template 4.n
 unpack  = .true.     ! unpack data
 lugi    = 0          ! no index file

 nullify(gfld_input%idsect)
 nullify(gfld_input%local)
 nullify(gfld_input%list_opt)
 nullify(gfld_input%igdtmpl)  ! holds the grid definition template information
 nullify(gfld_input%ipdtmpl)
 nullify(gfld_input%coord_list)
 nullify(gfld_input%idrtmpl)
 nullify(gfld_input%bmap)     ! holds the bitmap
 nullify(gfld_input%fld)      ! holds the data

!---------------------------------------------------------------------------
! degrib the data.  non-zero "iret" indicates a problem during degrib.
!---------------------------------------------------------------------------

 allocate(input_bitmap(mi))
 allocate(input_u_data(mi))
 allocate(input_v_data(mi))

!---------------------------------------------------------------------------
! read u-wind record.
!---------------------------------------------------------------------------

 j = 0
 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld_input, iret)

 if (iret /= 0) stop

!---------------------------------------------------------------------------
! does input data have a bitmap?
!---------------------------------------------------------------------------

 if (gfld_input%ibmap==0) then  ! input data has bitmap
   ibi          = 1             ! tell ipolates to use bitmap
   input_bitmap = gfld_input%bmap
 else                           ! no bitmap, data everywhere
   ibi          = 0             ! tell ipolates there is no bitmap
   input_bitmap = .true.
 endif

 input_u_data = gfld_input%fld  ! the input u-wind data
 
!---------------------------------------------------------------------------
! read v-wind record.
!---------------------------------------------------------------------------

 j = 1
 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld_input, iret)

 if (iret /= 0) stop

 input_v_data = gfld_input%fld  ! the input v-wind data

 call baclose (iunit, iret)

!---------------------------------------------------------------------------
! setup arguments for ipolatev (vector interpolation) call.
!---------------------------------------------------------------------------

 km       = 1                         ! number of records to interpolate
 ip       = 0                         ! bilinear interpolation
 ipopt    = 0                         ! options for bilinear:
 ipopt(1) = 75                        ! set minimum mask to 75%

!---------------------------------------------------------------------------
! interpolate to four random station points.
!---------------------------------------------------------------------------

 mo = 4
 no = mo

!---------------------------------------------------------------------------
! when interpolating to random station points, need to pass to ipolatev
! their latitude, longitude and the sines and cosines of the vector
! rotation angles.  the vector rotation is defined:
!
! ugrid=crot*uearth-sort*vearth
! vgrid=srot*uearth+cort*vearth
!---------------------------------------------------------------------------

 allocate (output_rlat(mo))
 allocate (output_rlon(mo))
 allocate (output_srot(mo))
 allocate (output_crot(mo))
 allocate (output_u_data(mo))
 allocate (output_v_data(mo))
 allocate (output_bitmap(mo))

 output_rlat(1) = 45.0
 output_rlon(1) = -100.0
 output_rlat(2) = 35.0
 output_rlon(2) = -100.0
 output_rlat(3) = 40.0
 output_rlon(3) = -90.0
 output_rlat(4) = 35.0
 output_rlon(4) = -120.0

 output_srot = 0.0   ! no turning of wind
 output_crot = 1.0

!---------------------------------------------------------------------------
! call ipolatev to interpolate vector data.  non-zero "iret" indicates
! a problem.
!---------------------------------------------------------------------------
 
 call ipolatev(ip, ipopt, gfld_input%igdtnum, gfld_input%igdtmpl, &
               gfld_input%igdtlen, igdtnumo, igdtmplo, igdtleno, &
               mi, mo, km, ibi, input_bitmap, input_u_data, input_v_data, &
               no, output_rlat, output_rlon, output_crot, output_srot, &
               ibo, output_bitmap, output_u_data, output_v_data, iret)

 if (iret /= 0) stop

 do k = 1, mo
   print*,'station point ',k,' latitude ',output_rlat(k),' longitude ', &
   output_rlon(k), ' u-wind ', output_u_data(k), ' v-wind ', output_v_data(k)
 enddo

 end program example_2
\endcode

## Examples: Spectral Processing & Transformation

Example 1. Interpolate heights and winds from a latlon grid
           to two antipodal polar stereographic grids.
           Subprograms GETGB and PUTGB from w3lib are referenced.

\code{fortran}
c  unit number 11 is the input latlon grib file
c  unit number 31 is the input latlon grib index file
c  unit number 51 is the output northern polar stereographic grib file
c  unit number 52 is the output southern polar stereographic grib file
c  nominal spectral truncation is r40
c  maximum input gridsize is 360x181
c  maximum number of levels wanted is 12
      parameter(lug=11,lui=31,lun=51,lus=52)
      parameter(iromb=1,maxwv=40,jf=360*181,kx=12)
      integer kp5(kx),kp6(kx),kp7(kx)
      integer kpo(kx)
      data kpo/1000,850,700,500,400,300,250,200,150,100,70,50/
c height
      km=12
      kp5=7
      kp6=100
      kp7=kpo
      call gs65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c winds
      km=12
      kp5=33
      kp6=100
      kp7=kpo
      call gv65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c
      stop
      end
c
      subroutine gs65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c  interpolates a scalar field using spectral transforms.
      integer kp5(km),kp6(km),kp7(km)
c  output grids are 65x65 (381 km true at latitide 60).
c  nh grid oriented at 280E; sh grid oriented at 100E.
      parameter(nph=32,nps=2*nph+1,npq=nps*nps)
      parameter(true=60.,xmesh=381.e3,orient=280.)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi) 
      real gn(npq,km),gs(npq,km)
      integer jpds(25),jgds(22),kpds(25,km),kgds(22,km)
      logical lb(jf)
      real f(jf,km)
c
      g2=((1.+sin(abs(true)/dpr))*rerth/xmesh)**2
      r2=2*nph**2
      rlatn1=dpr*asin((g2-r2)/(g2+r2))
      rlonn1=mod(orient+315,360.)
      rlats1=-rlatn1
      rlons1=mod(rlonn1+270,360.)
      jpds=-1
      do k=1,km
        jpds(5)=kp5(k)
        jpds(6)=kp6(k)
        jpds(7)=kp7(k)
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,f(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
      enddo
      idrt=kgds(1,1)
      imax=kgds(2,1)
      jmax=kgds(3,1)
c
      call sptruns(iromb,maxwv,idrt,imax,jmax,km,nps,
     &             0,0,0,jf,0,0,0,0,true,xmesh,orient,f,gn,gs)
c
      do k=1,km
        kpds(3,k)=27
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlatn1*1.e3)
        kgds(5,k)=nint(rlonn1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(orient*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=0
        kgds(11,k)=64
        call putgb(lun,npq,kpds(1,k),kgds(1,k),lb,gn(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=28
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlats1*1.e3)
        kgds(5,k)=nint(rlons1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(mod(orient+180,360.)*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=128
        kgds(11,k)=64
        call putgb(lus,npq,kpds(1,k),kgds(1,k),lb,gs(1,k),iret)
      enddo
c
      end
c
      subroutine gv65(lug,lui,lun,lus,jf,km,kp5,kp6,kp7,iromb,maxwv)
c  interpolates a vector field using spectral transforms.
      integer kp5(km),kp6(km),kp7(km)
c  output grids are 65x65 (381 km true at latitide 60).
c  nh grid oriented at 280E; sh grid oriented at 100E.
c  winds are rotated to be relative to grid coordinates.
      parameter(nph=32,nps=2*nph+1,npq=nps*nps)
      parameter(true=60.,xmesh=381.e3,orient=280.)
      parameter(rerth=6.3712e6)
      parameter(pi=3.14159265358979,dpr=180./pi) 
      real un(npq,km),vn(npq,km),us(npq,km),vs(npq,km)
      integer jpds(25),jgds(22),kpds(25,km),kgds(22,km)
      logical lb(jf)
      real u(jf,km),v(jf,km)
c
      g2=((1.+sin(abs(true)/dpr))*rerth/xmesh)**2
      r2=2*nph**2
      rlatn1=dpr*asin((g2-r2)/(g2+r2))
      rlonn1=mod(orient+315,360.)
      rlats1=-rlatn1
      rlons1=mod(rlonn1+270,360.)
      jpds=-1
      do k=1,km
        jpds(5)=kp5(k)
        jpds(6)=kp6(k)
        jpds(7)=kp7(k)
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,u(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
        jpds=kpds(:,k)
        jgds=kgds(:,k)
        jpds(5)=jpds(5)+1
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,v(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
      enddo
      idrt=kgds(1,1)
      imax=kgds(2,1)
      jmax=kgds(3,1)
c
      call sptrunsv(iromb,maxwv,idrt,imax,jmax,km,nps,
     &              0,0,0,jf,0,0,0,0,true,xmesh,orient,u,v,
     &              .true.,un,vn,us,vs,.false.,dum,dum,dum,dum,
     &              .false.,dum,dum,dum,dum)
c
      do k=1,km
        kpds(3,k)=27
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlatn1*1.e3)
        kgds(5,k)=nint(rlonn1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(orient*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=0
        kgds(11,k)=64
        kpds(5,k)=kp5(k)
        call putgb(lun,npq,kpds(1,k),kgds(1,k),lb,un(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=27
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlatn1*1.e3)
        kgds(5,k)=nint(rlonn1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(orient*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=0
        kgds(11,k)=64
        kpds(5,k)=kp5(k)+1
        call putgb(lun,npq,kpds(1,k),kgds(1,k),lb,vn(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=28
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlats1*1.e3)
        kgds(5,k)=nint(rlons1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(mod(orient+180,360.)*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=128
        kgds(11,k)=64
        kpds(5,k)=kp5(k)
        call putgb(lus,npq,kpds(1,k),kgds(1,k),lb,us(1,k),iret)
      enddo
      do k=1,km
        kpds(3,k)=28
        kgds(1,k)=5
        kgds(2,k)=nps
        kgds(3,k)=nps
        kgds(4,k)=nint(rlats1*1.e3)
        kgds(5,k)=nint(rlons1*1.e3)
        kgds(6,k)=8
        kgds(7,k)=nint(mod(orient+180,360.)*1.e3)
        kgds(8,k)=nint(xmesh)
        kgds(9,k)=nint(xmesh)
        kgds(10,k)=128
        kgds(11,k)=64
        kpds(5,k)=kp5(k)+1
        call putgb(lus,npq,kpds(1,k),kgds(1,k),lb,vs(1,k),iret)
      enddo
c
      end
\endcode

Example 2. Spectrally truncate winds in place on a latlon grid.

\code{fortran}
c  unit number 11 is the input latlon grib file
c  unit number 31 is the input latlon grib index file
c  unit number 51 is the output latlon grib file
c  nominal spectral truncation is r40
c  maximum input gridsize is 360x181
c  maximum number of levels wanted is 12
      parameter(lug=11,lui=31,luo=51)
      parameter(iromb=1,maxwv=40,jf=360*181,kx=12)
      integer kp5(kx),kp6(kx),kp7(kx)
      integer kpo(kx)
      data kpo/1000,850,700,500,400,300,250,200,150,100,70,50/
c winds
      km=12
      kp5=33
      kp6=100
      kp7=kpo
      call gvr40(lug,lui,luo,jf,km,kp5,kp6,kp7,iromb,maxwv)
c
      stop
      end
c
      subroutine gvr40(lug,lui,luo,jf,km,kp5,kp6,kp7,iromb,maxwv)
c  interpolates a vector field using spectral transforms.
      integer kp5(km),kp6(km),kp7(km)
      integer jpds(25),jgds(22),kpds(25,km),kgds(22,km)
      logical lb(jf)
      real u(jf,km),v(jf,km)
c
      jpds=-1
      do k=1,km
        jpds(5)=kp5(k)
        jpds(6)=kp6(k)
        jpds(7)=kp7(k)
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,u(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
        jpds=kpds(:,k)
        jgds=kgds(:,k)
        jpds(5)=jpds(5)+1
        j=0
        call getgb(lug,lui,jf,j,jpds,jgds,kf,j,kpds(1,k),kgds(1,k),
     &             lb,v(1,k),iret)
        if(iret.ne.0) call exit(1)
        if(mod(kpds(4,k)/64,2).eq.1) call exit(2)
      enddo
      idrt=kgds(1,1)
      imax=kgds(2,1)
      jmax=kgds(3,1)
c
      call sptrunv(iromb,maxwv,idrt,imax,jmax,idrt,imax,jmax,km,
     &             0,0,0,jf,0,0,jf,0,u,v,.true.,u,v,
     &             .false.,dum,dum,.false.,dum,dum)
c
      do k=1,km
        kpds(5,k)=kp5(k)
        call putgb(luo,kf,kpds(1,k),kgds(1,k),lb,u(1,k),iret)
      enddo
      do k=1,km
        kpds(5,k)=kp5(k)+1
        call putgb(luo,kf,kpds(1,k),kgds(1,k),lb,v(1,k),iret)
      enddo
c
      end
\endcode

Example 3. Compute latlon temperatures from spectral temperatures and
           compute latlon winds from spectral divergence and vorticity.

\code{fortran}
c  unit number 11 is the input sigma file
c  unit number 51 is the output latlon file
c  nominal spectral truncation is t62
c  output gridsize is 144x73
c  number of levels is 28
      parameter(iromb=0,maxwv=62)
      parameter(idrt=0,im=144,jm=73)
      parameter(levs=28)
      parameter(mx=(maxwv+1)*((iromb+1)*maxwv+2)/2)
      real t(mx,levs),d(mx,levs),z(mx,levs)
      real tg(im,jm,km),ug(im,jm,km),vg(im,jm,km)
c  temperature
      do k=1,4
        read(11)
      enddo
      do k=1,levs
        read(11) (t(m,k),m=1,mx)
      enddo
      call sptran(iromb,maxwv,idrt,im,jm,levs,0,0,0,0,0,0,0,0,1,
     &            t,tg(1,1,1),tg(1,jm,1),1)
      call sptran(
      do k=1,levs
        write(51) ((tg(i,j,k),i=1,im),j=1,jm)
      enddo
c  winds
      do k=1,levs
        read(11) (d(m,k),m=1,mx)
        read(11) (z(m,k),m=1,mx)
      enddo
      call sptranv(iromb,maxwv,idrt,im,jm,levs,0,0,0,0,0,0,0,0,1,
     &             d,z,ug(1,1,1),ug(1,jm,1),vg(1,1,1),vg(1,jm,1),1)
      do k=1,levs
        write(51) ((ug(i,j,k),i=1,im),j=1,jm)
        write(51) ((vg(i,j,k),i=1,im),j=1,jm)
      enddo
      end
\endcode
