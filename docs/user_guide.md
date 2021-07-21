# Documentation of the general interpolation library - iplib

## Introduction

The NCEP general interpolation library (ip2lib) contains Fortran 90 subprograms
to be used for interpolating between nearly all grids used at NCEP.
The library is particularly efficient when interpolating many fields at one time.

There are currently six interpolation methods available in the library:
bilinear, bicubic, neighbor, budget, spectral and neighbor-budget.

Some of the methods have interpolation sub-options.  A few methods
have restrictions on the type of input or output grids.
Also, several methods can perform interpolation on fields with bitmaps
(i.e. some points on the input grid may be undefined).  In this case,
the bitmap is interpolated to the output grid.  Only valid input points
are used to interpolate to valid output points.  An output bitmap will also be
created to locate invalid data where the output grid extends outside the domain
of the input grid. 

The driver routine for interpolating scalars is ipolates, while the routine
for interpolating vectors is ipolatev.  The interpolation method is chosen
via the first argument of these routines (variable IP).  Sub-options are
set via the IPOPT array.

Bilinear interpolation is chosen by setting IP=0.  This method has two 
sub-options.  (1) The percent of valid input data required to make output
data (the default is 50%). (2) If valid input data is not found near an
a spiral search may be performed.  The spiral search is only
an option for scalar data.  The bilinear method also has no restrictions
and can interpolate with bitmaps.

Bicubic interpolation is chosen by setting IP=1.  This method has two sub-options,
(1) A monotonic constraint option for straight bicubic or
for constraining the output value to be within the range of the four
surrounding input values.  (2) The percent of valid input data
required to make output data, which defaults to 50%.  Note: the bicubic method
cannot interpolate data with bitmaps.

Neighbor interpolation is chosen by setting IP=2.  Neighbor interpolation
means that the output value is set to the nearest input value.  It would be
appropriate for interpolating integer fields such as vegetation index.
This method has one sub-option: If valid input data is not found near an 
an output point, a spiral search is optionally performed.  The neighbor
method has no restrictions and can interpolate with bitmaps.

Budget interpolation is chosen by setting IP=3.  Budget interpolation
means a low-order interpolation method that quasi-conserves area averages.
It would be appropriate for interpolating budget fields such as precipitation.
This method assumes that the field really represents box averages where each
box extends halfway to its neighboring grid point in each direction.
The method actually averages bilinearly interpolated values in a square array
of points distributed within each output grid box.  This method can
interpolate data with bitmaps.  There are several sub-options:
-  (1) The number of points in the radius of the square array may be set.
   The default is 2, meaning that 25 sample points will be averaged
   for each output value.  
-  (2) The respective averaging weights for the radius points are adjustable.
   The default is for all weights equal to 1, giving an unweighted average.
-  (3) Optionally, one may assume the boxes stretch nearly all the way to each of 
   the neighboring grid points and the weights are the adjoint of the bilinear 
   interpolation weights.
-  (4) The percent of valid input data required to make output data is
   adjustable.  The default is 50%.
-  (5) In cases where there is no or insufficient valid input data,
   a spiral search may be invoked to search for the nearest valid data.
   search square (scalar interpolation only). 

Spectral interpolation is chosen by setting IP=4.  This method has two sub-options,
to (1) set the spectral shape (triangular or rhomboidal) and (2) set the 
spectral truncation.  The input grid must be a global cylindrical grid
(either Gaussian or equidistant).  This method cannot interpolate data
with bitmaps. Unless the output grid is a global cylindrical
grid, a polar stereographic grid centered at the pole, or a Mercator grid,
this method can be quite expensive.

Neighbor-budget interpolation is chosen by setting IP=6.  This method
computes weighted averages of neighbor points arranged in a square box
centered around each output grid point and stretching nearly halfway
to each of the neighboring grid points. The main difference with 
the budget interpolation (IP=3) is neighbor vs bilinear interpolation
of the square box of points.  There are the following sub-options:
-  (1) The number of points in the radius of the square array may be set.
   The default is 2, meaning that 25 sample points will be averaged
   for each output value.  
-  (2) The respective averaging weights for the radius points are adjustable.
   The default is for all weights equal to 1, giving an unweighted average.
-  (3) The percent of valid input data required to make output data is
   adjustable.  The default is 50%.

The library can handle two-dimensional vector fields as well as scalar fields.
The input and output vectors are rotated if necessary so that they are
either resolved relative to their defined grid in the direction of
increasing x and y coordinates or resolved relative to eastward and northward
directions on the earth.  The rotation is determined by the grid definitions.
Vectors are generally interpolated (by all methods but spectral interpolation)
by moving the relevant input vectors along a great circle to the output point,
keeping their orientations with respect to the great circle constant, before
independently interpolating the respective components.  This ensures that vector
interpolation will be consistent over the whole globe including the poles.

The input and output grids are defined by their respective GRIB2 grid definition
template and template number as decoced by the NCEP G2 library.  There
are six map projections recognized by the library:

Grid def. template #    | Map projection
--------------------    | ---------------------------------
       00               | Equidistant cyclindrical
       01               | Rotated equidistant cylindrical
       10               | Mercator cyclindrical
       20               | Polar stereographic azimuthal
       30               | Lambert conformal conical
       40               | Gaussian equidistant cyclindrical

If the output grid definition template number is negative, then the
output data may be just a set of station points.  In this case, the user must pass
the number of points to be output along with their latitudes and longitudes.
For vector interpolation, the vector rotations parameters must also be passed.
On the other hand, for non-negative output data representation types,
the number of output grid points and their latitudes and longitudes
(and the vector rotation parameters for vector interpolation) are all
returned by the interpolation subprograms.

If an output equidistant cylindrical grid contains multiple pole points, then
the pole points are forced to be self-consistent.  That is, scalar fields
are obliged to be constant at the pole and vector components are obliged
to exhibit a wavenumber one variation at the pole.

Generally, only regular grids can be interpolated in this library.  However,
the thinned WAFS grids may be expanded to a regular grid (or vice versa)
using subprograms ipxwafs/2/3.  Eta data (with Arakawa "E" staggering)
on the "H" or "V" grid may be expanded to a filled regular grid (or vice versa)
using subprogram ipxetas.

The return code issued by an interpolation subprogram determines whether
it ran successfully or how it failed.  Check nonzero return codes
against the docblock of the respective subprogram.

Developers are encouraged to create additional interpolation methods or
to create additional map projection "wizards" for ip2lib.

Questions may be directed to: NCEP.List.EMC.nceplibs.Developers@noaa.gov

## Entry point list

Scalar field interpolation subprograms

   Name       |Function
   ----       |------------------------------------------------------------------
   IPOLATES   |IREDELL'S POLATE FOR SCALAR FIELDS
   POLATES0   |INTERPOLATE SCALAR FIELDS (BILINEAR)
   POLATES1   |INTERPOLATE SCALAR FIELDS (BICUBIC)
   POLATES2   |INTERPOLATE SCALAR FIELDS (NEIGHBOR)
   POLATES3   |INTERPOLATE SCALAR FIELDS (BUDGET)
   POLATES4   |INTERPOLATE SCALAR FIELDS (SPECTRAL)
   POLATES6   |INTERPOLATE SCALAR FIELDS (NEIGHBOR-BUDGET)
   POLFIXS    |MAKE MULTIPLE POLE SCALAR VALUES CONSISTENT

Vector field interpolation subprograms

   Name       |Function
   ----       |------------------------------------------------------------------
   IPOLATEV   |IREDELL'S POLATE FOR VECTOR FIELDS
   POLATEV0   |INTERPOLATE VECTOR FIELDS (BILINEAR)
   POLATEV1   |INTERPOLATE VECTOR FIELDS (BICUBIC)
   POLATEV2   |INTERPOLATE VECTOR FIELDS (NEIGHBOR)
   POLATEV3   |INTERPOLATE VECTOR FIELDS (BUDGET)
   POLATEV4   |INTERPOLATE VECTOR FIELDS (SPECTRAL)
   POLATEV6   |INTERPOLATE VECTOR FIELDS (NEIGHBOR-BUDGET)
   MOVECT     |MOVE A VECTOR ALONG A GREAT CIRCLE
   POLFIXV    |MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT

Grid description section decoders

   Name       |Function
   ----       |------------------------------------------------------------------
   GDSWZD                          |GRID DESCRIPTION SECTION WIZARD
   GDSWZD_C                        |'C' WRAPPER FOR CALLING GDSWZD
   GDSWZD_EQUID_CYLIND             |GDS WIZARD FOR EQUIDISTANT CYCLINDRICAL
   GDSWZD_MERCATOR                 |GDS WIZARD FOR MERCATOR CYCLINDRICAL
   GDSWZD_LAMBERT_CONF             |GDS WIZARD FOR LAMBERT CONFORMAL CONICAL
   GDSWZD_GAUSSIAN                 |GDS WIZARD FOR GAUSSIAN CYCLINDRICAL
   GDSWZD_POLAR_STEREO             |GDS WIZARD FOR POLAR STEREOGRAPHIC
   GDSWZD_ROT_EQUID_CYLIND_EGRID   |GDS WIZARD FOR ROTATED EQUIDISTANT CYCLINDRICAL "E" STAGGER.
   GDSWZD_ROT_EQUID_CYLIND         |GDS WIZARD FOR ROTATED EQUIDISTANT CYCLINDRICAL NON "E" STAGGER.
   IJKGDS0/1                       |RETURN FIELD POSITION FOR A GIVEN GRID POINT
   
Transform subprograms for special irregular grids

   Name       |Function
   ----       |------------------------------------------------------------------
   IPXWAFS/2/3 |  EXPAND OR CONTRACT WAFS GRIDS

## How to invoke ip2lib: examples

<pre>
***********************************************************************
Example 1.  Read a grib 2 file of scalar data on a global regular
            1-deg lat/lon grid and call ipolates to interpolate
            it to NCEP standard grid 218, a lambert conformal grid.
            Uses the NCEP G2 library to degrib the data.
***********************************************************************

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
! the output grid specs.  this is ncep grid 218, a lambert conformal
! grid.  the grid definition information is stored in section 3
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
! open the grib 2 file containing data to be interpolated.  for this
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

***********************************************************************
Example 2.  Read a grib 2 file of u/v wind data on a global regular
            1-deg lat/lon grid and call ipolatev to interpolate
            it to four random station points.  Uses the NCEP
            G2 library to degrib the data.
***********************************************************************

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
</pre>
