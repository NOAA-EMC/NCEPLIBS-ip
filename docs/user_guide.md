@mainpage

## Documentation of the general interpolation library iplib

January, 2014

The NCEP general interpolation library (iplib) contains Fortran 90 subprograms
to be used for interpolating between nearly all grids used at NCEP.
The library is particularly efficient when interpolating many fields at one time.
The library has been extensively tested with AIX and Intel Fortran compilers.

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
for interpolating vectors is ipolatv.  The interpolation method is chosen
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
  (1) The number of points in the radius of the square array may be set.
   The default is 2, meaning that 25 sample points will be averaged
   for each output value.  
  (2) The respective averaging weights for the radius points are adjustable.
   The default is for all weights equal to 1, giving an unweighted average.
  (3) Optionally, one may assume the boxes stretch nearly all the way to each of 
   the neighboring grid points and the weights are the adjoint of the bilinear 
   interpolation weights.
  (4) The percent of valid input data required to make output data is
   adjustable.  The default is 50%.
  (5) In cases where there is no or insufficient valid input data,
   a spiral search may be invoked to search for the nearest valid data.
   search squre (scalar interpolation only). 

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
  (1) The number of points in the radius of the square array may be set.
   The default is 2, meaning that 25 sample points will be averaged
   for each output value.  
  (2) The respective averaging weights for the radius points are adjustable.
   The default is for all weights equal to 1, giving an unweighted average.
  (3) The percent of valid input data required to make output data is
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

The input and output grids are defined by their respective GRIB grid description
sections (GDS) passed in integer form KGDS as decoded by subprogram w3fi63
in w3ncolib or by subprogram makgds in this library.  That is, the interpolation
subprograms can readily interpolate from a GRIB field that is unpacked
by subprogram w3fi63; the interpolation subprograms can also readily
interpolate to an NCEP pre-defined grid that is expanded into KGDS form by 
subprogram makgds (which in turn calls w3fi71).  There are currently seven grid
projections recognized in the library.  The projections are respectively
equidistant cylindrical (KGDS(1)=000), Mercator cylindrical (KGDS(1)=001),
Lambert conformal conical (KGDS(1)=003), Gaussian cylindrical (KGDS(1)=004), 
polar stereographic azimuthal (KGDS(1)=005), semi-staggered "E" grid 
on rotated equidistant cyclindrical (KGDS(1)=203), and non-staggered
grid on rotated equidistant cyclindrical grid (KGDS(1)=205).

If the output data representation type is negative (KGDSO(1)<0), then the
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
the thinned WAFS grids and the staggered eta grids can be interpolated by using
transform subprograms (ipxwafs and ipxetas, respectively) in this library that
will either expand the irregular grid to a regular grid or contract a regular
grid to an irregular grid as necessary.

The return code issued by an interpolation subprogram determines whether
it ran successfully or how it failed.  Check nonzero return codes
against the docblock of the respective subprogram.

Developers are encouraged to create additional interpolation methods or
to create additional map projection "wizards" for iplib.

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
   GDSWZD     |GRID DESCRIPTION SECTION WIZARD
   GDSWZD_C   |'C' WRAPPER FOR CALLING GDSWZD
   GDSWZD00   |GDS WIZARD FOR EQUIDISTANT CYCLINDRICAL
   GDSWZD01   |GDS WIZARD FOR MERCATOR CYCLINDRICAL
   GDSWZD03   |GDS WIZARD FOR LAMBERT CONFORMAL CONICAL
   GDSWZD04   |GDS WIZARD FOR GAUSSIAN CYCLINDRICAL
   GDSWZD05   |GDS WIZARD FOR POLAR STEREOGRAPHIC
   GDSWZDCB   |GDS WIZARD FOR ROTATED EQUIDISTANT CYCLINDRICAL "E" STAGGER.
   GDSWZDCD   |GDS WIZARD FOR ROTATED EQUIDISTANT CYCLINDRICAL NON "E" STAGGER.
   GAUSSLAT   |COMPUTE GAUSSIAN LATITUDES
   IJKGDS0/1  |RETURN FIELD POSITION FOR A GIVEN GRID POINT
   MAKGDS     |MAKE OR BREAK A GRID DESCRIPTION SECTION

Transform subprograms for special irregular grids

   Name        |Function
   ----        |------------------------------------------------------------------
   IPXWAFS/2/3 |  EXPAND OR CONTRACT WAFS GRIDS

## How to inoke iplib: examples

### Example 1

Interpolate from an arbitrary input grid (probably 1x1) to NCEP grid
27 (65x65 northern polar stereographic).  Interpolate heights
bilinearly and winds bicubically.  Interpolate soil moisture and
precipitation using bitmaps with the budget method.  Encode the soil
moisture bitmap in GRIB.  Subprograms GETGB and PUTGB from w3ncolib
are referenced.

<pre>
c example of using ipolate package.
c see documentation of ipolates and ipolatev
c for further possible options.
      integer ipopt(20)
      integer jpds(25),jgds(22),kpdsi(25),kgdsi(22),kpdso(25),kgdso(22)
      parameter(ji=360*181,ig=27,jo=65*65,km=4)
      real rlat(jo),rlon(jo),crot(jo),srot(jo)
      integer ibi(km),ibo(km)
      logical li(ji,km),lo(jo,km)
      real hi(ji,km),ri(ji),ui(ji),vi(ji)
      real ho(jo,km),ro(jo),uo(jo),vo(jo)
      character gdso(42)
      integer lev(km)
      data lev/1000,500,250,100/
c define 65x65 grid
      call makgds(ig,kgdso,gdso,lengds,iret)
      if(iret.ne.0) call exit(iret)
       kgdso(4)=-20826! fix w3fi71 error
      print *,'kgdso=',kgdso
      ipopt=0
c interpolate 4 levels of height to 65x65 bilinearly
      do k=1,km
        jpds=-1
        jpds(5)=7
        jpds(6)=100
        jpds(7)=lev(k)
        call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &             li(1,k),hi(1,k),iret)
        if(iret.ne.0) call exit(iret)
        call putgb(50,ki,kpdsi,kgdsi,li(1,k),hi(1,k),iret)
        if(iret.ne.0) call exit(iret)
        ibi(k)=mod(kpdsi(4)/64,2)
        print *,'ibi(k)=',ibi(k)
      enddo
      call ipolates(0,ipopt,kgdsi,kgdso,ji,jo,km,ibi,li,hi,
     &              ko,rlat,rlon,ibo,lo,ho,iret)
      if(iret.ne.0) call exit(iret)
      kpdso=kpdsi
      kpdso(3)=ig
      do k=1,km
        kpdso(7)=lev(k)
        call putgb(51,ko,kpdso,kgdso,lo(1,k),ho(1,k),iret)
        if(iret.ne.0) call exit(iret)
      enddo
c interpolate precipitation to 65x65 with budget method
c (zero rain is masked out)
      jpds=-1
      jpds(5)=61
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,ri,iret)
      if(iret.ne.0) call exit(iret)
      call putgb(50,ki,kpdsi,kgdsi,li,ri,iret)
      if(iret.ne.0) call exit(iret)
      ipopt(1)=-1
      ipopt(2)=-1
      li(1:ki,1)=ri(1:ki).gt.0.
      call ipolates(3,ipopt,kgdsi,kgdso,ji,jo,1,1,li,ri,
     &              ko,rlat,rlon,ibo,lo,ro,iret)
      if(iret.ne.0) call exit(iret)
      kpdso=kpdsi
      kpdso(3)=ig
      call putgb(51,ko,kpdso,kgdso,lo,ro,iret)
      if(iret.ne.0) call exit(iret)
c interpolate soil moisture to 65x65 with budget method
      jpds=-1
      jpds(5)=144
      jpds(6)=112
      jpds(7)=10
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,ri,iret)
      if(iret.ne.0) call exit(iret)
      call putgb(50,ki,kpdsi,kgdsi,li,ri,iret)
      if(iret.ne.0) call exit(iret)
      ibi(1)=mod(kpdsi(4)/64,2)
      ipopt(1)=-1
      ipopt(2)=-1
      call ipolates(3,ipopt,kgdsi,kgdso,ji,jo,1,ibi,li,ri,
     &              ko,rlat,rlon,ibo,lo,ro,iret)
      if(iret.ne.0) call exit(iret)
      kpdso=kpdsi
      kpdso(3)=ig
      kpdso(4)=128+64*ibo(1)
      call putgb(51,ko,kpdso,kgdso,lo,ro,iret)
      if(iret.ne.0) call exit(iret)
c interpolate 200 mb winds to 65x65 bicubically
      jpds=-1
      jpds(5)=33
      jpds(6)=100
      jpds(7)=200
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,ui,iret)
      if(iret.ne.0) call exit(iret)
      call putgb(50,ki,kpdsi,kgdsi,li,ui,iret)
      if(iret.ne.0) call exit(iret)
      jpds(5)=34
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,vi,iret)
      if(iret.ne.0) call exit(iret)
      call putgb(50,ki,kpdsi,kgdsi,li,vi,iret)
      if(iret.ne.0) call exit(iret)
      ipopt(1)=0
      call ipolatev(1,ipopt,kgdsi,kgdso,ji,jo,1,0,li,ui,vi,
     &              ko,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
      if(iret.ne.0) call exit(iret)
      kpdso=kpdsi
      kpdso(3)=ig
      kpdso(5)=33
      call putgb(51,ko,kpdso,kgdso,lo,uo,iret)
      if(iret.ne.0) call exit(iret)
      kpdso(5)=34
      call putgb(51,ko,kpdso,kgdso,lo,vo,iret)
      if(iret.ne.0) call exit(iret)
      stop
      end
</pre>

### Example 2

Interpolate winds from an arbitrary input grid (probably 1x1) to 4
station points while truncating spectrally to R30.

<pre>
c read and unpack the 500 mb winds, truncate to R30,
c and interpolate to 4 corners of a box
      integer ipopt(20)
      integer jpds(25),jgds(22),kpdsi(25),kgdsi(22),kgdso(22)
      parameter(jf=360*181,kp=4)
      real rlat(kp),rlon(kp),crot(kp),srot(kp)
      logical lgi(jf),lgo(kp)
      real ui(jf),vi(jf),uo(kp),vo(kp)
      jpds=-1
      jpds(5)=33
      jpds(6)=100
      jpds(7)=500
      call getgb(11,31,jf,0,jpds,jgds,kf,kr,kpdsi,kgdsi,
     &           lgi,ui,iret)
      jpds(5)=34
      call getgb(11,31,jf,0,jpds,jgds,kf,kr,kpdsi,kgdsi,
     &           lgi,vi,iret)
      kgdso=-1
      rlat(1)=20.
      rlat(2)=20.
      rlat(3)=10.
      rlat(4)=10.
      rlon(1)=-50.
      rlon(2)=-40.
      rlon(3)=-50.
      rlon(4)=-40.
      crot=1.
      srot=0.
      ipopt(1)=1
      ipopt(2)=30
      uo=-999
      vo=-999
      call ipolatev(4,ipopt,kgdsi,kgdso,jf,kp,1,0,lgi,ui,vi,
     &              kp,rlat,rlon,crot,srot,ibo,lgo,uo,vo,iret)
      print '(2(2x,2f8.2))',(uo(k),vo(k),k=1,kp)
      end
</pre>

### Example 3

Interpolate winds from an arbitrary input grid (probably 1x1)
bilinearly to 3 station points.

<pre>
c read and unpack 4 levels of heights and winds
c and interpolate to 3 sonde sites.
      integer ipopt(20)
      integer jpds(25),jgds(22),kpdsi(25),kgdsi(22),kgdso(22)
      parameter(ji=360*181,km=4,jo=3)
      real rlat(jo),rlon(jo),crot(jo),srot(jo)
      integer ibi(km),ibo(km)
      logical li(ji,km),lo(jo,km)
      real hi(ji,km),ui(ji,km),vi(ji,km),ho(jo,km),uo(jo,km),vo(jo,km)
      integer lev(km)
      data lev/1000,500,250,100/
c define output locations
      kgdso=-1
      rlat(1)=22.2
      rlat(2)=33.3
      rlat(3)=44.4
      rlon(1)=-50.
      rlon(2)=-40.
      rlon(3)=-30.
      crot=1.
      srot=0.
c heights
      do k=1,km
        jpds=-1
        jpds(5)=7
        jpds(6)=100
        jpds(7)=lev(k)
        call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &             li(1,k),hi(1,k),iret)
        if(iret.ne.0) call exit(iret)
        ibi(k)=mod(kpdsi(4)/64,2)
      enddo
      call ipolates(0,ipopt,kgdsi,kgdso,ji,jo,km,ibi,li,hi,
     &              jo,rlat,rlon,ibo,lo,ho,iret)
      if(iret.ne.0) call exit(iret)
c winds
      do k=1,km
        jpds=-1
        jpds(5)=33
        jpds(6)=100
        jpds(7)=lev(k)
        call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &             li(1,k),ui(1,k),iret)
        if(iret.ne.0) call exit(iret)
        jpds=-1
        jpds(5)=34
        jpds(6)=100
        jpds(7)=lev(k)
        call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &             li(1,k),vi(1,k),iret)
        if(iret.ne.0) call exit(iret)
        ibi(k)=mod(kpdsi(4)/64,2)
      enddo
      call ipolatev(0,ipopt,kgdsi,kgdso,ji,jo,km,ibi,li,ui,vi,
     &              jo,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
      if(iret.ne.0) call exit(iret)
      print '((i8,3(2x,3f8.2)))',
     & (lev(k),(ho(j,k),uo(j,k),vo(j,k),j=1,jo),k=1,km)
      end
</pre>

Example 4. Interpolate 850 mb height and winds from the staggered meso-eta
           to a regional 0.25 degree grid.
<pre>
      integer ipopt(20)
      integer jpds(25),jgds(22),kpdsi(25),kgdsi(22),kpdso(25),kgdso(22)
      integer kgdsi2(22)
      parameter(ji=361*271,ig=255,jo=121*81)
      real rlat(jo),rlon(jo),crot(jo),srot(jo)
      logical li(ji),lo(jo)
      real fi(ji),fi2(ji),fo(jo)
      real vi(ji),vi2(ji),vo(jo)
      character gdso(400)
      kgdso=0
      kgdso(1)=0
      kgdso(2)=121
      kgdso(3)=81
      kgdso(4)=30000
      kgdso(5)=-90000
      kgdso(6)=128
      kgdso(7)=50000
      kgdso(8)=-60000
      kgdso(9)=250
      kgdso(10)=250
      kgdso(11)=64
      kgdso(19)=0
      kgdso(20)=255
      if(iret.ne.0) call exit(iret)
      ipopt=0
      jpds=-1
      jpds(6)=100
      jpds(7)=850
      jpds(5)=7
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,fi,iret)
      if(iret.ne.0) call exit(iret)
      call ipxetas(1,ji,ji,1,kgdsi,fi,kgdsi2,fi2,iret)
      if(iret.ne.0) call exit(iret)
      call ipolates(0,ipopt,kgdsi2,kgdso,ji,jo,1,0,li,fi2,
     &              ko,rlat,rlon,ibo,lo,fo,iret)
      if(iret.ne.0) call exit(iret)
      kpdso=kpdsi
      kpdso(3)=ig
      kpdso(4)=128+64*ibo
      kpdso(22)=1
      call putgb(51,ko,kpdso,kgdso,lo,fo,iret)
      if(iret.ne.0) call exit(iret)
      jpds(5)=33
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,fi,iret)
      if(iret.ne.0) call exit(iret)
      call ipxetas(2,ji,ji,1,kgdsi,fi,kgdsi2,fi2,iret)
      if(iret.ne.0) call exit(iret)
      jpds(5)=34
      call getgb(11,31,ji,0,jpds,jgds,ki,kr,kpdsi,kgdsi,
     &           li,vi,iret)
      if(iret.ne.0) call exit(iret)
      call ipxetas(2,ji,ji,1,kgdsi,vi,kgdsi2,vi2,iret)
      if(iret.ne.0) call exit(iret)
      call ipolatev(0,ipopt,kgdsi2,kgdso,ji,jo,1,0,li,fi2,vi2,
     &              ko,rlat,rlon,crot,srot,ibo,lo,fo,vo,iret)
      if(iret.ne.0) call exit(iret)
      kpdso=kpdsi
      kpdso(3)=ig
      kpdso(4)=128+64*ibo
      kpdso(22)=1
      kpdso(5)=33
      call putgb(51,ko,kpdso,kgdso,lo,fo,iret)
      if(iret.ne.0) call exit(iret)
      kpdso(5)=34
      call putgb(51,ko,kpdso,kgdso,lo,vo,iret)
      if(iret.ne.0) call exit(iret)
      end
</pre>
