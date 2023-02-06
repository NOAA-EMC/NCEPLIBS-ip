!> @file
!> @brief C wrapper for routine gdswzd().
!> @author JOVIC @date 2016-04-10

!> @brief Module that contains C wrapper for routine gdswzd().
!> @author JOVIC @date 2016-04-10
module gdswzd_c_mod
  use gdswzd_mod
  implicit none

contains

  !> C wrapper for routine gdswzd.
  !>
  !> Use this routine to call 'gdswzd' from a C or C++ program. Takes
  !> advantage of the fortran/c interoperability standard.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2016-04-10 | Jovic | initial version.
  !> 2016-04-26 | G. Gayno | update for grib 2
  !>
  !> ### Example: Mixed precsion iplib call
  !>
  !> @code
  !>  include "iplib.h"
  !>
  !>  int *igdtmpl;
  !>  int igdtnum, igdtlen, iopt, npts, nret;
  !>  double fill;
  !>  double *xpts, *ypts, *rlon, *rlat;
  !>  double *crot, *srot, *xlon, *xlat, *ylon, *ylat, *area;
  !>
  !>  gdswzd(igdtnum, igdtmpl, igdtlen, iopt, npts,
  !>         fill, xpts, ypts, rlon, rlat,
  !>         &nret, crot, srot, xlon,
  !>         xlat, ylon, ylat, area);
  !> @endcode
  !>
  !> @param[in] igdtnum grid definition template number. Corresponds
  !> to the gfld%igdtnum component of the ncep g2 library gridmod data
  !> structure.
  !> - 00 equidistant cylindrical
  !> - 01 rotated equidistant cylindrical. "e" and non-"e" staggered
  !> - 10 mercator cyclindrical
  !> - 20 polar stereographic azimuthal
  !> - 30 lambert conformal conical
  !> - 40 gaussian equidistant cyclindrical
  !> @param[in] igdtmpl (igdtlen) grid definition template array.
  !> corresponds to the gfld%igdtmpl component of the ncep g2 library
  !> gridmod data structure for section three.
  !> all projections:
  !> - 1 shape of earth, octet 15
  !> - 2 scale factor of spherical earth radius, octet 16
  !> - 3 scaled value of radius of spherical earth, octets 17-20
  !> - 4 scale factor of major axis of elliptical earth, octet 21
  !> - 5 scaled value of major axis of elliptical earth, octets 22-25
  !> - 6 scale factor of minor axis of elliptical earth, octet 26
  !> - 7 scaled value of minor axis of elliptical earth, octets 27-30
  !> equidistant cyclindrical:
  !> - 8 number of points along a parallel, octs 31-34
  !> - 9 number of points along a meridian, octs 35-38
  !> - 10 basic angle of initial production domain, octets 39-42.
  !> - 11 subdivisions of basic angle, octets 43-46
  !> - 12 latitude of first grid point, octets 47-50
  !> - 13 longitude of first grid point, octets 51-54
  !> - 14 resolution and component flags, octet 55
  !> - 15 latitude of last grid point, octets 56-59
  !> - 16 longitude of last grid point, octets 60-63
  !> - 17 i-direction increment, octets 64-67
  !> - 18 j-direction increment, octets 68-71
  !> - 19 scanning mode, octet 72
  !> mercator cyclindrical:
  !> - 8  number of points along a parallel, octs 31-34
  !> - 9  number of points along a meridian, octs 35-38
  !> - 10 latitude of first point, octets 39-42
  !> - 11 longitude of first point, octets 43-46
  !> - 12 resolution and component flags, octet 47
  !> - 13 tangent latitude, octets 48-51
  !> - 14 latitude of last point, octets 52-55
  !> - 15 longitude of last point, octets 56-59
  !> - 16 scanning mode flags, octet 60
  !> - 17 orientation of grid, octets 61-64
  !> - 18 longitudinal grid length, octets 65-68
  !> - 19 latitudinal grid length, octets 69-72
  !> Lambert conformal conical:
  !> - 8  number of points along x-axis, octs 31-34
  !> - 9  number of points along y-axis, octs 35-38
  !> - 10 latitude of first point, octets 39-42
  !> - 11 longitude of first point, octets 43-46
  !> - 12 resolution of component flag, octet 47
  !> - 13 latitude where grid lengths specified, octets 48-51
  !> - 14 longitude of meridian that is parallel to y-axis, octets 52-55
  !> - 15 x-direction grid length, octets 56-59
  !> - 16 y-direction grid length, octets 60-63
  !> - 17 projection center flag, octet 64
  !> - 18 scanning mode, octet 65
  !> - 19 first tangent latitude from pole, octets 66-69
  !> - 20 second tangent latitude from pole, octets 70-73
  !> - 21 latitude of south pole of projection, octets 74-77
  !> - 22 longitude of south pole of projection, octets 78-81
  !> gaussian cylindrical:
  !> - 8  number of points along a parallel, octs 31-34
  !> - 9  number of points along a meridian, octs 35-38
  !> - 10 basic angle of initial production domain, octets 39-42
  !> - 11 subdivisions of basic angle, octets 43-46
  !> - 12 latitude of first grid point, octets 47-50
  !> - 13 longitude of first grid point, octets 51-54
  !> - 14 resolution and component flags, octet 55
  !> - 15 latitude of last grid point, octets 56-59
  !> - 16 longitude of last grid point, octets 60-63
  !> - 17 i-direction increment, octets 64-67
  !> - 18 number of parallels between pole and equator, octets 68-71
  !> - 19 scanning mode, octet 72
  !> polar stereographic azimuthal:
  !> - 8  number of points along x-axis, octets 31-34
  !> - 9  number of points along y-axis, octets 35-38
  !> - 10 latitude of first grid point, octets 39-42
  !> - 11 longitude of first grid point, octets 43-46
  !> - 12 resolution and component flags, octet 47
  !> - 13 true latitude, octets 48-51
  !> - 14 orientation longitude, octets 52-55
  !> - 15 x-direction grid length, octets 56-59
  !> - 16 y-direction grid length, octets 60-63
  !> - 17 projection center flag, octet 64
  !> - 18 scanning mode flags, octet 65
  !> rotated equidistant cyclindrical:
  !> - 8  number of points along a parallel, octs 31-34
  !> - 9  number of points along a meridian, octs 35-38
  !> - 10 basic angle of initial production domain, octets 39-42
  !> - 11 subdivisions of basic angle, octets 43-46
  !> - 12 latitude of first grid point, octets 47-50
  !> - 13 longitude of first grid point, octets 51-54
  !> - 14 resolution and component flags, octet 55
  !> - 15 latitude of last grid point, octets 56-59
  !> - 16 longitude of last grid point, octets 60-63
  !> - 17 i-direction increment, octets 64-67
  !> - 18 j-direction increment, octets 68-71
  !> - 19 scanning mode, octet 72
  !> - 20 latitude of southern pole of projection, octets 73-76
  !> - 21 longitude of southern pole of projection, octets 77-80
  !> - 22 angle of rotation of projection, octs 81-84
  !> @param[in] igdtlen number of elements of the grid definition
  !> template array. Corresponds to the gfld%igdtlen component of the
  !> ncep g2 library gridmod data structure.
  !> @param[in] iopt option flag
  !> - 0 to compute earth coords of all the grid points
  !> - 1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param[in] npts integer maximum number of coordinates
  !> @param[in] fill real fill value to set invalid output data (must
  !> be impossible value; suggested value: -9999.)
  !> @param[inout] xpts real (npts) grid x point coordinates if iopt>0
  !> @param[inout] ypts real (npts) grid y point coordinates if iopt>0
  !> @param[inout] rlon real (npts) earth longitudes in degrees e if
  !> iopt<0 (acceptable range: -360. to 360.)
  !> @param[inout] rlat real (npts) earth latitudes in degrees n if
  !> iopt<0 (acceptable range: -90. to 90.)
  !> @param[out] nret number of valid points computed (-1 if
  !> projection unrecognized)
  !> @param[out] crot (npts) clockwise vector rotation cosines
  !> @param[out] srot (npts) clockwise vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth, vgrid=srot*uearth+crot*vearth)
  !> @param[out] xlon (npts) dx/dlon in 1/degrees
  !> @param[out] xlat (npts) dx/dlat in 1/degrees
  !> @param[out] ylon (npts) dy/dlon in 1/degrees
  !> @param[out] ylat (npts) dy/dlat in 1/degrees
  !> @param[out] AREA (npts) area weights in m**2 (Proportional to the
  !>
  !> @author JOVIC @date 2016-04-10
  SUBROUTINE GDSWZD_C(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET,CROT,SROT, &
       XLON,XLAT,YLON,YLAT,AREA) BIND(C, NAME='gdswzd')
    USE, INTRINSIC :: ISO_C_BINDING

    USE GDSWZD_MOD

    IMPLICIT NONE

    INTEGER(KIND=C_INT),  INTENT(IN)         :: IGDTMPL(IGDTLEN)
    INTEGER(KIND=C_INT),  VALUE, INTENT(IN)  :: IGDTNUM, IGDTLEN
    INTEGER(KIND=C_INT),  VALUE, INTENT(IN)  :: IOPT, NPTS
    INTEGER(KIND=C_INT),  INTENT(OUT)        :: NRET

#if (LSIZE==D)
    REAL(KIND=C_DOUBLE), VALUE, INTENT(IN)   :: FILL
    REAL(KIND=C_DOUBLE), INTENT(INOUT)       :: XPTS(NPTS),YPTS(NPTS),RLON(NPTS),RLAT(NPTS)
    REAL(KIND=C_DOUBLE), INTENT(OUT)         :: CROT(NPTS),SROT(NPTS),XLON(NPTS),XLAT(NPTS)
    REAL(KIND=C_DOUBLE), INTENT(OUT)         :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
#else
    REAL(KIND=C_FLOAT), VALUE, INTENT(IN)    :: FILL
    REAL(KIND=C_FLOAT), INTENT(INOUT)        :: XPTS(NPTS),YPTS(NPTS),RLON(NPTS),RLAT(NPTS)
    REAL(KIND=C_FLOAT), INTENT(OUT)          :: CROT(NPTS),SROT(NPTS),XLON(NPTS),XLAT(NPTS)
    REAL(KIND=C_FLOAT), INTENT(OUT)          :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
#endif

    CALL GDSWZD(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
         XPTS,YPTS,RLON,RLAT,NRET, &
         CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

  END SUBROUTINE GDSWZD_C

  !> C wrapper for routine gdswzd.
  !> Use this routine to call 'gdswzd' from a C or C++ program.
  !> Takes advantage of the fortran/c interoperability standard.
  !>
  !> ### Example: mixed precsion iplib call
  !>@code
  !>  include "iplib.h"
  !>
  !>  int kgds[200];
  !>  int iopt, npts, nret;
  !>  double fill;
  !>  double *xpts, *ypts, *rlon, *rlat;
  !>  double *crot, *srot, *xlon, *xlat, *ylon, *ylat, *area;
  !>
  !>  gdswzd(kgds, iopt, npts, fill,
  !>         xpts, ypts, rlon, rlat,
  !>         &nret, crot, srot, xlon,
  !>         xlat, ylon, ylat, area);
  !> @endcode
  !>
  !> @param kgds (200) gds parameters as decoded by w3fi63
  !> @param iopt option flag
  !> - 0 to compute earth coords of all the grid points
  !> - 1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param npts maximum number of coordinates
  !> @param fill fill value to set invalid output data (must be
  !> impossible value; suggested value: -9999.)
  !> @param xpts (npts) grid x point coordinates if iopt>0
  !> @param ypts (npts) grid y point coordinates if iopt>0
  !> @param[out] rlon (npts) earth longitudes in degrees e if iopt<0
  !> (acceptable range: -360. to 360.)
  !> @param[out] rlat (npts) earth latitudes in degrees n if iopt<0
  !> (acceptable range: -90. to 90.)
  !> @param[out] nret number of valid points computed (-1 if
  !> projection unrecognized)
  !> @param[out] crot (npts) clockwise vector rotation cosines
  !> @param[out] srot (npts) clockwise vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] xlon (npts) dx/dlon in 1/degrees
  !> @param[out] xlat (npts) dx/dlat in 1/degrees
  !> @param[out] ylon (npts) dy/dlon in 1/degrees
  !> @param[out] ylat (npts) dy/dlat in 1/degrees
  !> @param[out] area (npts) area weights in m**2 (proportional to the
  !> square of the map factor in the case of conformal projections.)
  !>
  !>  @author JOVIC @date 2016-04-10
  SUBROUTINE GDSWZD_C_grib1(KGDS,IOPT,NPTS,FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA) BIND(C, NAME='gdswzd_grib1')
    USE, INTRINSIC :: ISO_C_BINDING

    USE GDSWZD_MOD

    IMPLICIT NONE

    INTEGER(KIND=C_INT),  INTENT(IN)         :: KGDS(200)
    INTEGER(KIND=C_INT),  VALUE, INTENT(IN)  :: IOPT
    INTEGER(KIND=C_INT),  VALUE, INTENT(IN)  :: NPTS
    INTEGER(KIND=C_INT),  INTENT(OUT)        :: NRET

#if (LSIZE==D)
    REAL(KIND=C_DOUBLE), VALUE, INTENT(IN)   :: FILL
    REAL(KIND=C_DOUBLE), INTENT(INOUT)       :: XPTS(NPTS),YPTS(NPTS),RLON(NPTS),RLAT(NPTS)
    REAL(KIND=C_DOUBLE), INTENT(OUT)         :: CROT(NPTS),SROT(NPTS),XLON(NPTS),XLAT(NPTS)
    REAL(KIND=C_DOUBLE), INTENT(OUT)         :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
#else
    REAL(KIND=C_FLOAT), VALUE, INTENT(IN)    :: FILL
    REAL(KIND=C_FLOAT), INTENT(INOUT)        :: XPTS(NPTS),YPTS(NPTS),RLON(NPTS),RLAT(NPTS)
    REAL(KIND=C_FLOAT), INTENT(OUT)          :: CROT(NPTS),SROT(NPTS),XLON(NPTS),XLAT(NPTS)
    REAL(KIND=C_FLOAT), INTENT(OUT)          :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
#endif

    CALL GDSWZD(KGDS,IOPT,NPTS,FILL,XPTS,YPTS,RLON,RLAT,NRET, &
         CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

  END SUBROUTINE GDSWZD_C_GRIB1

end module gdswzd_c_mod
