!> @file
!> @brief Determine earth radius and shape.
!> @author Gayno @date 2015-07-14

!> Determine earth radius and shape.
!>
!> Determine the radius and shape of the earth from
!> the grib 2 grid definition template array - section 3.
!>
!> @param[in] igdtmpl integer (igdtlen) grid definition template
!> array. Corresponds to the gfld%igdtmpl component of the ncep g2
!> library gridmod data structure. For all map projections recognized
!> by iplib, the entries use by this routine are:
!> - 1 shape of earth, section 3, octet 15 
!> - 2 scale factor of spherical earth radius, octet 16
!> - 3 scaled value of radius of spherical earth, octets 17-20
!> - 4 scale factor of major axis of elliptical earth, octet 21
!> - 5 scaled value of major axis of elliptical earth, octets 22-25
!> - 6 scale factor of minor axis of elliptical earth, octet 26
!> - 7 scaled value of minor axis of elliptical earth, octets 27-30
!> @param[in] igdtlen integer number of elements of the grid
!> definition template array. Corresponds to the gfld%igdtlen
!> component of the ncep g2 library gridmod data structure.
!> @param[out] radius real earth radius in meters. For ellipitical
!> earths, this is the semi major axis. See "map projectsions - a
!> working manual" by Snyder (1987) for details.
!> @param[out] eccen_squared real earth eccentricity squared
!>
!> @author Gayno @date 2015-07-14
SUBROUTINE EARTH_RADIUS(IGDTMPL, IGDTLEN, RADIUS, ECCEN_SQUARED)
  IMPLICIT NONE
  !
  INTEGER,                INTENT(IN   ) :: IGDTLEN
  INTEGER,                INTENT(IN   ) :: IGDTMPL(IGDTLEN)
  !
  REAL,                   INTENT(  OUT) :: ECCEN_SQUARED
  REAL,                   INTENT(  OUT) :: RADIUS
  !
  REAL                                  :: FLAT
  REAL                                  :: MAJOR_AXIS, MINOR_AXIS
  !
  SELECT CASE (IGDTMPL(1))
  CASE (0)
     RADIUS        = 6367470.0
     ECCEN_SQUARED = 0.0
  CASE (1)  ! USER SPECIFIED SPHERICAL
     RADIUS        = FLOAT(IGDTMPL(3))/FLOAT(10**IGDTMPL(2))
     ECCEN_SQUARED = 0.0
  CASE (2)  ! IAU 1965
     RADIUS        = 6378160.0      ! SEMI MAJOR AXIS
     FLAT          = 1.0/297.0      ! FLATTENING
     ECCEN_SQUARED = (2.0*FLAT) - (FLAT**2)
  CASE (3)  ! USER SPECIFIED ELLIPTICAL (KM)
     MAJOR_AXIS    = FLOAT(IGDTMPL(5))/FLOAT(10**IGDTMPL(4))
     MAJOR_AXIS    = MAJOR_AXIS * 1000.0
     MINOR_AXIS    = FLOAT(IGDTMPL(7))/FLOAT(10**IGDTMPL(6))
     MINOR_AXIS    = MINOR_AXIS * 1000.0
     ECCEN_SQUARED = 1.0 - (MINOR_AXIS**2 / MAJOR_AXIS**2)
     RADIUS        = MAJOR_AXIS
  CASE (4)  ! IAG-GRS80 MODEL
     RADIUS        = 6378137.0      ! SEMI MAJOR AXIS
     FLAT          = 1.0/298.2572   ! FLATTENING
     ECCEN_SQUARED = (2.0*FLAT) - (FLAT**2)
  CASE (5)  ! WGS84 DATUM
     RADIUS        = 6378137.0      ! SEMI MAJOR AXIS
     ECCEN_SQUARED = 0.00669437999013
  CASE (6)
     RADIUS        = 6371229.0
     ECCEN_SQUARED = 0.0
  CASE (7)  ! USER SPECIFIED ELLIPTICAL (M)
     MAJOR_AXIS    = FLOAT(IGDTMPL(5))/FLOAT(10**IGDTMPL(4))
     MINOR_AXIS    = FLOAT(IGDTMPL(7))/FLOAT(10**IGDTMPL(6))
     ECCEN_SQUARED = 1.0 - (MINOR_AXIS**2 / MAJOR_AXIS**2)
     RADIUS        = MAJOR_AXIS
  CASE (8)
     RADIUS        = 6371200.0
     ECCEN_SQUARED = 0.0
  CASE DEFAULT
     RADIUS        = -9999.
     ECCEN_SQUARED = -9999.
  END SELECT
  !
  RETURN
  !
END SUBROUTINE EARTH_RADIUS
