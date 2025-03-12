!> @file
!> @brief GDS wizard for mercator cylindrical.
!>
!> @author Iredell @date 96-04-10

!> @brief GDS wizard for mercator cylindrical.
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.10 -
!> Mercator](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-10.shtml).
!>
!> @author Iredell @date 96-04-10
module ip_mercator_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use ip_constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_mercator_grid

  type, extends(ip_grid) :: ip_mercator_grid
     real :: rlat1 !< Latitude of first grid point. Section 3, octets 39-42.
     real :: rlon1 !< Longitude of first grid point. Section 3, octets 43-46.
     real :: rlon2 !< Longitude of last grid point. Section 3, octets 56-59.
     real :: rlati !< Latitude at which the Mercator projection intersects the Earth. Section 3, octets 48-51.
     real :: hi !< Scan mode in the 'i' direction. Section 3, octet 60.
     real :: dlon !< Longitudinal direction grid length. Section 3, octets 65-68.
     real :: dphi !< Latitudinal direction grid length. Section 3, octets 69-72.
   contains
     !> Initializes a gaussian grid given a grib1_descriptor object.
     !> @return N/A @copydoc ip_mercator_grid_mod::init_grib1
     procedure :: init_grib1
     !> Initializes a gaussian grid given a grib2_descriptor object.
     !> @return N/A @copydoc ip_mercator_grid_mod::init_grib2
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
     !> for Gaussian grids.
     !> @return N/A @copydoc ip_mercator_grid_mod::gdswzd_mercator
     procedure :: gdswzd => gdswzd_mercator
  end type ip_mercator_grid

  REAL :: DLON !< Longitudinal direction grid length.
  REAL :: DPHI !< Latitudinal direction grid length.
  REAL :: RERTH !< Radius of the Earth.

CONTAINS

  !> Initializes a mercator grid given a grib1_descriptor object.
  !>
  !> @param[inout] self ip_mercator_grid object.
  !> @param[in] g1_desc GRIB1 descriptor.
  !>
  !> @author Iredell @date 96-04-10  
  subroutine init_grib1(self, g1_desc)
    class(ip_mercator_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan, jscan
    real :: dy, hj

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.0

      self%IM=KGDS(2)
      self%JM=KGDS(3)

      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      self%RLATI=KGDS(9)*1.E-3

      ISCAN=MOD(KGDS(11)/128,2)
      JSCAN=MOD(KGDS(11)/64,2)

      DY=KGDS(13)
      self%HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DPHI=HJ*DY/(self%RERTH*COS(self%RLATI/DPR))

      ! defaults
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0

      self%iwrap = nint(360 / abs(self%dlon))
      if (self%im < self%iwrap) self%iwrap = 0
    end associate

  end subroutine init_grib1

  !> Init GRIB2.
  !>
  !> @param[inout] self ip_mercator_grid object.
  !> @param[in] g2_desc GRIB2 descriptor.
  !>
  !> @author Iredell @date 96-04-10  
  subroutine init_grib2(self, g2_desc)
    class(ip_mercator_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscan, jscan
    real :: hj, dy

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%RLAT1=FLOAT(IGDTMPL(10))*1.0E-6
      self%RLON1=FLOAT(IGDTMPL(11))*1.0E-6
      self%RLON2=FLOAT(IGDTMPL(15))*1.0E-6
      self%RLATI=FLOAT(IGDTMPL(13))*1.0E-6

      ISCAN=MOD(IGDTMPL(16)/128,2)
      JSCAN=MOD(IGDTMPL(16)/64,2)

      DY=FLOAT(IGDTMPL(19))*1.0E-3
      self%HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DPHI=HJ*DY/(self%RERTH*COS(self%RLATI/DPR))

      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
      self%nscan=mod(igdtmpl(16) / 32,2)
      self%nscan_field_pos = self%nscan

      self%iwrap = nint(360 / abs(self%dlon))
      if(self%im < self%iwrap) self%iwrap = 0

    end associate
  end subroutine init_grib2

  !> GDS wizard for mercator cylindrical.
  !>
  !> This routine decodes the grib 2 grid definition template (passed
  !> in integer form as decoded by the ncep g2 library) and returns
  !> one of the following:
  !> - (iopt=+1) earth coordinates of selected grid coordinates
  !> - (iopt=-1) grid coordinates of selected earth coordinates
  !>
  !> Works for mercator cylindrical projections.
  !>
  !> If the selected coordinates are more than one gridpoint beyond
  !> the the edges of the grid domain, then the relevant output
  !> elements are set to fill values.
  !>
  !> The actual number of valid points computed is returned too.
  !>
  !> Optionally, the vector rotations, map jacobians and the grid box
  !> areas may be returned. To compute the vector rotations, the
  !> optional arguments 'srot' and 'crot' must be present. To compute
  !> the map jacobians, the optional arguments 'xlon', 'xlat', 'ylon',
  !> 'ylat' must be present. to compute the grid box areas, the
  !> optional argument 'area' must be present.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | iredell | Initial
  !> 96-10-01 | iredell | protected against unresolvable points
  !> 97-10-20 | iredell | include map options
  !> 2015-01-21 | gayno | merger of gdswiz01() and gdswzd01(). Make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. Make part of a module. move vector rotation, map jacobian and grid box area computations to separate subroutines.
  !> 2015-07-13 | gayno | convert to grib 2. replace grib 1 kgds array with grib 2 grid definition template array. Rename.
  !> 2018-07-20 | wesley | add threads.
  !>
  !> @param[in] self grid descriptor.
  !> @param[in] iopt option flag
  !> - 1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param[in] npts maximum number of coordinates
  !> @param[in] fill fill value to set invalid output data (must be
  !> impossible value; suggested value: -9999.)
  !> @param[inout] xpts (npts) grid x point coordinates if iopt>0
  !> @param[inout] ypts (npts) grid y point coordinates if iopt>0
  !> @param[inout] rlon (npts) earth longitudes in degrees e if iopt<0
  !> (acceptable range: -360. to 360.)
  !> @param[inout] rlat (npts) earth latitudes in degrees n if iopt<0
  !> (acceptable range: -90. to 90.)
  !> @param[out] nret number of valid points computed
  !> @param[out] crot optional (npts) clockwise vector rotation cosines
  !> @param[out] srot optional (npts) clockwise vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] xlon optional (npts) dx/dlon in 1/degrees
  !> @param[out] xlat optional (npts) dx/dlat in 1/degrees
  !> @param[out] ylon optional (npts) dy/dlon in 1/degrees
  !> @param[out] ylat optional (npts) dy/dlat in 1/degrees
  !> @param[out] area optional (npts) area weights in m**2
  !> (proportional to the square of the map factor)
  !>
  !> @author Iredell @date 96-04-10
  SUBROUTINE GDSWZD_MERCATOR(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE
    !
    class(ip_mercator_grid), intent(in) :: self
    INTEGER,           INTENT(IN   ) :: IOPT, NPTS
    INTEGER,           INTENT(  OUT) :: NRET
    !
    REAL,              INTENT(IN   ) :: FILL
    REAL,              INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,              INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,    INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,    INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,    INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                          :: IM, JM, N
    !
    LOGICAL                          :: LROT, LMAP, LAREA
    !
    REAL                             :: HI
    REAL                             :: RLAT1, RLON1, RLON2, RLATI
    REAL                             :: XMAX, XMIN, YMAX, YMIN
    REAL                             :: YE
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1
    RLON2=self%rlon2
    RLATI=self%rlati

    HI=self%hi

    DLON=self%dlon
    DPHI=self%dphi
    rerth = self%rerth

    YE=1-LOG(TAN((RLAT1+90)/2/DPR))/DPHI
    XMIN=0
    XMAX=IM+1
    IF(IM.EQ.NINT(360/ABS(DLON))) XMAX=IM+2
    YMIN=0
    YMAX=JM+1
    NRET=0
    IF(PRESENT(CROT).AND.PRESENT(SROT))THEN
       LROT=.TRUE.
    ELSE
       LROT=.FALSE.
    ENDIF
    IF(PRESENT(XLON).AND.PRESENT(XLAT).AND.PRESENT(YLON).AND.PRESENT(YLAT))THEN
       LMAP=.TRUE.
    ELSE
       LMAP=.FALSE.
    ENDIF
    IF(PRESENT(AREA))THEN
       LAREA=.TRUE.
    ELSE
       LAREA=.FALSE.
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLON(N)=MOD(RLON1+DLON*(XPTS(N)-1)+3600,360.)
             RLAT(N)=2*ATAN(EXP(DPHI*(YPTS(N)-YE)))*DPR-90
             NRET=NRET+1
             IF(LROT)  CALL MERCATOR_VECT_ROT(CROT(N),SROT(N))
             IF(LMAP)  CALL MERCATOR_MAP_JACOB(RLAT(N),XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL MERCATOR_GRID_AREA(RLAT(N),AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LT.90) THEN
             XPTS(N)=1+HI*MOD(HI*(RLON(N)-RLON1)+3600,360.)/DLON
             YPTS(N)=YE+LOG(TAN((RLAT(N)+90)/2/DPR))/DPHI
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT)  CALL MERCATOR_VECT_ROT(CROT(N),SROT(N))
                IF(LMAP)  CALL MERCATOR_MAP_JACOB(RLAT(N),XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL MERCATOR_GRID_AREA(RLAT(N),AREA(N))
             ELSE
                XPTS(N)=FILL
                YPTS(N)=FILL
             ENDIF
          ELSE
             XPTS(N)=FILL
             YPTS(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_MERCATOR

  !> Vector rotation fields for mercator cylindrical grids.
  !>
  !> This subprogram computes the vector rotation sines and cosines
  !> for a mercator cylindrical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "mercator_vect_rot".
  !>
  !> @param[in] crot clockwise vector rotation cosines (real)
  !> @param[in] srot clockwise vector rotation sines (real)
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE MERCATOR_VECT_ROT(CROT,SROT)
    IMPLICIT NONE

    REAL,                INTENT(  OUT) :: CROT, SROT

    CROT=1.0
    SROT=0.0

  END SUBROUTINE MERCATOR_VECT_ROT

  !> Map jacobians for mercator cylindrical grids.
  !>
  !> This subprogram computes the map jacobians for a mercator
  !> cylindrical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "mercator_map_jacob"
  !>
  !> @param[in] rlat latitude in degrees (real)
  !> @param[out] xlon dx/dlon in 1/degrees (real)
  !> @param[out] xlat dx/dlat in 1/degrees (real)
  !> @param[out] ylon dy/dlon in 1/degrees (real)
  !> @param[out] ylat dy/dlat in 1/degrees (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE MERCATOR_MAP_JACOB(RLAT,XLON,XLAT,YLON,YLAT)
    IMPLICIT NONE

    REAL,                INTENT(IN   ) :: RLAT
    REAL,                INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    XLON=1./DLON
    XLAT=0.
    YLON=0.
    YLAT=1./DPHI/COS(RLAT/DPR)/DPR

  END SUBROUTINE MERCATOR_MAP_JACOB

  !> Grid box area for mercator cylindrical grids.
  !>
  !> This subprogram computes the grid box area for a mercator
  !> cylindrical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "mercator_grid_area"
  !>
  !> @param[in] rlat latitude of grid point in degrees (real)
  !> @param[out] area area weights in m**2 (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE MERCATOR_GRID_AREA(RLAT,AREA)
    IMPLICIT NONE

    REAL,              INTENT(IN   ) :: RLAT
    REAL,              INTENT(  OUT) :: AREA

    AREA=RERTH**2*COS(RLAT/DPR)**2*DPHI*DLON/DPR

  END SUBROUTINE MERCATOR_GRID_AREA

end module ip_mercator_grid_mod

