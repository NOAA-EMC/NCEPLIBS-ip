!> @file
!> @brief GDS wizard for lambert conformal conical.
!>
!> @author Iredell @date 96-04-10

!> @brief Lambert conformal grib decoder and grid coordinate
!! transformations.
!!
!! Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.30
!! Lambert
!! conformal](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-30.shtml).
!!
!> @author Iredell @date 96-04-10
module ip_lambert_conf_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  use ip_constants_mod
  implicit none

  private
  public :: ip_lambert_conf_grid

  type, extends(ip_grid) :: ip_lambert_conf_grid
     real :: rlat1 !< La1― latitude of first grid point. GRIB2, Section 3.30, octet 39-42.
     real :: rlon1 !< Lo1― longitude of first grid point. GRIB2, Section 3.30, octet 43-46.
     real :: rlati1 !< First latitude from the pole at which the secant cone cuts the sphere. GRIB2, Section 3, octets 66-69.
     real :: rlati2 !< Second latitude from the pole at which the scant cone cuts the sphere. GRIB2, Section 3, octets 70-73.
     real :: orient !<  Longitude of meridian parallel to y-axis along which latitude increases at the latitude increases. GRIB2, Section 3, octets 52-55.
     real :: dxs !< x-direction grid length adjusted for scan mode. GRIB2, Section 3, octets 56-59.
     real :: dys !< y-direction grid length adjusted for scan model. GRIB2, Section 3, octets 60-63.
     real :: h !< Hemisphere flag. 1-NH, minus 1-SH.
     integer :: irot !< vector rotation flag. When "1", vectors are grid relative. When "0", vectors are earth relative. GRIB2, Section 3, octet 55.
   contains
     !> Initializes a gaussian grid given a grib1_descriptor object.
     !> @return N/A @copydoc ip_lambert_conf_grid_mod::init_grib1
     procedure :: init_grib1
     !> Initializes a gaussian grid given a grib2_descriptor object.
     !> @return N/A @copydoc ip_lambert_conf_grid_mod::init_grib2
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
     !> for Gaussian grids.
     !> @return N/A @copydoc ip_lambert_conf_grid_mod::gdswzd_lambert_conf
     procedure :: gdswzd => gdswzd_lambert_conf
  end type ip_lambert_conf_grid

  INTEGER :: IROT !< vector rotation flag. When "1", vectors are grid relative. When "0", vectors are earth relative. GRIB2, Section 3, octet 55.
  REAL :: AN !< Cone factor
  REAL :: DXS !< x-direction grid length adjusted for scan mode. GRIB2, Section 3, octets 56-59.
  REAL :: DYS !< y-direction grid length adjusted for scan model. GRIB2, Section 3, octets 60-63.
  REAL :: H !<  Hemisphere flag. 1-NH, minus 1-SH.
  REAL :: RERTH !< Radius of the earth. GRIB2, Section 3, octets 15-30.
  REAL :: TINYREAL=TINY(1.0) !< Smallest positive real value (use for equality comparisons)

contains

  !> Initializes a Lambert Conformal grid given a grib1_descriptor object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Iredell @date 96-04-10  
  subroutine init_grib1(self, g1_desc)
    class(ip_lambert_conf_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    real :: dx, dy, hi, hj
    integer :: iproj, iscan, jscan

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.0

      self%IM=KGDS(2)
      self%JM=KGDS(3)

      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3

      self%IROT=MOD(KGDS(6)/8,2)
      self%ORIENT=KGDS(7)*1.E-3

      DX=KGDS(8)
      DY=KGDS(9)

      IPROJ=MOD(KGDS(10)/128,2)
      ISCAN=MOD(KGDS(11)/128,2)
      JSCAN=MOD(KGDS(11)/64,2)

      self%RLATI1=KGDS(12)*1.E-3
      self%RLATI2=KGDS(13)*1.E-3
      self%H=(-1.)**IPROJ

      HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)
      self%DXS=DX*HI
      self%DYS=DY*HJ

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
    end associate

  end subroutine init_grib1

  !> Initializes a Lambert Conformal grid given a grib2_descriptor object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g2_desc A grib2_descriptor
  !!
  !! @author Iredell @date 96-04-10  
  subroutine init_grib2(self, g2_desc)
    class(ip_lambert_conf_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real :: dx, dy, hi, hj
    integer :: iproj, iscan, jscan


    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%RLAT1=FLOAT(IGDTMPL(10))*1.0E-6
      self%RLON1=FLOAT(IGDTMPL(11))*1.0E-6

      self%IROT=MOD(IGDTMPL(12)/8,2)
      self%ORIENT=FLOAT(IGDTMPL(14))*1.0E-6

      DX=FLOAT(IGDTMPL(15))*1.0E-3
      DY=FLOAT(IGDTMPL(16))*1.0E-3

      IPROJ=MOD(IGDTMPL(17)/128,2)
      ISCAN=MOD(IGDTMPL(18)/128,2)
      JSCAN=MOD(IGDTMPL(18)/64,2)

      self%RLATI1=FLOAT(IGDTMPL(19))*1.0E-6
      self%RLATI2=FLOAT(IGDTMPL(20))*1.0E-6

      self%H=(-1.)**IPROJ
      HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)
      self%DXS=DX*HI
      self%DYS=DY*HJ

      self%nscan = mod(igdtmpl(18) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
    end associate
  end subroutine init_grib2

  !> GDS wizard for lambert conformal conical.
  !>
  !> This subprogram decodes the grib 2 grid definition template
  !> (passed in integer form as decoded by the ncep g2 library) and
  !> returns one of the following:
  !> - (iopt=+1) earth coordinates of selected grid coordinates
  !> - (iopt=-1) grid coordinates of selected earth coordinates
  !>
  !> Works for lambert conformal conical projections.
  !>
  !> If the selected coordinates are more than one gridpoint beyond
  !> the the edges of the grid domain, then the relevant output
  !> elements are set to fill values.
  !>
  !> The actual number of valid points computed is returned too.
  !>
  !> Optionally, the vector rotations, map jacobians and grid box
  !> areas for this grid may be returned as well.
  !>
  !> To compute the vector rotations, the optional arguments 'srot'
  !> and 'crot' must be present. To compute the map jacobians, the
  !> optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be
  !> present. To compute the grid box areas the optional argument
  !> 'area' must be present.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | iredell | Initial.
  !> 96-10-01 | iredell | protected against unresolvable points
  !> 97-10-20 | iredell | include map options
  !> 1999-04-27 | gilbert | corrected minor error calculating variable an for the secant projection case (rlati1.ne.rlati2).
  !> 2012-08-14 | gayno | fix problem with sh grids. Ensure grid box area always positive.
  !> 2015-01-21 | gayno | merger of gdswiz03() and gdswzd03(). Make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. Make part of a module. Move vector rotation, map jacobian and grid box area computations to separate subroutines.
  !> 2015-07-13 | gayno | Convert to grib 2. Replace grib 1 kgds array with grib 2 grid definition template array. Rename routine.
  !> 2018-07-20 | wesley | add threads.
  !>
  !> @param[in] self ip_lambert_conf_grid object.
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
  SUBROUTINE GDSWZD_LAMBERT_CONF(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE
    !
    class(ip_lambert_conf_grid), intent(in) :: self
    INTEGER,        INTENT(IN   ) :: IOPT, NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,           INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                       :: IM, JM, N
    !
    LOGICAL                       :: LROT, LMAP, LAREA
    !
    REAL                          :: ANTR, DI, DJ
    REAL                          :: DLON1
    REAL                          :: DE, DE2, DR2
    REAL                          :: ORIENT, RLAT1, RLON1
    REAL                          :: RLATI1, RLATI2
    REAL                          :: XMAX, XMIN, YMAX, YMIN, XP, YP
    REAL                          :: DLON, DR
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL

    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1

    IROT=self%irot
    ORIENT=self%orient

    RLATI1=self%rlati1
    RLATI2=self%rlati2

    H=self%h
    DXS=self%dxs
    DYS=self%dys

    rerth = self%rerth

    IF(ABS(RLATI1-RLATI2).LT.TINYREAL) THEN
       AN=SIN(RLATI1/DPR)
    ELSE
       AN=LOG(COS(RLATI1/DPR)/COS(RLATI2/DPR))/ &
            LOG(TAN((90-RLATI1)/2/DPR)/TAN((90-RLATI2)/2/DPR))
    ENDIF
    DE=RERTH*COS(RLATI1/DPR)*TAN((RLATI1+90)/2/DPR)**AN/AN
    IF(ABS(H*RLAT1-90).LT.TINYREAL) THEN
       XP=1
       YP=1
    ELSE
       DR=DE/TAN((RLAT1+90)/2/DPR)**AN
       DLON1=MOD(RLON1-ORIENT+180+3600,360.)-180
       XP=1-SIN(AN*DLON1/DPR)*DR/DXS
       YP=1+COS(AN*DLON1/DPR)*DR/DYS
    ENDIF
    ANTR=1/(2*AN)
    DE2=DE**2
    XMIN=0
    XMAX=IM+1
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
    ! TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       !$OMP PARALLEL DO PRIVATE(N,DI,DJ,DR2,DR,DLON) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             DI=H*(XPTS(N)-XP)*DXS
             DJ=H*(YPTS(N)-YP)*DYS
             DR2=DI**2+DJ**2
             DR=SQRT(DR2)
             IF(DR2.LT.DE2*1.E-6) THEN
                RLON(N)=0.
                RLAT(N)=H*90.
             ELSE
                RLON(N)=MOD(ORIENT+1./AN*DPR*ATAN2(DI,-DJ)+3600,360.)
                RLAT(N)=(2*DPR*ATAN((DE2/DR2)**ANTR)-90)
             ENDIF
             NRET=NRET+1
             DLON=MOD(RLON(N)-ORIENT+180+3600,360.)-180
             IF(LROT)  CALL LAMBERT_CONF_VECT_ROT(DLON,CROT(N),SROT(N))
             IF(LMAP)  CALL LAMBERT_CONF_MAP_JACOB(RLAT(N),FILL, DLON, DR, &
                  XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL LAMBERT_CONF_GRID_AREA(RLAT(N),FILL,DR,AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N,DR,DLON) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(ABS(RLON(N)).LT.(360.+TINYREAL).AND.ABS(RLAT(N)).LT.(90.+TINYREAL).AND. &
               ABS(H*RLAT(N)+90).GT.TINYREAL) THEN
             DR=H*DE*TAN((90-RLAT(N))/2/DPR)**AN
             DLON=MOD(RLON(N)-ORIENT+180+3600,360.)-180
             XPTS(N)=XP+H*SIN(AN*DLON/DPR)*DR/DXS
             YPTS(N)=YP-H*COS(AN*DLON/DPR)*DR/DYS
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT)  CALL LAMBERT_CONF_VECT_ROT(DLON,CROT(N),SROT(N))
                IF(LMAP)  CALL LAMBERT_CONF_MAP_JACOB(RLAT(N),FILL,DLON,DR, &
                     XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL LAMBERT_CONF_GRID_AREA(RLAT(N),FILL,DR,AREA(N))
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
  END SUBROUTINE GDSWZD_LAMBERT_CONF

  !> Vector rotation fields for lambert conformal conical.
  !>
  !> This subprogram computes the vector rotation sines and
  !> cosines for a lambert conformal conical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "lambert_conf_vect_rot"
  !> 2018-07-20 | wesley | pass in dlon for threading.
  !>
  !> @param[in] dlon from orientation longitude (real)
  !> @param[out] crot vector rotation cosines (real)
  !> @param[out] srot vector rotation sines (real)
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE LAMBERT_CONF_VECT_ROT(DLON,CROT,SROT)
    IMPLICIT NONE
    REAL,           INTENT(   IN) :: DLON
    REAL,           INTENT(  OUT) :: CROT, SROT

    IF(IROT.EQ.1) THEN
       CROT=COS(AN*DLON/DPR)
       SROT=SIN(AN*DLON/DPR)
    ELSE
       CROT=1.
       SROT=0.
    ENDIF

  END SUBROUTINE LAMBERT_CONF_VECT_ROT

  !> Map jacobians for lambert conformal conical.
  !>
  !> This subprogram computes the map jacobians for a lambert
  !> conformal conical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | Gayno | initial version
  !> 2015-09-17 | Gayno | rename as "lambert_conf_map_jacob"
  !> 2018-07-20 | Wesley | pass dlon and dr for threading.
  !>
  !> @param[in] rlat grid point latitude in degrees (real)
  !> @param[in] fill fill value for undefined points (real)
  !> @param[in] dlon distance from orientation longitude (real)
  !> @param[in] dr distance from pole point (real)
  !> @param[out] xlon dx/dlon in 1/degrees (real)
  !> @param[out] xlat dx/dlat in 1/degrees (real)
  !> @param[out] ylon dy/dlon in 1/degrees (real)
  !> @param[out] ylat dy/dlat in 1/degrees (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE LAMBERT_CONF_MAP_JACOB(RLAT,FILL,DLON,DR,XLON,XLAT,YLON,YLAT)
    IMPLICIT NONE

    REAL,           INTENT(IN   ) :: RLAT, FILL, DLON, DR
    REAL,           INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    REAL                          :: CLAT

    CLAT=COS(RLAT/DPR)
    IF(CLAT.LE.0.OR.DR.LE.0) THEN
       XLON=FILL
       XLAT=FILL
       YLON=FILL
       YLAT=FILL
    ELSE
       XLON=H*COS(AN*DLON/DPR)*AN/DPR*DR/DXS
       XLAT=-H*SIN(AN*DLON/DPR)*AN/DPR*DR/DXS/CLAT
       YLON=H*SIN(AN*DLON/DPR)*AN/DPR*DR/DYS
       YLAT=H*COS(AN*DLON/DPR)*AN/DPR*DR/DYS/CLAT
    ENDIF

  END SUBROUTINE LAMBERT_CONF_MAP_JACOB

  !> Grid box area for lambert conformal conical.
  !>
  !> This subprogram computes the grid box area for a lambert
  !> conformal conical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | Gayno | initial version
  !> 2015-09-17 | Gayno | rename as "lambert_conf_grid_area"
  !> 2018-07-20 | Wesley | pass in dr for threading.
  !>
  !> @param[in] rlat latitude of grid point in degrees (real)
  !> @param[in] fill fill value for undefined points (real)
  !> @param[in] dr distance from pole point (real)
  !> @param[out] area area weights in m**2 (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE LAMBERT_CONF_GRID_AREA(RLAT,FILL,DR,AREA)
    IMPLICIT NONE

    REAL,           INTENT(IN   ) :: RLAT
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(IN   ) :: DR
    REAL,           INTENT(  OUT) :: AREA

    REAL                          :: CLAT

    CLAT=COS(RLAT/DPR)
    IF(CLAT.LE.0.OR.DR.LE.0) THEN
       AREA=FILL
    ELSE
       AREA=RERTH**2*CLAT**2*ABS(DXS)*ABS(DYS)/(AN*DR)**2
    ENDIF

  END SUBROUTINE LAMBERT_CONF_GRID_AREA

end module ip_lambert_conf_grid_mod

