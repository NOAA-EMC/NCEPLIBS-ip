!> @file
!> @brief GDS wizard for polar stereographic azimuthal.
!>
!> @author Iredell @date 96-04-10

!> @brief GDS wizard for polar stereographic azimuthal.
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.20 Polar
!> stereographic
!> projection](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-20.shtml).
!>
!> @author Iredell @date 96-04-10
module ip_polar_stereo_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use ip_constants_mod, only: DPR, PI, PI2, PI4, RERTH_WGS84, E2_WGS84
  use earth_radius_mod
  implicit none

  private
  public :: ip_polar_stereo_grid

  type, extends(ip_grid) :: ip_polar_stereo_grid
     logical :: elliptical !< When true/false, computations are based on an elliptical/spherical earth.
     real :: rlat1 !< Latitude of the first grid point.
     real :: rlon1 !< Longitude of the first grid point.
     real :: orient !< Orientation longitude.
     real :: h !< Hemisphere flag. 0 - NH; 1 - SH.
     real :: dxs !< 'x'-direction grid length, adjusted by the scanning mode.
     real :: dys !< 'y'-direction grid length, adjusted by the scanning mode.
     real :: slatr !< Standard latitude of grid in radians.
     !> Rotation flag. When '0' the u/v vector components are relative
     !> to north/east. When '1' the u/v vector components are grid
     !> relative.
     integer :: irot
   contains
     !> Initializes a grid given a grib1_descriptor object.
     !> @return N/A @copydoc ip_polar_stereo_grid_mod::init_grib1
     procedure :: init_grib1
     !> Initializes a grid given a grib2_descriptor object.
     !> @return N/A @copydoc ip_polar_stereo_grid_mod::init_grib2
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates
     !> (iopt = -1).
     !> @return N/A @copydoc ip_polar_stereo_grid_mod::gdswzd_polar_stereo
     procedure :: gdswzd => gdswzd_polar_stereo
  end type ip_polar_stereo_grid

  INTEGER :: IROT !< Local copy of irot.
  REAL :: DE2 !< Square of DE.
  REAL :: DXS !< Local copy of dxs.
  REAL :: DYS !< Local copy of dys.
  REAL :: E2 !< Eccentricity squared.
  REAL :: RERTH !< Radius of the Earth.
  REAL :: H !< Local copy of h.
  REAL :: ORIENT !< Local copy of orient.
  REAL :: TINYREAL=TINY(1.0) !< Smallest positive real value (use for equality comparisons)

CONTAINS

  !> Initializes a polar stereographic grid given a grib1_descriptor
  !! object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Iredell @date 96-04-10  
  subroutine init_grib1(self, g1_desc)
    class(ip_polar_stereo_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    REAL, PARAMETER :: SLAT=60.0  ! standard latitude according grib1 standard

    real :: dx, dy, hi, hj
    integer :: iproj, iscan, jscan

    associate(kgds => g1_desc%gds)
      self%ELLIPTICAL=MOD(KGDS(6)/64,2).EQ.1

      if (.not. self%elliptical) then
         self%rerth = 6.3712E6
         self%eccen_squared = 0d0
      else
         self%rerth = RERTH_WGS84
         self%eccen_squared = E2_WGS84 !wgs84 datum
      end if

      self%IM=KGDS(2)
      self%JM=KGDS(3)

      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3

      self%IROT=MOD(KGDS(6)/8,2)

      self%SLATR=SLAT/DPR

      self%ORIENT=KGDS(7)*1.E-3

      DX=KGDS(8)
      DY=KGDS(9)

      IPROJ=MOD(KGDS(10)/128,2)
      ISCAN=MOD(KGDS(11)/128,2)
      JSCAN=MOD(KGDS(11)/64,2)

      self%H=(-1.)**IPROJ
      HI=(-1.)**ISCAN
      HJ=(-1.)**(1-JSCAN)

      IF(ABS(self%H+1.).LT.TINYREAL) self%ORIENT=self%ORIENT+180.

      self%DXS=DX*HI
      self%DYS=DY*HJ

      self%iwrap= 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
    end associate

  end subroutine init_grib1

  !> Initializes a polar stereographic grid given a grib2_descriptor
  !! object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g2_desc A grib2_descriptor
  !!
  !! @author Iredell @date 96-04-10  
  subroutine init_grib2(self, g2_desc)
    class(ip_polar_stereo_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real :: slat, dx, dy, hi, hj
    integer :: iproj, iscan, jscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%ELLIPTICAL = self%eccen_squared > 0.0

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%RLAT1=FLOAT(IGDTMPL(10))*1.E-6
      self%RLON1=FLOAT(IGDTMPL(11))*1.E-6

      self%IROT=MOD(IGDTMPL(12)/8,2)

      SLAT=FLOAT(ABS(IGDTMPL(13)))*1.E-6
      self%SLATR=SLAT/DPR

      self%ORIENT=FLOAT(IGDTMPL(14))*1.E-6

      DX=FLOAT(IGDTMPL(15))*1.E-3
      DY=FLOAT(IGDTMPL(16))*1.E-3

      IPROJ=MOD(IGDTMPL(17)/128,2)
      ISCAN=MOD(IGDTMPL(18)/128,2)
      JSCAN=MOD(IGDTMPL(18)/64,2)

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

  !> GDS wizard for polar stereographic azimuthal
  !>
  !> This subprogram decodes the grib 2 grid definition template
  !> (passed in integer form as decoded by the ncep g2 library) and
  !> returns one of the following:
  !> - (iopt=+1) earth coordinates of selected grid coordinates
  !> - (iopt=-1) grid coordinates of selected earth coordinates
  !>
  !> Works for polar stereographic azimuthal projections.
  !>
  !> If the selected coordinates are more than one gridpoint beyond
  !> the the edges of the grid domain, then the relevant output
  !> elements are set to fill values.
  !>
  !> The actual number of valid points computed is returned too.
  !>
  !> Optionally, the vector rotations, map jacobians, and grid box
  !> areas may be returned as well. Routine works for both spherical
  !> and elliptical earths with the exception of the map jacobians and
  !> grid box areas, which are only computed for spherical earths.
  !>
  !> To compute the vector rotations, the optional arguments 'srot'
  !> and 'crot' must be present. To compute the map jacobians, the
  !> optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be
  !> present.  to compute the grid box areas, the optional argument
  !> 'area' must be present.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | iredell | Initial
  !> 97-10-20 | iredell | include map options
  !> 09-05-13 | gayno | ensure area always positive
  !> 2015-01-21 | gayno | merger of gdswiz05 and gdswzd05. make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. make part of a module. move vector rotation, map jacobian and grid box area computations to separate subroutines. include option for elliptical earths.
  !> 2015-07-13 | gayno | convert to grib 2. replace grib 1 kgds array with grib 2 grid definition template array. rename routine.
  !> 2018-07-20 | wesley | add threading.
  !>
  !> @param[in] self grid
  !> @param[in] iopt option flag
  !> - 1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param[in] npts maximum number of coordinates
  !> @param[in] fill fill value to set invalid output data
  !> (must be impossible value; suggested value: -9999.)
  !> @param[inout] xpts (npts) grid x point coordinates if iopt>0
  !> @param[inout] ypts (npts) grid y point coordinates if iopt>0
  !> @param[inout] rlon (npts) earth longitudes in degrees e if iopt<0
  !> (acceptable range: -360. to 360.)
  !> @param[inout] rlat (npts) earth latitudes in degrees n if iopt<0
  !> (acceptable range: -90. to 90.)
  !> @param[out] nret number of valid points computed
  !> @param[out] crot optional (npts) clockwise vector rotation cosines
  !> @param[out] srot optional (npts) clockwise vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth;
  !> vgrid=srot*uearth+crot*vearth)
  !> @param[out] xlon optional (npts) dx/dlon in 1/degrees
  !> @param[out] xlat optional (npts) dx/dlat in 1/degrees
  !> @param[out] ylon optional (npts) dy/dlon in 1/degrees
  !> @param[out] ylat optional (npts) dy/dlat in 1/degrees
  !> @param[out] area optional (npts) area weights in m**2
  !> (proportional to the square of the map factor)
  !>
  !> @author Iredell @date 96-04-10
  SUBROUTINE GDSWZD_POLAR_STEREO(self,IOPT,NPTS, &
       FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE
    !

    class(ip_polar_stereo_grid), intent(in) :: self
    INTEGER,          INTENT(IN   ) :: IOPT, NPTS
    INTEGER,          INTENT(  OUT) :: NRET
    !
    REAL,             INTENT(IN   ) :: FILL
    REAL,             INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,             INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,   INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,   INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,   INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                         :: IM, JM
    INTEGER                         :: ITER, N
    !
    LOGICAL                         :: ELLIPTICAL, LROT, LMAP, LAREA
    !
    REAL                            :: ALAT, ALAT1, ALONG, DIFF
    REAL                            :: DI, DJ, DE
    REAL                            :: DR, E, E_OVER_2
    REAL                            :: MC, SLATR
    REAL                            :: RLAT1, RLON1, RHO, T, TC
    REAL                            :: XMAX, XMIN, YMAX, YMIN
    REAL                            :: XP, YP, DR2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL

    elliptical = self%elliptical
    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1

    IROT=self%irot
    SLATR=self%slatr
    ORIENT=self%orient

    H=self%h
    DXS=self%dxs
    DYS=self%dys

    rerth = self%rerth
    e2 = self%eccen_squared
    !
    ! FIND X/Y OF POLE
    IF (.NOT.ELLIPTICAL) THEN
       DE=(1.+SIN(SLATR))*RERTH
       DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
       XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/DXS
       YP=1+COS((RLON1-ORIENT)/DPR)*DR/DYS
       DE2=DE**2
    ELSE
       E=SQRT(E2)
       E_OVER_2=E*0.5
       ALAT=H*RLAT1/DPR
       ALONG = (RLON1-ORIENT)/DPR
       T=TAN(PI4-ALAT/2.)/((1.-E*SIN(ALAT))/  &
            (1.+E*SIN(ALAT)))**(E_OVER_2)
       TC=TAN(PI4-SLATR/2.)/((1.-E*SIN(SLATR))/  &
            (1.+E*SIN(SLATR)))**(E_OVER_2)
       MC=COS(SLATR)/SQRT(1.0-E2*(SIN(SLATR)**2))
       RHO=RERTH*MC*T/TC
       YP = 1.0 + RHO*COS(H*ALONG)/DYS
       XP = 1.0 - RHO*SIN(H*ALONG)/DXS
    ENDIF ! ELLIPTICAL
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
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       IF(.NOT.ELLIPTICAL)THEN
          !$OMP PARALLEL DO PRIVATE(N,DI,DJ,DR2) REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                DI=(XPTS(N)-XP)*DXS
                DJ=(YPTS(N)-YP)*DYS
                DR2=DI**2+DJ**2
                IF(DR2.LT.DE2*1.E-6) THEN
                   RLON(N)=0.
                   RLAT(N)=H*90.
                ELSE
                   RLON(N)=MOD(ORIENT+H*DPR*ATAN2(DI,-DJ)+3600,360.)
                   RLAT(N)=H*DPR*ASIN((DE2-DR2)/(DE2+DR2))
                ENDIF
                NRET=NRET+1
                IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
                IF(LMAP) CALL POLAR_STEREO_MAP_JACOB(RLON(N),RLAT(N),DR2, &
                     XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL POLAR_STEREO_GRID_AREA(RLAT(N),DR2,AREA(N))
             ELSE
                RLON(N)=FILL
                RLAT(N)=FILL
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ELSE ! ELLIPTICAL
          !$OMP PARALLEL DO PRIVATE(N,DI,DJ,RHO,T,ALONG,ALAT1,ALAT,DIFF) &
          !$OMP& REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND.  &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                DI=(XPTS(N)-XP)*DXS
                DJ=(YPTS(N)-YP)*DYS
                RHO=SQRT(DI*DI+DJ*DJ)
                T=(RHO*TC)/(RERTH*MC)
                IF(ABS(YPTS(N)-YP)<0.01)THEN
                   IF(DI>0.0) ALONG=ORIENT+H*90.0
                   IF(DI<=0.0) ALONG=ORIENT-H*90.0
                ELSE
                   ALONG=ORIENT+H*ATAN(DI/(-DJ))*DPR
                   IF(DJ>0) ALONG=ALONG+180.
                END IF
                ALAT1=PI2-2.0*ATAN(T)
                DO ITER=1,10
                   ALAT = PI2 - 2.0*ATAN(T*(((1.0-E*SIN(ALAT1))/  &
                        (1.0+E*SIN(ALAT1)))**(E_OVER_2)))
                   DIFF = ABS(ALAT-ALAT1)*DPR
                   IF (DIFF < 0.000001) EXIT
                   ALAT1=ALAT
                ENDDO
                RLAT(N)=H*ALAT*DPR
                RLON(N)=ALONG
                IF(RLON(N)<0.0) RLON(N)=RLON(N)+360.
                IF(RLON(N)>360.0) RLON(N)=RLON(N)-360.0
                NRET=NRET+1
                IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
             ELSE
                RLON(N)=FILL
                RLAT(N)=FILL
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF ! ELLIPTICAL
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       IF(.NOT.ELLIPTICAL)THEN
          !$OMP PARALLEL DO PRIVATE(N,DR,DR2) REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(ABS(RLON(N)).LT.(360.+TINYREAL).AND.ABS(RLAT(N)).LT.(90.+TINYREAL).AND. &
                  ABS(H*RLAT(N)+90).GT.TINYREAL) THEN
                DR=DE*TAN((90-H*RLAT(N))/2/DPR)
                DR2=DR**2
                XPTS(N)=XP+H*SIN((RLON(N)-ORIENT)/DPR)*DR/DXS
                YPTS(N)=YP-COS((RLON(N)-ORIENT)/DPR)*DR/DYS
                IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                     YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                   NRET=NRET+1
                   IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
                   IF(LMAP) CALL POLAR_STEREO_MAP_JACOB(RLON(N),RLAT(N),DR2, &
                        XLON(N),XLAT(N),YLON(N),YLAT(N))
                   IF(LAREA) CALL POLAR_STEREO_GRID_AREA(RLAT(N),DR2,AREA(N))
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
       ELSE  ! ELLIPTICAL CASE
          !$OMP PARALLEL DO PRIVATE(N,ALAT,ALONG,T,RHO) REDUCTION(+:NRET) SCHEDULE(STATIC)
          DO N=1,NPTS
             IF(ABS(RLON(N)).LT.(360+TINYREAL).AND.ABS(RLAT(N)).LT.(90+TINYREAL).AND.  &
                  ABS(H*RLAT(N)+90).GT.TINYREAL) THEN
                ALAT = H*RLAT(N)/DPR
                ALONG = (RLON(N)-ORIENT)/DPR
                T=TAN(PI4-ALAT*0.5)/((1.-E*SIN(ALAT))/  &
                     (1.+E*SIN(ALAT)))**(E_OVER_2)
                RHO=RERTH*MC*T/TC
                XPTS(N)= XP + RHO*SIN(H*ALONG) / DXS
                YPTS(N)= YP - RHO*COS(H*ALONG) / DYS
                IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND.  &
                     YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                   NRET=NRET+1
                   IF(LROT) CALL POLAR_STEREO_VECT_ROT(RLON(N),CROT(N),SROT(N))
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
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_POLAR_STEREO

  !> Vector rotation fields for polar stereographic grids.
  !>
  !>
  !> This subprogram computes the vector rotation sines and
  !> cosines for a polar stereographic azimuthal grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "polar_stereo_vect_rot"
  !>
  !> @param[in] rlon grid point longitude in degrees (real)
  !> @param[in] crot clockwise vector rotation cosines (real)
  !> @param[in] srot clockwise vector rotation sines (real)
  !> (ugrid=crot*uearth-srot*vearth;
  !> vgrid=srot*uearth+crot*vearth)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE POLAR_STEREO_VECT_ROT(RLON, CROT, SROT)
    IMPLICIT NONE

    REAL,             INTENT(IN   ) :: RLON
    REAL,             INTENT(  OUT) :: CROT, SROT

    IF(IROT.EQ.1) THEN
       CROT=H*COS((RLON-ORIENT)/DPR)
       SROT=SIN((RLON-ORIENT)/DPR)
    ELSE
       CROT=1.
       SROT=0.
    ENDIF

  END SUBROUTINE POLAR_STEREO_VECT_ROT

  !> Map jacobians for polar stereographic grids.
  !>
  !> This subprogram computes the map jacobians for
  !> a polar stereographic azimuthal grid (spherical
  !> earth).
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "polar_stereo_map_jacob"
  !> 2018-07-20 | wesley | pass in dr2 for threading.
  !>
  !> @param[in] rlon longitude in degrees (real)
  !> @param[in] rlat latitude in degrees (real)
  !> @param[in] dr2 squared distance from pole (real)
  !> @param[out] xlon dx/dlon in 1/degrees (real)
  !> @param[out] xlat dx/dlat in 1/degrees (real)
  !> @param[out] ylon dy/dlon in 1/degrees (real)
  !> @param[out] ylat dy/dlat in 1/degrees (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE POLAR_STEREO_MAP_JACOB(RLON,RLAT,DR2,XLON,XLAT,YLON,YLAT)
    IMPLICIT NONE

    REAL,             INTENT(IN   ) :: RLON, RLAT, DR2
    REAL,             INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    REAL                            :: CLAT, DE, DR

    IF(DR2.LT.DE2*1.E-6) THEN
       DE=SQRT(DE2)
       XLON=0.
       XLAT=-SIN((RLON-ORIENT)/DPR)/DPR*DE/DXS/2
       YLON=0.
       YLAT=H*COS((RLON-ORIENT)/DPR)/DPR*DE/DYS/2
    ELSE
       DR=SQRT(DR2)
       CLAT=COS(RLAT/DPR)
       XLON=H*COS((RLON-ORIENT)/DPR)/DPR*DR/DXS
       XLAT=-SIN((RLON-ORIENT)/DPR)/DPR*DR/DXS/CLAT
       YLON=SIN((RLON-ORIENT)/DPR)/DPR*DR/DYS
       YLAT=H*COS((RLON-ORIENT)/DPR)/DPR*DR/DYS/CLAT
    ENDIF

  END SUBROUTINE POLAR_STEREO_MAP_JACOB

  !> Grid box area for polar stereographic grids.
  !>
  !> This subprogram computes the grid box area for
  !> a polar stereographic azimuthal grid (spherical
  !> earth).
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "polar_stereo_grid_area".
  !> 2018-07-20 | wesley | pass in dr2 for threading.
  !>
  !> @param[in] rlat latitude of grid point in degrees (real)
  !> @param[in] dr2 squared distance from pole (real)
  !> @param[out] area area weights in m**2 (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE POLAR_STEREO_GRID_AREA(RLAT, DR2, AREA)
    IMPLICIT NONE

    REAL,             INTENT(IN   ) :: RLAT, DR2
    REAL,             INTENT(  OUT) :: AREA

    REAL                            :: CLAT

    IF(DR2.LT.DE2*1.E-6) THEN
       AREA=RERTH**2*ABS(DXS)*ABS(DYS)*4/DE2
    ELSE
       CLAT=COS(RLAT/DPR)
       AREA=RERTH**2*CLAT**2*ABS(DXS)*ABS(DYS)/DR2
    ENDIF

  END SUBROUTINE POLAR_STEREO_GRID_AREA

end module ip_polar_stereo_grid_mod
