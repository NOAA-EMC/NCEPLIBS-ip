module ip_lambert_conf_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  use constants_mod
  implicit none

  private
  public :: ip_lambert_conf_grid

  type, extends(ip_grid) :: ip_lambert_conf_grid
     real :: rlat1, rlon1, rlati1, rlati2, orient
     real :: dxs, dys, h
     integer :: irot
   contains
     procedure :: init_grib1
     procedure :: init_grib2
     procedure :: gdswzd => gdswzd_lambert_conf
  end type ip_lambert_conf_grid


  INTEGER                       :: IROT

  REAL                          :: AN, DXS, DYS, H
  REAL                          :: RERTH

contains

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

  SUBROUTINE GDSWZD_LAMBERT_CONF(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GDSWZD_LAMBERT_CONF   GDS WIZARD FOR LAMBERT CONFORMAL CONICAL
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION
    !           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
    !           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
    !             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
    !             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
    !           WORKS FOR LAMBERT CONFORMAL CONICAL PROJECTIONS.
    !           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
    !           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
    !           OUTPUT ELEMENTS ARE SET TO FILL VALUES.
    !           THE ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
    !           OPTIONALLY, THE VECTOR ROTATIONS, MAP JACOBIANS AND
    !           GRID BOX AREAS FOR THIS GRID MAY BE RETURNED AS WELL.
    !           TO COMPUTE THE VECTOR ROTATIONS, THE OPTIONAL ARGUMENTS 
    !           'SROT' AND 'CROT'  MUST BE PRESENT.  TO COMPUTE THE MAP
    !           JACOBIANS, THE OPTIONAL ARGUMENTS 'XLON', 'XLAT', 
    !           'YLON', 'YLAT' MUST BE PRESENT. TO COMPUTE THE GRID BOX 
    !           AREAS THE OPTIONAL ARGUMENT 'AREA' MUST BE PRESENT.
    !
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   96-10-01  IREDELL  PROTECTED AGAINST UNRESOLVABLE POINTS
    !   97-10-20  IREDELL  INCLUDE MAP OPTIONS
    ! 1999-04-27  GILBERT  CORRECTED MINOR ERROR CALCULATING VARIABLE AN
    !                      FOR THE SECANT PROJECTION CASE (RLATI1.NE.RLATI2).
    ! 2012-08-14  GAYNO    FIX PROBLEM WITH SH GRIDS.  ENSURE GRID BOX
    !                      AREA ALWAYS POSITIVE.
    ! 2015-01-21  GAYNO    MERGER OF GDSWIZ03 AND GDSWZD03.  MAKE
    !                      CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
    !                      OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE.
    !                      MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
    !                      BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !                      RENAME ROUTINE AS "GDSWZD_LAMBERT_CONF".
    ! 2018-07-20  WESLEY   ADD THREADS.
    !
    ! USAGE:    CALL GDSWZD_LAMBERT_CONF(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,
    !    &                               FILL,XPTS,YPTS,RLON,RLAT,NRET,
    !    &                               CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                MUST BE "30" FOR LAMBERT CONFORMAL GRIDS.
    !     IGDTMPL  - INTEGER (IGDTLEN) GRID DEFINITION TEMPLATE ARRAY.
    !                CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE FOR SECTION
    !                THREE:
    !                 (1):  SHAPE OF EARTH, OCTET 15
    !                 (2):  SCALE FACTOR OF SPHERICAL EARTH RADIUS,
    !                       OCTET 16
    !                 (3):  SCALED VALUE OF RADIUS OF SPHERICAL EARTH,
    !                       OCTETS 17-20
    !                 (4):  SCALE FACTOR OF MAJOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTET 21
    !                 (5):  SCALED VALUE OF MAJOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTETS 22-25
    !                 (6):  SCALE FACTOR OF MINOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTET 26
    !                 (7):  SCALED VALUE OF MINOR AXIS OF ELLIPTICAL EARTH,
    !                       OCTETS 27-30
    !                 (8):  NUMBER OF POINTS ALONG X-AXIS, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG Y-AXIS, OCTS 35-38
    !                 (10): LATITUDE OF FIRST POINT, OCTETS 39-42
    !                 (11): LONGITUDE OF FIRST POINT, OCTETS 43-46
    !                 (12): RESOLUTION OF COMPONENT FLAG, OCTET 47
    !                 (13): LATITUDE WHERE GRID LENGTHS SPECIFIED, 
    !                       OCTETS 48-51
    !                 (14): LONGITUDE OF MERIDIAN THAT IS PARALLEL TO
    !                       Y-AXIS, OCTETS 52-55
    !                 (15): X-DIRECTION GRID LENGTH, OCTETS 56-59
    !                 (16): Y-DIRECTION GRID LENGTH, OCTETS 60-63
    !                 (17): PROJECTION CENTER FLAG, OCTET 64
    !                 (18): SCANNING MODE, OCTET 65
    !                 (19): FIRST TANGENT LATITUDE FROM POLE, OCTETS 66-69
    !                 (20): SECOND TANGENT LATITUDE FROM POLE, OCTETS 70-73
    !                 (21): LATITUDE OF SOUTH POLE OF PROJECTION, 
    !                       OCTETS 74-77
    !                 (22): LONGITUDE OF SOUTH POLE OF PROJECTION, 
    !                       OCTETS 78-81
    !     IGDTLEN  - INTEGER NUMBER OF ELEMENTS (22) OF THE GRID DEFINITION
    !                TEMPLATE ARRAY.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IOPT     - INTEGER OPTION FLAG
    !                (+1 TO COMPUTE EARTH COORDS OF SELECTED GRID COORDS)
    !                (-1 TO COMPUTE GRID COORDS OF SELECTED EARTH COORDS)
    !     NPTS     - INTEGER MAXIMUM NUMBER OF COORDINATES
    !     FILL     - REAL FILL VALUE TO SET INVALID OUTPUT DATA
    !                (MUST BE IMPOSSIBLE VALUE; SUGGESTED VALUE: -9999.)
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT>0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT>0
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT<0
    !                (ACCEPTABLE RANGE: -360. TO 360.)
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT<0
    !                (ACCEPTABLE RANGE: -90. TO 90.)
    !
    !   OUTPUT ARGUMENT LIST:
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT<0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT<0
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT>0
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT>0
    !     NRET     - INTEGER NUMBER OF VALID POINTS COMPUTED
    !     CROT     - REAL, OPTIONAL (NPTS) CLOCKWISE VECTOR ROTATION COSINES
    !     SROT     - REAL, OPTIONAL (NPTS) CLOCKWISE VECTOR ROTATION SINES
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !     XLON     - REAL, OPTIONAL (NPTS) DX/DLON IN 1/DEGREES
    !     XLAT     - REAL, OPTIONAL (NPTS) DX/DLAT IN 1/DEGREES
    !     YLON     - REAL, OPTIONAL (NPTS) DY/DLON IN 1/DEGREES
    !     YLAT     - REAL, OPTIONAL (NPTS) DY/DLAT IN 1/DEGREES
    !     AREA     - REAL, OPTIONAL (NPTS) AREA WEIGHTS IN M**2
    !                (PROPORTIONAL TO THE SQUARE OF THE MAP FACTOR)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
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

    IF(RLATI1.EQ.RLATI2) THEN
       AN=SIN(RLATI1/DPR)
    ELSE
       AN=LOG(COS(RLATI1/DPR)/COS(RLATI2/DPR))/ &
            LOG(TAN((90-RLATI1)/2/DPR)/TAN((90-RLATI2)/2/DPR))
    ENDIF
    DE=RERTH*COS(RLATI1/DPR)*TAN((RLATI1+90)/2/DPR)**AN/AN
    IF(H*RLAT1.EQ.90) THEN
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
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90.AND. &
               H*RLAT(N).NE.-90) THEN
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
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE LAMBERT_CONF_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  LAMBERT_CONF_ERROR   ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.

    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL LAMBERT_CONF_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !
    !   INPUT ARGUMENT LIST:
    !     IOPT     - INTEGER OPTION FLAG
    !                (+1 TO COMPUTE EARTH COORDS OF SELECTED GRID COORDS)
    !                (-1 TO COMPUTE GRID COORDS OF SELECTED EARTH COORDS)
    !     NPTS     - INTEGER MAXIMUM NUMBER OF COORDINATES
    !     FILL     - REAL FILL VALUE TO SET INVALID OUTPUT DATA
    !                (MUST BE IMPOSSIBLE VALUE; SUGGESTED VALUE: -9999.)
    !   OUTPUT ARGUMENT LIST:
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT<0
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT<0
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT>0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT>0
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: IOPT, NPTS
    !
    REAL,    INTENT(IN   ) :: FILL
    REAL,    INTENT(  OUT) :: RLAT(NPTS),RLON(NPTS)
    REAL,    INTENT(  OUT) :: XPTS(NPTS),YPTS(NPTS)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(IOPT>=0) THEN
       RLON=FILL
       RLAT=FILL
    ENDIF
    IF(IOPT<=0) THEN
       XPTS=FILL
       YPTS=FILL
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE LAMBERT_CONF_ERROR
  !
  SUBROUTINE LAMBERT_CONF_VECT_ROT(DLON,CROT,SROT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  LAMBERT_CONF_VECT_ROT   VECTOR ROTATION FIELDS FOR
    !                                      LAMBERT CONFORMAL CONICAL
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
    !           COSINES FOR A LAMBERT CONFORMAL CONICAL GRID
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "LAMBERT_CONF_VECT_ROT"
    ! 2018-07-20  WESLEY   PASS IN DLON FOR THREADING.
    !
    ! USAGE:    CALL LAMBERT_CONF_VECT_ROT(DLON,CROT,SROT)
    !
    !   INPUT ARGUMENT LIST:
    !     DLON     - DISTANCE FROM ORIENTATION LONGITUDE (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     CROT     - CLOCKWISE VECTOR ROTATION COSINES (REAL)
    !     SROT     - CLOCKWISE VECTOR ROTATION SINES (REAL)
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
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
  !
  SUBROUTINE LAMBERT_CONF_MAP_JACOB(RLAT,FILL,DLON,DR,XLON,XLAT,YLON,YLAT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  LAMBERT_CONF_MAP_JACOB  MAP JACOBIANS FOR
    !                                      LAMBERT CONFORMAL CONICAL
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
    !           A LAMBERT CONFORMAL CONICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "LAMBERT_CONF_MAP_JACOB"
    ! 2018-07-20  WESLEY   PASS DLON AND DR FOR THREADING.
    !
    ! USAGE:  CALL LAMBERT_CONF_MAP_JACOB(RLAT,FILL,DLON,DR,XLON,XLAT,YLON,YLAT)
    !
    !   INPUT ARGUMENT LIST:
    !     RLAT     - GRID POINT LATITUDE IN DEGREES (REAL)
    !     FILL     - FILL VALUE FOR UNDEFINED POINTS (REAL)
    !     DLON     - DISTANCE FROM ORIENTATION LONGITUDE (REAL)
    !     DR       - DISTANCE FROM POLE POINT (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     XLON     - DX/DLON IN 1/DEGREES (REAL)
    !     XLAT     - DX/DLAT IN 1/DEGREES (REAL)
    !     YLON     - DY/DLON IN 1/DEGREES (REAL)
    !     YLAT     - DY/DLAT IN 1/DEGREES (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$

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
  !
  SUBROUTINE LAMBERT_CONF_GRID_AREA(RLAT,FILL,DR,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  LAMBERT_CONF_GRID_AREA  GRID BOX AREA FOR
    !                                      LAMBERT CONFORMAL CONICAL
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
    !           A LAMBERT CONFORMAL CONICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "LAMBERT_CONF_GRID_AREA"
    ! 2018-07-20  WESLEY   PASS IN DR FOR THREADING.
    !
    ! USAGE:  CALL LAMBERT_CONF_GRID_AREA(RLAT,FILL,DR,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     RLAT     - LATITUDE OF GRID POINT IN DEGREES (REAL)
    !     FILL     - FILL VALUE FOR UNDEFINED POINTS (REAL)
    !     DR       - DISTANCE FROM POLE POINT (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     AREA     - AREA WEIGHTS IN M**2 (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$

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

