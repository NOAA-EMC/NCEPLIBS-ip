module ip_mercator_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_mercator_grid

  type, extends(ip_grid) :: ip_mercator_grid
     real :: rlat1, rlon1, rlon2, rlati, hi, dlon, dphi
   contains
     procedure :: init_grib1
     procedure :: init_grib2
     procedure :: gdswzd => gdswzd_mercator
  end type ip_mercator_grid

  REAL                         :: DLON
  REAL                         :: DPHI
  REAL                         :: RERTH

CONTAINS

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

  SUBROUTINE GDSWZD_MERCATOR(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GDSWZD_MERCATOR   GDS WIZARD FOR MERCATOR CYLINDRICAL
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS ROUTINE DECODES THE GRIB 2 GRID DEFINITION
    !           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
    !           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
    !             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
    !             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
    !           WORKS FOR MERCATOR CYLINDRICAL PROJECTIONS.
    !           IF THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
    !           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
    !           OUTPUT ELEMENTS ARE SET TO FILL VALUES. THE
    !           ACTUAL NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
    !           OPTIONALLY, THE VECTOR ROTATIONS, MAP JACOBIANS AND
    !           THE GRID BOX AREAS MAY BE RETURNED.  TO COMPUTE THE
    !           VECTOR ROTATIONS, THE OPTIONAL ARGUMENTS 'SROT' AND 'CROT'
    !           MUST BE PRESENT.  TO COMPUTE THE MAP JACOBIANS, THE
    !           OPTIONAL ARGUMENTS 'XLON', 'XLAT', 'YLON', 'YLAT' MUST BE
    !           PRESENT. TO COMPUTE THE GRID BOX AREAS, THE OPTIONAL 
    !           ARGUMENT 'AREA' MUST BE PRESENT.
    !
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   96-10-01  IREDELL  PROTECTED AGAINST UNRESOLVABLE POINTS
    !   97-10-20  IREDELL  INCLUDE MAP OPTIONS
    ! 2015-01-21  GAYNO    MERGER OF GDSWIZ01 AND GDSWZD01.  MAKE
    !                      CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
    !                      OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE.
    !                      MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
    !                      BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !                      RENAME AS "GDSWZD_MERCATOR".
    ! 2018-07-20  WESLEY   ADD THREADS.
    !
    ! USAGE:    CALL GDSWZD_MERCATOR(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS, &
    !                                FILL,XPTS,YPTS,RLON,RLAT,NRET, &
    !                                CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                MUST BE "10" FOR MERCATOR GRIDS.
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
    !                 (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
    !                 (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
    !                 (10): LATITUDE OF FIRST POINT, OCTETS 39-42
    !                 (11): LONGITUDE OF FIRST POINT, OCTETS 43-46
    !                 (12): RESOLUTION AND COMPONENT FLAGS, OCTET 47
    !                 (13): TANGENT LATITUDE, OCTETS 48-51
    !                 (14): LATITUDE OF LAST POINT, OCTETS 52-55
    !                 (15): LONGITUDE OF LAST POINT, OCTETS 56-59
    !                 (16): SCANNING MODE FLAGS, OCTET 60
    !                 (17): ORIENTATION OF GRID, OCTETS 61-64
    !                 (18): LONGITUDINAL GRID LENGTH, OCTETS 65-68
    !                 (19): LATITUDINAL GRID LENGTH, OCTETS 69-72
    !     IGDTLEN  - INTEGER NUMBER OF ELEMENTS (19) OF THE GRID DEFINITION
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
    REAL                             :: DY, HI
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
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE MERCATOR_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  MERCATOR_ERROR   ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.

    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL MERCATOR_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !
    !   INPUT ARGUMENT LIST:
    !     IOPT     - INTEGER OPTION FLAG
    !                (+1 TO COMPUTE EARTH COORDS OF SELECTED GRID COORDS)
    !                (-1 TO COMPUTE GRID COORDS OF SELECTED EARTH COORDS)
    !     NPTS     - INTEGER MAXIMUM NUMBER OF COORDINATES
    !     FILL     - REAL FILL VALUE TO SET INVALID OUTPUT DATA
    !                (MUST BE IMPOSSIBLE VALUE; SUGGESTED VALUE: -9999.)
    !   OUTPUT ARGUMENT LIST:
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IF IOPT<0
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IF IOPT<0
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
  END SUBROUTINE MERCATOR_ERROR
  !
  SUBROUTINE MERCATOR_VECT_ROT(CROT,SROT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  MERCATOR_VECT_ROT   VECTOR ROTATION FIELDS FOR
    !                                  MERCATOR CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
    !           COSINES FOR A MERCATOR CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "MERCATOR_VECT_ROT".
    !
    ! USAGE:    CALL MERCATOR_VECT_ROT(CROT,SROT)
    !
    !   INPUT ARGUMENT LIST:
    !     NONE
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
    !
    IMPLICIT NONE

    REAL,                INTENT(  OUT) :: CROT, SROT

    CROT=1.0
    SROT=0.0

  END SUBROUTINE MERCATOR_VECT_ROT
  !
  SUBROUTINE MERCATOR_MAP_JACOB(RLAT,XLON,XLAT,YLON,YLAT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  MERCATOR_MAP_JACOB  MAP JACOBIANS FOR
    !                                  MERCATOR CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
    !           A MERCATOR CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "MERCATOR_MAP_JACOB"
    !
    ! USAGE:  CALL MERCATOR_MAP_JACOB(RLAT,XLON,XLAT,YLON,YLAT)
    !
    !   INPUT ARGUMENT LIST:
    !     RLAT     - LATITUDE IN DEGREES (REAL)
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

    REAL,                INTENT(IN   ) :: RLAT
    REAL,                INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    XLON=1./DLON
    XLAT=0.
    YLON=0.
    YLAT=1./DPHI/COS(RLAT/DPR)/DPR

  END SUBROUTINE MERCATOR_MAP_JACOB
  !
  SUBROUTINE MERCATOR_GRID_AREA(RLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  MERCATOR_GRID_AREA  GRID BOX AREA FOR
    !                                  MERCATOR CYLINDRICAL GRIDS
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
    !           A MERCATOR CYLINDRICAL GRID.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "MERCATOR_GRID_AREA"
    !
    ! USAGE:  CALL MERCATOR_GRID_AREA(RLAT,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     RLAT     - LATITUDE OF GRID POINT IN DEGREES (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     AREA     - AREA WEIGHTS IN M**2 (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$

    IMPLICIT NONE

    REAL,              INTENT(IN   ) :: RLAT
    REAL,              INTENT(  OUT) :: AREA

    AREA=RERTH**2*COS(RLAT/DPR)**2*DPHI*DLON/DPR

  END SUBROUTINE MERCATOR_GRID_AREA

end module ip_mercator_grid_mod

