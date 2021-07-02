module ip_rot_equid_cylind_grid_mod
  use iso_fortran_env, only: real64
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_rot_equid_cylind_grid

  integer, parameter :: kd = real64

  type, extends(ip_grid) :: ip_rot_equid_cylind_grid
     real(kd) :: clat0, dlats, dlons, rlon0, slat0, wbd, sbd
     integer :: irot
   contains
     procedure :: init_grib1
     procedure :: init_grib2
     procedure :: gdswzd => gdswzd_rot_equid_cylind
  end type ip_rot_equid_cylind_grid

  INTEGER                                :: IROT

  REAL(KIND=KD)                          :: RERTH
  REAL(KIND=KD)                          :: CLAT0, DLATS, DLONS
  REAL(KIND=KD)                          :: RLON0, SLAT0


CONTAINS

  subroutine init_grib1(self, g1_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    real(kd) :: hs, hs2, slat1, slat2, slatr, clon1, clon2, clat1, clat2, clatr, clonr, rlonr, rlatr
    integer :: iscale

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6_KD
      self%eccen_squared = 0d0

      RLAT1=KGDS(4)*1.E-3_KD
      RLON1=KGDS(5)*1.E-3_KD
      RLAT0=KGDS(7)*1.E-3_KD
      self%RLON0=KGDS(8)*1.E-3_KD
      RLAT2=KGDS(12)*1.E-3_KD
      RLON2=KGDS(13)*1.E-3_KD

      self%IROT=MOD(KGDS(6)/8,2)
      self%IM=KGDS(2)
      self%JM=KGDS(3)

      SLAT1=SIN(RLAT1/DPR)
      CLAT1=COS(RLAT1/DPR)
      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      HS=SIGN(1._KD,MOD(RLON1-self%RLON0+180+3600,360._KD)-180)
      CLON1=COS((RLON1-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT1-self%SLAT0*CLAT1*CLON1
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT1*CLON1+self%SLAT0*SLAT1)/CLATR
      RLATR=DPR*ASIN(SLATR)
      RLONR=HS*DPR*ACOS(CLONR)

      self%WBD=RLONR
      self%SBD=RLATR
      SLAT2=SIN(RLAT2/DPR)
      CLAT2=COS(RLAT2/DPR)
      HS2=SIGN(1._KD,MOD(RLON2-self%RLON0+180+3600,360._KD)-180)
      CLON2=COS((RLON2-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT2-self%SLAT0*CLAT2*CLON2
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT2*CLON2+self%SLAT0*SLAT2)/CLATR
      NBD=DPR*ASIN(SLATR)
      EBD=HS2*DPR*ACOS(CLONR)
      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
    end associate

  end subroutine init_grib1

  subroutine init_grib2(self, g2_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    integer :: iscale
    integer :: i_offset_odd, i_offset_even, j_offset

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

      CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,self%rerth,self%eccen_squared)

      I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
      J_OFFSET=MOD(IGDTMPL(19)/2,2)

      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6

      RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
      RLAT0=RLAT0+90.0_KD

      self%RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)

      RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
      RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)

      self%IROT=MOD(IGDTMPL(14)/8,2)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      self%WBD=RLON1
      IF (self%WBD > 180.0) self%WBD = self%WBD - 360.0
      self%SBD=RLAT1

      NBD=RLAT2
      EBD=RLON2

      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)

      IF(I_OFFSET_ODD==1) self%WBD=self%WBD+(0.5_KD*self%DLONS)
      IF(J_OFFSET==1) self%SBD=self%SBD+(0.5_KD*self%DLATS)

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
      self%nscan = mod(igdtmpl(19) / 32, 2)
      self%nscan_field_pos = self%nscan
    end associate
  end subroutine init_grib2



  SUBROUTINE GDSWZD_ROT_EQUID_CYLIND(self,IOPT,NPTS, &
       FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GDSWZD_ROT_EQUID_CYLIND  GDS WIZARD FOR ROTATED 
    !                                       EQUIDISTANT CYLINDRICAL
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2007-NOV-15
    !
    ! ABSTRACT: THIS SUBPROGRAM DECODES THE GRIB 2 GRID DEFINITION
    !           TEMPLATE (PASSED IN INTEGER FORM AS DECODED BY THE
    !           NCEP G2 LIBRARY) AND RETURNS ONE OF THE FOLLOWING:
    !             (IOPT=+1) EARTH COORDINATES OF SELECTED GRID COORDINATES
    !             (IOPT=-1) GRID COORDINATES OF SELECTED EARTH COORDINATES
    !           WORKS FOR NON-"E" STAGGERED ROTATED EQUIDISTANT CYLINDRICAL 
    !           PROJECTIONS. THE SCAN MODE (SECTION 3, OCTET 72, BITS 5-6)
    !           DETERMINE WHETHER THIS IS AN "H" OR "V" GRID.  IF
    !           THE SELECTED COORDINATES ARE MORE THAN ONE GRIDPOINT
    !           BEYOND THE THE EDGES OF THE GRID DOMAIN, THEN THE RELEVANT
    !           OUTPUT ELEMENTS ARE SET TO FILL VALUES.   THE ACTUAL 
    !           NUMBER OF VALID POINTS COMPUTED IS RETURNED TOO.
    !           OPTIONALLY, THE VECTOR ROTATIONS, THE MAP JACOBIANS AND
    !           THE GRID BOX AREAS MAY BE RETURNED AS WELL.  TO COMPUTE
    !           THE VECTOR ROTATIONS, THE OPTIONAL ARGUMENTS 'SROT' AND 'CROT'
    !           MUST BE PRESENT.  TO COMPUTE THE MAP JACOBIANS, THE
    !           OPTIONAL ARGUMENTS 'XLON', 'XLAT', 'YLON', 'YLAT' MUST 
    !           BE PRESENT. TO COMPUTE THE GRID BOX AREAS, THE OPTIONAL 
    !           ARGUMENT 'AREA' MUST BE PRESENT.
    !
    ! PROGRAM HISTORY LOG:
    ! 2010-JAN-15  GAYNO     BASED ON ROUTINES GDSWZDCB AND GDSWZDCA
    ! 2015-JAN-21  GAYNO     MERGER OF GDSWIZCD AND GDSWZDCD.  MAKE
    !                        CROT,SORT,XLON,XLAT,YLON,YLAT AND AREA
    !                        OPTIONAL ARGUMENTS.  MAKE PART OF A MODULE.
    !                        MOVE VECTOR ROTATION, MAP JACOBIAN AND GRID
    !                        BOX AREA COMPUTATIONS TO SEPARATE SUBROUTINES.
    ! 2015-JUL-13  GAYNO     CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAY
    !                        WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAY.
    !                        RENAME AS "GDSWZD_ROT_EQUID_CYLIND."
    ! 2018-07-20   WESLEY    ADD THREADS.
    !
    ! USAGE:   CALL GDSWZD_ROT_EQUID_CYLIND(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,
    !     &                                 FILL,XPTS,YPTS,RLON,RLAT,NRET,
    !     &                                 CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    !
    !   INPUT ARGUMENT LIST:
    !     IGDTNUM  - INTEGER GRID DEFINITION TEMPLATE NUMBER.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
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
    !                 (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN,
    !                       OCTETS 39-42
    !                 (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
    !                 (12): LATITUDE OF FIRST GRID POINT IN X/Y SPACE
    !                       (BEFORE ROTATION), OCTETS 47-50
    !                 (13): LONGITUDE OF FIRST GRID POINT IN X/Y
    !                       SPACE (BEFORE ROTATION), OCTETS 51-54
    !                 (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
    !                 (15): LATITUDE OF LAST GRID POINT IN X/Y SPACE
    !                       (BEFORE ROTATION), OCTETS 56-59
    !                 (16): LONGITUDE OF LAST GRID POINT IN X/Y SPACE
    !                       (BEFORE ROTATION), OCTETS 60-63
    !                 (17): I-DIRECTION INCREMENT, OCTETS 64-67
    !                 (18): J-DIRECTION INCREMENT, OCTETS 68-71
    !                 (19): SCANNING MODE, OCTET 72
    !                 (20): LATITUDE OF SOUTHERN POLE OF PROJECTION,
    !                       OCTETS 73-76
    !                 (21): LONGITUDE OF SOUTHERN POLE OF PROJECTION,
    !                       OCTETS 77-80
    !                 (22): ANGLE OF ROTATION OF PROJECTION, OCTS 81-84
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
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE

    class(ip_rot_equid_cylind_grid), intent(in) :: self
    INTEGER,                 INTENT(IN   ) :: IOPT, NPTS
    INTEGER,                 INTENT(  OUT) :: NRET
    !
    REAL,                    INTENT(IN   ) :: FILL
    REAL,                    INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,                    INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                                :: IM,JM,N
    !
    LOGICAL                                :: LROT, LMAP, LAREA
    !
    REAL                                   :: DUM1,DUM2
    REAL(KIND=KD)                          :: HS
    REAL(KIND=KD)                          :: CLONR,CLATR,SLATR
    REAL(KIND=KD)                          :: CLAT,SLAT,CLON
    REAL(KIND=KD)                          :: RLATR,RLONR
    REAL(KIND=KD)                          :: WBD,SBD
    REAL                                   :: XMIN,XMAX,YMIN,YMAX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IS THE EARTH RADIUS DEFINED?
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  IS THIS AN "E"-STAGGER GRID?  ROUTINE CAN'T PROCESS THOSE.
    ! I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
    ! I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
    ! J_OFFSET=MOD(IGDTMPL(19)/2,2)
    ! IF(I_OFFSET_ODD/=I_OFFSET_EVEN) THEN
    !    CALL ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !    RETURN
    ! ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    RLON0=self%rlon0
    IROT=self%irot

    IM=self%im
    JM=self%jm

    SLAT0=self%slat0
    CLAT0=self%clat0

    WBD=self%wbd
    SBD=self%sbd

    DLATS=self%dlats
    DLONS=self%dlons

    XMIN=0
    XMAX=IM+1
    YMIN=0
    YMAX=JM+1
    NRET=0

    rerth = self%rerth
    IF(RERTH<0.)THEN
       CALL ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
       RETURN
    ENDIF

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
       !$OMP PARALLEL DO PRIVATE(N,RLONR,RLATR,HS,CLONR,SLATR,CLATR,SLAT,CLAT,CLON) &
       !$OMP& REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLONR=WBD+(XPTS(N)-1._KD)*DLONS
             RLATR=SBD+(YPTS(N)-1._KD)*DLATS
             IF(RLONR <= 0._KD) THEN
                HS=-1.0_KD
             ELSE
                HS=1.0_KD
             ENDIF
             CLONR=COS(RLONR/DPR)
             SLATR=SIN(RLATR/DPR)
             CLATR=COS(RLATR/DPR)
             SLAT=CLAT0*SLATR+SLAT0*CLATR*CLONR
             IF(SLAT.LE.-1) THEN
                CLAT=0.
                CLON=COS(RLON0/DPR)
                RLON(N)=0.
                RLAT(N)=-90.
             ELSEIF(SLAT.GE.1) THEN
                CLAT=0.
                CLON=COS(RLON0/DPR)
                RLON(N)=0.
                RLAT(N)=90.
             ELSE
                CLAT=SQRT(1-SLAT**2)
                CLON=(CLAT0*CLATR*CLONR-SLAT0*SLATR)/CLAT
                CLON=MIN(MAX(CLON,-1._KD),1._KD)
                RLON(N)=MOD(RLON0+HS*DPR*ACOS(CLON)+3600,360._KD)
                RLAT(N)=DPR*ASIN(SLAT)
             ENDIF
             NRET=NRET+1
             IF(LROT) CALL ROT_EQUID_CYLIND_VECT_ROT(RLON(N), CLATR, SLATR, &
                  CLAT, SLAT, CLON, CROT(N), SROT(N))
             IF(LMAP) CALL ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON(N), CLATR, &
                  CLAT, SLAT, CLON, XLON(N), XLAT(N), YLON(N), YLAT(N))
             IF(LAREA) CALL ROT_EQUID_CYLIND_GRID_AREA(CLATR, FILL, AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N,HS,CLON,SLAT,CLAT,SLATR,CLATR,CLONR,RLONR,RLATR) &
       !$OMP& REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             HS=SIGN(1._KD,MOD(RLON(N)-RLON0+180+3600,360._KD)-180)
             CLON=COS((RLON(N)-RLON0)/DPR)
             SLAT=SIN(RLAT(N)/DPR)
             CLAT=COS(RLAT(N)/DPR)
             SLATR=CLAT0*SLAT-SLAT0*CLAT*CLON
             IF(SLATR.LE.-1) THEN
                CLATR=0._KD
                RLONR=0.
                RLATR=-90.
             ELSEIF(SLATR.GE.1) THEN
                CLATR=0._KD
                RLONR=0.
                RLATR=90.
             ELSE
                CLATR=SQRT(1-SLATR**2)
                CLONR=(CLAT0*CLAT*CLON+SLAT0*SLAT)/CLATR
                CLONR=MIN(MAX(CLONR,-1._KD),1._KD)
                RLONR=HS*DPR*ACOS(CLONR)
                RLATR=DPR*ASIN(SLATR)
             ENDIF
             XPTS(N)=(RLONR-WBD)/DLONS+1._KD
             YPTS(N)=(RLATR-SBD)/DLATS+1._KD
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT) CALL ROT_EQUID_CYLIND_VECT_ROT(RLON(N), CLATR, SLATR, &
                     CLAT, SLAT, CLON, CROT(N), SROT(N))
                IF(LMAP) CALL ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON(N), CLATR, &
                     CLAT, SLAT, CLON, XLON(N), XLAT(N), YLON(N), YLAT(N))
                IF(LAREA) CALL ROT_EQUID_CYLIND_GRID_AREA(CLATR, FILL, AREA(N))
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
  END SUBROUTINE GDSWZD_ROT_EQUID_CYLIND
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  ROT_EQUID_CYLIND_ERROR   ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
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
  END SUBROUTINE ROT_EQUID_CYLIND_ERROR
  !
  SUBROUTINE ROT_EQUID_CYLIND_VECT_ROT(RLON, CLATR, SLATR, CLAT, SLAT, &
       CLON, CROT, SROT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  ROT_EQUID_CYLIND_VECT_ROT  VECTOR ROTATION FIELDS FOR
    !                                         ROTATED EQUIDISTANT CYLINDRICAL
    !                                         GRIDS - NON "E" STAGGER.
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE VECTOR ROTATION SINES AND
    !           COSINES FOR A ROTATED EQUIDISTANT CYLINDRICAL GRID -
    !           NON "E" STAGGER.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-07-19  GAYNO    RENAME AS "ROT_EQUID_CYLIND_VECT_ROT."
    ! 2018-07-20  WESLEY   PASS IN CLATR, SLATR, CLAT, SLAT, CLON
    !                      FOR THREADING.
    !
    ! USAGE: CALL ROT_EQUID_CYLIND_VECT_ROT(RLON, CLATR, SLATR, CLAT, SLAT, &
    !                                       CLON, CROT, SROT)
    !
    !   INPUT ARGUMENT LIST:
    !     RLON     - LONGITUDE IN DEGREES (REAL)
    !     CLATR    - COSINE OF ROTATED LATITUDE (REAL)
    !     SLATR    - SINE OF ROTATED LATITUDE (REAL)
    !     CLAT     - COSINE OF LATITUDE (REAL)
    !     SLAT     - SINE OF LATITUDE (REAL)
    !     CLON     - COSINE OF LONGITUDE (REAL)
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

    REAL(KIND=KD),    INTENT(IN   ) :: CLAT, CLATR, CLON, SLAT, SLATR
    REAL         ,    INTENT(IN   ) :: RLON
    REAL         ,    INTENT(  OUT) :: CROT, SROT

    REAL(KIND=KD)                   :: SLON

    IF(IROT.EQ.1) THEN
       IF(CLATR.LE.0._KD) THEN
          CROT=-SIGN(1._KD,SLATR*SLAT0)
          SROT=0.
       ELSE
          SLON=SIN((RLON-RLON0)/DPR)
          CROT=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
          SROT=SLAT0*SLON/CLATR
       ENDIF
    ELSE
       CROT=1.
       SROT=0.
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_VECT_ROT
  !
  SUBROUTINE ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON, CLATR, CLAT, &
       SLAT, CLON, XLON, XLAT, YLON, YLAT)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  ROT_EQUID_CYLIND_MAP_JACOB 
    !     MAP JACOBIANS FOR ROTATED EQUIDISTANT CYLINDRICAL
    !     GRIDS - NON "E" STAGGER.
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE MAP JACOBIANS FOR
    !           A ROTATED EQUIDISTANT CYLINDRICAL GRID -
    !           NON "E" STAGGER.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-09-17  GAYNO    RENAME AS "ROT_EQUID_CYLIND_MAP_JACOB".
    ! 2018-07-20  WESLEY   PASS IN CLATR, CLAT, SLAT, CLON TO ALLOW
    !                      THREADING.
    !
    ! USAGE:  CALL ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON, CLATR, CLAT, &
    !                               SLAT, CLON, XLON, XLAT, YLON, YLAT)
    !
    !   INPUT ARGUMENT LIST:
    !     CLATR    - COSINE OF UNROTATED LATITUDE (REAL)
    !     CLAT     - COSINE OF LATITUDE (REAL)
    !     SLAT     - SINE OF LATITUDE (REAL)
    !     CLON     - COSINE OF LATITUDE (REAL)
    !     FILL     - FILL VALUE FOR UNDEFINED POINTS (REAL)
    !     RLON     - LONGITUDE IN DEGREES (REAL)
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
    !
    IMPLICIT NONE

    REAL(KIND=KD),    INTENT(IN   ) :: CLATR, CLAT, SLAT, CLON
    REAL         ,    INTENT(IN   ) :: FILL, RLON
    REAL         ,    INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    REAL(KIND=KD)                   :: SLON, TERM1, TERM2

    IF(CLATR.LE.0._KD) THEN
       XLON=FILL
       XLAT=FILL
       YLON=FILL
       YLAT=FILL
    ELSE
       SLON=SIN((RLON-RLON0)/DPR)
       TERM1=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
       TERM2=SLAT0*SLON/CLATR
       XLON=TERM1*CLAT/(DLONS*CLATR)
       XLAT=-TERM2/(DLONS*CLATR)
       YLON=TERM2*CLAT/DLATS
       YLAT=TERM1/DLATS
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_MAP_JACOB
  !
  SUBROUTINE ROT_EQUID_CYLIND_GRID_AREA(CLATR, FILL, AREA)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  ROT_EQUID_CYLIND_GRID_AREA  GRID BOX AREA FOR
    !                                          ROTATED EQUIDISTANT CYLINDRICAL
    !                                          GRIDS - NON "E" STAGGER.
    !
    !   PRGMMR: GAYNO     ORG: W/NMC23       DATE: 2015-01-21
    !
    ! ABSTRACT: THIS SUBPROGRAM COMPUTES THE GRID BOX AREA FOR
    !           A ROTATED EQUIDISTANT CYLINDRICAL GRID -
    !           NON "E" STAGGER.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-01-21  GAYNO    INITIAL VERSION
    ! 2015-07-19  GAYNO    RENAME AS "ROT_EQUID_CYLIND_GRID_AREA."
    ! 2018-07-20  WESLEY   PASS IN CLATR FOR THREADING
    !
    ! USAGE:  CALL ROT_EQUID_CYLIND_GRID_AREA(FILL,AREA,CLATR)
    !
    !   INPUT ARGUMENT LIST:
    !     CLATR    - COSINE OF UNROTATED LATITUDE (REAL)
    !     FILL     - FILL VALUE FOR UNDEFINED POINTS (REAL)
    !
    !   OUTPUT ARGUMENT LIST:
    !     AREA     - AREA WEIGHTS IN M**2 (REAL)
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    !
    IMPLICIT NONE

    REAL(KIND=KD),    INTENT(IN   ) :: CLATR
    REAL,             INTENT(IN   ) :: FILL
    REAL,             INTENT(  OUT) :: AREA

    IF(CLATR.LE.0._KD) THEN
       AREA=FILL
    ELSE
       AREA=2._KD*(RERTH**2)*CLATR*(DLONS/DPR)*SIN(0.5_KD*DLATS/DPR)
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_GRID_AREA

end module ip_rot_equid_cylind_grid_mod

