module neighbor_budget_interp_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  implicit none

  private
  public :: interpolate_neighbor_budget

  interface interpolate_neighbor_budget
     module procedure interpolate_neighbor_budget_scalar
     module procedure interpolate_neighbor_budget_vector
  end interface interpolate_neighbor_budget

  INTEGER,                 SAVE  :: MIX=-1
  REAL,        ALLOCATABLE,SAVE  :: CROI(:),SROI(:)
  REAL,        ALLOCATABLE,SAVE  :: XPTI(:),YPTI(:)
  REAL,        ALLOCATABLE,SAVE  :: RLOI(:),RLAI(:)

  class(ip_grid), allocatable :: prev_grid_in

contains

  SUBROUTINE interpolate_neighbor_budget_scalar(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATES6   INTERPOLATE SCALAR FIELDS (BUDGET)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS BUDGET INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
    !           THE ALGORITHM SIMPLY COMPUTES (WEIGHTED) AVERAGES
    !           OF NEIGHBOR POINTS ARRANGED IN A SQUARE BOX
    !           CENTERED AROUND EACH OUTPUT GRID POINT AND STRETCHING
    !           NEARLY HALFWAY TO EACH OF THE NEIGHBORING GRID POINTS.
    !           OPTIONS ALLOW CHOICES OF NUMBER OF POINTS IN EACH RADIUS
    !           FROM THE CENTER POINT (IPOPT(1)) WHICH DEFAULTS TO 2
    !           (IF IPOPT(1)=-1) MEANING THAT 25 POINTS WILL BE AVERAGED;
    !           FURTHER OPTIONS ARE THE RESPECTIVE WEIGHTS FOR THE RADIUS
    !           POINTS STARTING AT THE CENTER POINT (IPOPT(2:2+IPOPT(1))
    !           WHICH DEFAULTS TO ALL 1 (IF IPOPT(1)=-1 OR IPOPT(2)=-1).
    !           ANOTHER OPTION IS THE MINIMUM PERCENTAGE FOR MASK,
    !           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
    !           (IPOPT(3+IPOPT(1)) WHICH DEFAULTS TO 50 (IF -1).
    !           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
    !           THE CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMI/O=01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMI/O=10) MERCATOR CYLINDRICAL
    !             (IGDTNUMI/O=20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMI/O=30) LAMBERT CONFORMAL CONICAL
    !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
    !           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
    !           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
    !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
    !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !   96-10-04  IREDELL  NEIGHBOR POINTS NOT BILINEAR INTERPOLATION
    ! 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
    ! 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      VERSION OF GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATES6(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,GI, &
    !                         NO,RLAT,RLON,IBO,LO,GO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1) IS NUMBER OF RADIUS POINTS
    !                (DEFAULTS TO 2 IF IPOPT(1)=-1);
    !                IPOPT(2:2+IPOPT(1)) ARE RESPECTIVE WEIGHTS
    !                (DEFAULTS TO ALL 1 IF IPOPT(1)=-1 OR IPOPT(2)=-1).
    !                IPOPT(3+IPOPT(1)) IS MINIMUM PERCENTAGE FOR MASK
    !                (DEFAULTS TO 50 IF IPOPT(3+IPOPT(1)=-1)
    !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
    !                  00 - EQUIDISTANT CYLINDRICAL
    !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
    !                       AND NON-"E" STAGGERED
    !                  10 - MERCATOR CYCLINDRICAL
    !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
    !                  30 - LAMBERT CONFORMAL CONICAL
    !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
    !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
    !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATES FOR COMPLETE DEFINITION.
    !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                SAME DEFINITION AS "IGDTNUMI".
    !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
    !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATES FOR COMPLETE DEFINITION.
    !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
    !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
    !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
    !     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
    !     GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES
    !     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
    !     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
    !     GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
    !     IRET     - INTEGER RETURN CODE
    !                0    SUCCESSFUL INTERPOLATION
    !                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
    !                3    UNRECOGNIZED OUTPUT GRID
    !                31   INVALID UNDEFINED OUTPUT GRID
    !                32   INVALID BUDGET METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   GDSWZD       GRID DESCRIPTION SECTION WIZARD
    !   POLFIXS      MAKE MULTIPLE POLE SCALAR VALUES CONSISTENT
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$

    class(ip_grid), intent(in) :: grid_in, grid_out

    INTEGER,        INTENT(IN   ) :: IBI(KM), IPOPT(20), KM, MI, MO
    INTEGER,        INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,      INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,      INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,           INTENT(IN   ) :: GI(MI,KM)
    REAL,           INTENT(  OUT) :: GO(MO,KM), RLAT(MO), RLON(MO)
    !
    REAL,       PARAMETER         :: FILL=-9999.
    !
    INTEGER                       :: IB, I1
    INTEGER                       :: JB, J1, K, LB, LSW, MP, N
    INTEGER                       :: N11(MO), NB, NB1, NB2, NB3, NB4, NV
    !
    REAL                          :: PMP,RLOB(MO),RLAB(MO)
    REAL                          :: WB, WO(MO,KM), XI, YI
    REAL                          :: XPTB(MO),YPTB(MO),XPTS(MO),YPTS(MO)

    logical :: to_station_points

    select type(grid_out)
    type is(ip_station_points_grid)
       to_station_points = .true.
       class default
       to_station_points = .false.
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(.not. to_station_points) THEN
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
    ELSE
       IRET=31
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    NB1=IPOPT(1)
    IF(NB1.EQ.-1) NB1=2
    IF(IRET.EQ.0.AND.NB1.LT.0) IRET=32
    LSW=1
    IF(IPOPT(1).EQ.-1.OR.IPOPT(2).EQ.-1) LSW=0
    IF(IRET.EQ.0.AND.LSW.EQ.1.AND.NB1.GT.15) IRET=32
    MP=IPOPT(3+IPOPT(1))
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01
    IF(IRET.EQ.0) THEN
       NB2=2*NB1+1
       NB3=NB2*NB2
       NB4=NB3
       IF(LSW.EQ.1) THEN
          NB4=IPOPT(2)
          DO IB=1,NB1
             NB4=NB4+8*IB*IPOPT(2+IB)
          ENDDO
       ENDIF
    ELSE
       NB2=0
       NB3=0
       NB4=0
    ENDIF
    DO K=1,KM
       DO N=1,NO
          GO(N,K)=0.
          WO(N,K)=0.
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX

    DO NB=1,NB3
       !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS
       JB=(NB-1)/NB2-NB1
       IB=NB-(JB+NB1)*NB2-NB1-1
       LB=MAX(ABS(IB),ABS(JB))
       WB=1
       IF(LSW.EQ.1) WB=IPOPT(2+LB)
       IF(WB.NE.0) THEN
          DO N=1,NO
             XPTB(N)=XPTS(N)+IB/REAL(NB2)
             YPTB(N)=YPTS(N)+JB/REAL(NB2)
          ENDDO
          CALL GDSWZD(grid_out, 1,NO,FILL,XPTB,YPTB,RLOB,RLAB,NV)
          CALL GDSWZD(grid_in,-1,NO,FILL,XPTB,YPTB,RLOB,RLAB,NV)
          IF(IRET.EQ.0.AND.NV.EQ.0.AND.LB.EQ.0) IRET=2
          DO N=1,NO
             XI=XPTB(N)
             YI=YPTB(N)
             IF(XI.NE.FILL.AND.YI.NE.FILL) THEN
                I1=NINT(XI)
                J1=NINT(YI)
                N11(N)=grid_in%field_pos(i1, j1)
             ELSE
                N11(N)=0
             ENDIF
          ENDDO
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !  INTERPOLATE WITH OR WITHOUT BITMAPS
          DO K=1,KM
             DO N=1,NO
                IF(N11(N).GT.0) THEN
                   IF(IBI(K).EQ.0.OR.LI(N11(N),K)) THEN
                      GO(N,K)=GO(N,K)+WB*GI(N11(N),K)
                      WO(N,K)=WO(N,K)+WB
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE OUTPUT BITMAPS AND FIELDS
    DO K=1,KM
       IBO(K)=IBI(K)
       DO N=1,NO
          LO(N,K)=WO(N,K).GE.PMP*NB4
          IF(LO(N,K)) THEN
             GO(N,K)=GO(N,K)/WO(N,K)
          ELSE
             IBO(K)=1
             GO(N,K)=0.
          ENDIF
       ENDDO
    ENDDO

    select type(grid_out)
    type is(ip_equid_cylind_grid)
       CALL POLFIXS(NO,MO,KM,RLAT,IBO,LO,GO)
    end select

  END SUBROUTINE INTERPOLATE_NEIGHBOR_BUDGET_SCALAR


  SUBROUTINE interpolate_neighbor_budget_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATEV6   INTERPOLATE VECTOR FIELDS (BUDGET)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS BUDGET INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
    !           THE ALGORITHM SIMPLY COMPUTES (WEIGHTED) AVERAGES
    !           OF NEIGHBOR POINTS ARRANGED IN A SQUARE BOX
    !           CENTERED AROUND EACH OUTPUT GRID POINT AND STRETCHING
    !           NEARLY HALFWAY TO EACH OF THE NEIGHBORING GRID POINTS.
    !           OPTIONS ALLOW CHOICES OF NUMBER OF POINTS IN EACH RADIUS
    !           FROM THE CENTER POINT (IPOPT(1)) WHICH DEFAULTS TO 2
    !           (IF IPOPT(1)=-1) MEANING THAT 25 POINTS WILL BE AVERAGED;
    !           FURTHER OPTIONS ARE THE RESPECTIVE WEIGHTS FOR THE RADIUS
    !           POINTS STARTING AT THE CENTER POINT (IPOPT(2:2+IPOPT(1))
    !           WHICH DEFAULTS TO ALL 1 (IF IPOPT(1)=-1 OR IPOPT(2)=-1).
    !           ANOTHER OPTION IS THE MINIMUM PERCENTAGE FOR MASK,
    !           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
    !           (IPOPT(3+IPOPT(1)) WHICH DEFAULTS TO 50 (IF -1).
    !           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
    !
    !           THE INPUT AND OUTPUT GRIDS ARE DEFINED BY THEIR GRIB 2 GRID
    !           DEFINITION TEMPLATE AS DECODED BY THE NCEP G2 LIBRARY.  THE
    !           CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMI/O=01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMI/O=10) MERCATOR CYLINDRICAL
    !             (IGDTNUMI/O=20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMI/O=30) LAMBERT CONFORMAL CONICAL
    !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
    !
    !           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
    !           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
    !           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
    !           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
    !           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
    !
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
    !           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
    !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
    !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 1999-04-08  IREDELL  SPLIT IJKGDS INTO TWO PIECES
    ! 2001-06-18  IREDELL  INCLUDE MINIMUM MASK PERCENTAGE OPTION
    ! 2002-01-17  IREDELL  SAVE DATA FROM LAST CALL FOR OPTIMIZATION
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      ROUTINE GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATEV6(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,UI,VI, &
    !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1) IS NUMBER OF RADIUS POINTS
    !                (DEFAULTS TO 2 IF IPOPT(1)=-1);
    !                IPOPT(2:2+IPOPT(1)) ARE RESPECTIVE WEIGHTS
    !                (DEFAULTS TO ALL 1 IF IPOPT(1)=-1 OR IPOPT(2)=-1).
    !                IPOPT(3+IPOPT(1)) IS MINIMUM PERCENTAGE FOR MASK
    !                (DEFAULTS TO 50 IF IPOPT(3+IPOPT(1)=-1)
    !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
    !                  00 - EQUIDISTANT CYLINDRICAL
    !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
    !                       AND NON-"E" STAGGERED
    !                  10 - MERCATOR CYCLINDRICAL
    !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
    !                  30 - LAMBERT CONFORMAL CONICAL
    !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
    !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
    !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATEV FOR COMPLETE DEFINITION.
    !     IGDTLENI - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - INPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     IGDTNUMO - INTEGER GRID DEFINITION TEMPLATE NUMBER - OUTPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. 
    !                SAME DEFINITION AS "IGDTNUMI"
    !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
    !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATEV FOR COMPLETE DEFINITION.
    !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
    !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
    !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS
    !     LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
    !     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
    !     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES
    !     CROT     - REAL (MO) VECTOR ROTATION COSINES
    !     SROT     - REAL (MO) VECTOR ROTATION SINES
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
    !     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
    !     UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
    !     VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
    !     IRET     - INTEGER RETURN CODE
    !                0    SUCCESSFUL INTERPOLATION
    !                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
    !                3    UNRECOGNIZED OUTPUT GRID
    !                31   INVALID UNDEFINED OUTPUT GRID
    !                32   INVALID BUDGET METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   GDSWZD        GRID DESCRIPTION SECTION WIZARD
    !   MOVECT        MOVE A VECTOR ALONG A GREAT CIRCLE
    !   POLFIXV       MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
    !   CHECK_GRIDS6V DETERMINE IF INPUT OR OUTPUT GRIDS HAVE CHANGED
    !                 BETWEEN CALLS TO THIS ROUTINE.
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    class(ip_grid), intent(in) :: grid_in, grid_out

    INTEGER,         INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,         INTENT(IN   ) :: KM, MI, MO
    INTEGER,         INTENT(  OUT) :: IRET, NO, IBO(KM)
    !
    LOGICAL*1,       INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,       INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,            INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,            INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,            INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,            INTENT(  OUT) :: CROT(MO),SROT(MO)
    !
    REAL,            PARAMETER     :: FILL=-9999.
    !
    INTEGER                        :: N11(MO)
    INTEGER                        :: IB, JB, I1, J1
    INTEGER                        :: K, LB, LSW, MP, N, NV
    INTEGER                        :: NB, NB1, NB2, NB3, NB4
    !
    LOGICAL                        :: SAME_GRID
    !
    REAL                           :: C11(MO),S11(MO)
    REAL                           :: CM11, SM11, PMP
    REAL                           :: U11, V11, UROT, VROT
    REAL                           :: WB, WO(MO,KM), XI, YI
    REAL                           :: RLOB(MO),RLAB(MO)
    REAL                           :: XPTS(MO),YPTS(MO)
    REAL                           :: XPTB(MO),YPTB(MO)

    logical :: to_station_points

    select type(grid_out)
    type is(ip_station_points_grid)
       to_station_points = .true.
       class default
       to_station_points = .false.
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(.not. to_station_points) THEN
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
    ELSE
       IRET=31
    ENDIF

    if (.not. allocated(prev_grid_in)) then
       allocate(prev_grid_in, source = grid_in)

       same_grid = .false.
    else
       same_grid = grid_in == prev_grid_in

       if (.not. same_grid) then
          deallocate(prev_grid_in)
          allocate(prev_grid_in, source = grid_in)
       end if
    end if

    IF(.NOT.SAME_GRID) THEN
       IF(MIX.NE.MI) THEN
          IF(MIX.GE.0) DEALLOCATE(XPTI,YPTI,RLOI,RLAI,CROI,SROI)
          ALLOCATE(XPTI(MI),YPTI(MI),RLOI(MI),RLAI(MI),CROI(MI),SROI(MI))
          MIX=MI
       ENDIF
       CALL GDSWZD(grid_in,0,MI,FILL,XPTI,YPTI, &
            RLOI,RLAI,NV,CROI,SROI)
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    NB1=IPOPT(1)
    IF(NB1.EQ.-1) NB1=2
    IF(IRET.EQ.0.AND.NB1.LT.0) IRET=32
    LSW=1
    IF(IPOPT(1).EQ.-1.OR.IPOPT(2).EQ.-1) LSW=0
    IF(IRET.EQ.0.AND.LSW.EQ.1.AND.NB1.GT.15) IRET=32
    MP=IPOPT(3+IPOPT(1))
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01
    IF(IRET.EQ.0) THEN
       NB2=2*NB1+1
       NB3=NB2*NB2
       NB4=NB3
       IF(LSW.EQ.1) THEN
          NB4=IPOPT(2)
          DO IB=1,NB1
             NB4=NB4+8*IB*IPOPT(2+IB)
          ENDDO
       ENDIF
    ELSE
       NB2=0
       NB3=0
       NB4=0
    ENDIF
    DO K=1,KM
       DO N=1,NO
          UO(N,K)=0
          VO(N,K)=0
          WO(N,K)=0.
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
    DO NB=1,NB3
       !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS AND ROTATIONS
       JB=(NB-1)/NB2-NB1
       IB=NB-(JB+NB1)*NB2-NB1-1
       LB=MAX(ABS(IB),ABS(JB))
       WB=1
       IF(LSW.EQ.1) WB=IPOPT(2+LB)
       IF(WB.NE.0) THEN
          DO N=1,NO
             XPTB(N)=XPTS(N)+IB/REAL(NB2)
             YPTB(N)=YPTS(N)+JB/REAL(NB2)
          ENDDO
          CALL GDSWZD(grid_out, 1,NO,FILL,XPTB,YPTB, &
               RLOB,RLAB,NV)
          CALL GDSWZD(grid_in,-1,NO,FILL,XPTB,YPTB, &
               RLOB,RLAB,NV)
          IF(IRET.EQ.0.AND.NV.EQ.0.AND.LB.EQ.0) IRET=2
          DO N=1,NO
             XI=XPTB(N)
             YI=YPTB(N)
             IF(XI.NE.FILL.AND.YI.NE.FILL) THEN
                I1=NINT(XI)
                J1=NINT(YI)
                N11(N)=grid_in%field_pos(i1, j1)
                IF(N11(N).GT.0) THEN
                   CALL MOVECT(RLAI(N11(N)),RLOI(N11(N)),RLAT(N),RLON(N),CM11,SM11)
                   C11(N)=CM11*CROI(N11(N))+SM11*SROI(N11(N))
                   S11(N)=SM11*CROI(N11(N))-CM11*SROI(N11(N))
                ENDIF
             ELSE
                N11(N)=0
             ENDIF
          ENDDO
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !  INTERPOLATE WITH OR WITHOUT BITMAPS
          DO K=1,KM
             DO N=1,NO
                IF(N11(N).GT.0) THEN
                   IF(IBI(K).EQ.0.OR.LI(N11(N),K)) THEN
                      U11=C11(N)*UI(N11(N),K)-S11(N)*VI(N11(N),K)
                      V11=S11(N)*UI(N11(N),K)+C11(N)*VI(N11(N),K)
                      UO(N,K)=UO(N,K)+WB*U11
                      VO(N,K)=VO(N,K)+WB*V11
                      WO(N,K)=WO(N,K)+WB
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO  ! NB LOOP
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE OUTPUT BITMAPS AND FIELDS
    DO K=1,KM
       IBO(K)=IBI(K)
       DO N=1,NO
          LO(N,K)=WO(N,K).GE.PMP*NB4
          IF(LO(N,K)) THEN
             UO(N,K)=UO(N,K)/WO(N,K)
             VO(N,K)=VO(N,K)/WO(N,K)
             UROT=CROT(N)*UO(N,K)-SROT(N)*VO(N,K)
             VROT=SROT(N)*UO(N,K)+CROT(N)*VO(N,K)
             UO(N,K)=UROT
             VO(N,K)=VROT
          ELSE
             IBO(K)=1
             UO(N,K)=0.
             VO(N,K)=0.
          ENDIF
       ENDDO
    ENDDO

    select type(grid_out)
    type is(ip_equid_cylind_grid)
       CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE INTERPOLATE_NEIGHBOR_BUDGET_VECTOR

end module neighbor_budget_interp_mod
