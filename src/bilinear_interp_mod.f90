!> @file
!! @brief Bilinear interpolation routines for scalars and vectors
!! @author Mark Iredell, Kyle Gerheiser

!> Bilinear interpolation routines for scalars and vectors
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
module bilinear_interp_mod
  use gdswzd_mod
  use ip_grids_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
  use polfix_mod
  implicit none

  private
  public :: interpolate_bilinear

  ! Save coeffecients between calls and only compute if grids have changed
  INTEGER,                    SAVE  :: NOX=-1,IRETX=-1
  INTEGER,        ALLOCATABLE,SAVE  :: NXY(:,:,:)
  REAL,           ALLOCATABLE,SAVE  :: RLATX(:),RLONX(:)
  REAL,           ALLOCATABLE,SAVE  :: CROTX(:),SROTX(:)
  REAL,           ALLOCATABLE,SAVE  :: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)

  class(ip_grid), allocatable :: prev_grid_in, prev_grid_out

  interface interpolate_bilinear
     module procedure interpolate_bilinear_scalar
     module procedure interpolate_bilinear_vector
  end interface interpolate_bilinear

contains

  !> THIS SUBPROGRAM PERFORMS BILINEAR INTERPOLATION
  !! FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
  !! OPTIONS ALLOW VARYING THE MINIMUM PERCENTAGE FOR MASK,
  !! I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
  !! (IPOPT(1)) WHICH DEFAULTS TO 50 (IF IPOPT(1)=-1).
  !! ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
  !! IF NO INPUT DATA IS FOUND NEAR THE OUTPUT POINT, A SPIRAL
  !! SEARCH MAY BE INVOKED BY SETTING IPOPT(2)> 0.
  !! NO SEARCHING IS DONE IF OUTPUT POINT IS OUTSIDE THE INPUT GRID.
  !! AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
  !! AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
  !! ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
  !! IF IGDTNUMO<0, IN WHICH CASE THE NUMBER OF POINTS
  !! AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
  !! INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
  !! OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
  !! EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
  !!
  !! THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
  !! @param[in] ipopt interpolation options
  !! IPOPT(1) IS MINIMUM PERCENTAGE FOR MASK (DEFAULTS TO 50 IF IPOPT(1)=-1)
  !! IPOPT(2) IS WIDTH OF SQUARE TO EXAMINE IN SPIRAL SEARCH (DEFAULTS TO NO SEARCH IF IPOPT(2)=-1)
  !! @param[in] grid_in Input grid
  !! @param[in] grid_out Output grid
  !! @param[in]  MI SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1 OR DIMENSION OF INPUT GRID FIELDS IF KM=1
  !! @param[out] MO SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1 OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
  !! @param[in]  km NUMBER OF FIELDS TO INTERPOLATE
  !! @param[in]  IBI INPUT BITMAP FLAGS
  !! @param[in]  LI INPUT BITMAPS (IF SOME IBI(K)=1)
  !! @param[in]  GI INPUT FIELDS TO INTERPOLATE
  !! @param[in,out] NO  NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO<0)
  !! @param[in,out] RLAT OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO<0)
  !! @param[in,out] RLON OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO<0)
  !! @param[out] IBO OUTPUT BITMAP FLAGS
  !! @param[out] LO OUTPUT BITMAPS (ALWAYS OUTPUT)
  !! @param[out] GO OUTPUT FIELDS INTERPOLATED
  !! @param[out] IRET RETURN CODE
  !! 0 SUCCESSFUL INTERPOLATION, 2 UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP, 3 UNRECOGNIZED OUTPUT GRID
  subroutine interpolate_bilinear_scalar(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,               INTENT(IN   ) :: IPOPT(20)
    INTEGER,               INTENT(IN   ) :: MI,MO,KM
    INTEGER,               INTENT(IN   ) :: IBI(KM)
    INTEGER,               INTENT(INOUT) :: NO
    INTEGER,               INTENT(  OUT) :: IRET, IBO(KM)
    !
    LOGICAL*1,             INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,             INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,                  INTENT(IN   ) :: GI(MI,KM)
    REAL,                  INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT) :: GO(MO,KM)
    !
    REAL,                  PARAMETER     :: FILL=-9999.
    !
    INTEGER                              :: IJX(2),IJY(2)
    INTEGER                              :: MP,N,I,J,K
    INTEGER                              :: NK,NV
    INTEGER                              :: MSPIRAL,I1,J1,IXS,JXS
    INTEGER                              :: MX,KXS,KXT,IX,JX,NX
    !
    LOGICAL                              :: SAME_GRIDI, SAME_GRIDO
    !
    REAL                                 :: WX(2),WY(2)
    REAL                                 :: XPTS(MO),YPTS(MO)
    REAL                                 :: PMP,XIJ,YIJ,XF,YF,G,W

    logical :: to_station_points

    IRET=0
    MP=IPOPT(1)
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01
    MSPIRAL=MAX(IPOPT(2),0)

    if (.not. allocated(prev_grid_in) .or. .not. allocated(prev_grid_out)) then
       allocate(prev_grid_in, source = grid_in)
       allocate(prev_grid_out, source = grid_out)

       same_gridi = .false.
       same_grido = .false.
    else
       same_gridi = grid_in == prev_grid_in
       same_grido = grid_out == prev_grid_out

       if (.not. same_gridi .or. .not. same_grido) then
          deallocate(prev_grid_in)
          deallocate(prev_grid_out)

          allocate(prev_grid_in, source = grid_in)
          allocate(prev_grid_out, source = grid_out)
       end if
    end if

    select type(grid_out)
    type is(ip_station_points_grid)
       to_station_points = .true.
       class default
       to_station_points = .false.
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SAVE OR SKIP WEIGHT COMPUTATION
    IF(IRET==0.AND.(to_station_points.OR..NOT.SAME_GRIDI.OR..NOT.SAME_GRIDO))THEN
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
       IF(.not. to_station_points) THEN
          CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
          IF(NO.EQ.0) IRET=3
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)

       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  ALLOCATE AND SAVE GRID DATA
       IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,NXY,WXY)
          ALLOCATE(RLATX(NO),RLONX(NO),NXY(2,2,NO),WXY(2,2,NO))
          NOX=NO
       ENDIF
       IRETX=IRET
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE WEIGHTS
       IF(IRET.EQ.0) THEN
          !$OMP PARALLEL DO PRIVATE(N,XIJ,YIJ,IJX,IJY,XF,YF,J,I,WX,WY) SCHEDULE(STATIC)
          DO N=1,NO
             RLONX(N)=RLON(N)
             RLATX(N)=RLAT(N)
             XIJ=XPTS(N)
             YIJ=YPTS(N)
             IF(XIJ.NE.FILL.AND.YIJ.NE.FILL) THEN
                IJX(1:2)=FLOOR(XIJ)+(/0,1/)
                IJY(1:2)=FLOOR(YIJ)+(/0,1/)
                XF=XIJ-IJX(1)
                YF=YIJ-IJY(1)
                WX(1)=(1-XF)
                WX(2)=XF
                WY(1)=(1-YF)
                WY(2)=YF
                DO J=1,2
                   DO I=1,2
                      NXY(I,J,N)=grid_in%field_pos(ijx(i), ijy(j))
                      WXY(I,J,N)=WX(I)*WY(J)
                   ENDDO
                ENDDO
             ELSE
                NXY(:,:,N)=0
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE OVER ALL FIELDS
    IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
       IF(.not. to_station_points) THEN
          NO=NOX
          DO N=1,NO
             RLON(N)=RLONX(N)
             RLAT(N)=RLATX(N)
          ENDDO
       ENDIF
       !$OMP PARALLEL DO &
       !$OMP PRIVATE(NK,K,N,G,W,J,I) &
       !$OMP PRIVATE(I1,J1,IXS,JXS,MX,KXS,KXT,IX,JX,NX) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          G=0
          W=0
          DO J=1,2
             DO I=1,2
                IF(NXY(I,J,N).GT.0)THEN
                   IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
                      G=G+WXY(I,J,N)*GI(NXY(I,J,N),K)
                      W=W+WXY(I,J,N)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          LO(N,K)=W.GE.PMP
          IF(LO(N,K)) THEN
             GO(N,K)=G/W
          ELSEIF(MSPIRAL.GT.0.AND.XPTS(N).NE.FILL.AND.YPTS(N).NE.FILL) THEN
             I1=NINT(XPTS(N))
             J1=NINT(YPTS(N))
             IXS=SIGN(1.,XPTS(N)-I1)
             JXS=SIGN(1.,YPTS(N)-J1)
             SPIRAL : DO MX=1,MSPIRAL**2
                KXS=SQRT(4*MX-2.5)
                KXT=MX-(KXS**2/4+1)
                SELECT CASE(MOD(KXS,4))
                CASE(1)
                   IX=I1-IXS*(KXS/4-KXT)
                   JX=J1-JXS*KXS/4
                CASE(2)
                   IX=I1+IXS*(1+KXS/4)
                   JX=J1-JXS*(KXS/4-KXT)
                CASE(3)
                   IX=I1+IXS*(1+KXS/4-KXT)
                   JX=J1+JXS*(1+KXS/4)
                CASE DEFAULT
                   IX=I1-IXS*KXS/4
                   JX=J1+JXS*(KXS/4-KXT)
                END SELECT
                NX=grid_in%field_pos(ix, jx)
                IF(NX.GT.0.)THEN
                   IF(LI(NX,K).OR.IBI(K).EQ.0)THEN
                      GO(N,K)=GI(NX,K)
                      LO(N,K)=.TRUE.
                      EXIT SPIRAL
                   ENDIF
                ENDIF
             ENDDO SPIRAL
             IF(.NOT.LO(N,K))THEN
                IBO(K)=1
                GO(N,K)=0.
             ENDIF
          ELSE
             GO(N,K)=0.
          ENDIF
       ENDDO
       DO K=1,KM
          IBO(K)=IBI(K)
          IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
       ENDDO
       select type(grid_out)
       type is(ip_equid_cylind_grid)
          CALL POLFIXS(NO,MO,KM,RLAT,IBO,LO,GO)
       end select
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ELSE
       IF(IRET.EQ.0) IRET=IRETX
       IF(.not. to_station_points) NO=0
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  end subroutine interpolate_bilinear_scalar

  SUBROUTINE interpolate_bilinear_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATEV0   INTERPOLATE VECTOR FIELDS (BILINEAR)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS BILINEAR INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
    !           OPTIONS ALLOW VARYING THE MINIMUM PERCENTAGE FOR MASK,
    !           I.E. PERCENT VALID INPUT DATA REQUIRED TO MAKE OUTPUT DATA,
    !           (IPOPT(1)) WHICH DEFAULTS TO 50 (IF IPOPT(1)=-1).
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
    !           ON THE OTHER HAND, THE DATA MAY BE INTERPOLATED TO A SET OF
    !           STATION POINTS IF IGDTNUMO<0, IN WHICH CASE THE NUMBER
    !           OF POINTS AND THEIR LATITUDES AND LONGITUDES MUST BE
    !           INPUT ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !
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
    ! 2007-05-22  IREDELL  EXTRAPOLATE UP TO HALF A GRID CELL
    ! 2007-10-30  IREDELL  SAVE WEIGHTS AND THREAD FOR PERFORMANCE
    ! 2012-06-26  GAYNO    FIX OUT-OF-BOUNDS ERROR.  SEE NCEPLIBS
    !                      TICKET #9.
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH MERGED VERSION
    !                      OF GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATEV0(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,UI,VI, &
    !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1) IS MINIMUM PERCENTAGE FOR MASK
    !                (DEFAULTS TO 50 IF IPOPT(1)=-1)
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
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE. IGDTNUMO<0
    !                MEANS INTERPOLATE TO RANDOM STATION POINTS.
    !                OTHERWISE, SAME DEFINITION AS "IGDTNUMI".
    !     IGDTMPLO - INTEGER (IGDTLENO) GRID DEFINITION TEMPLATE ARRAY -
    !                OUTPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
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
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO<0)
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO<0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO<0)
    !     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF IGDTNUMO<0)
    !     SROT     - REAL (MO) VECTOR ROTATION SINES (IF IGDTNUMO<0)
    !                (UGRID=CROT*UEARTH-SROT*VEARTH;
    !                 VGRID=SROT*UEARTH+CROT*VEARTH)
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO>=0)
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     CROT     - REAL (MO) VECTOR ROTATION COSINES (IF IGDTNUMO>=0)
    !     SROT     - REAL (MO) VECTOR ROTATION SINES (IF IGDTNUMO>=0)
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
    !
    ! SUBPROGRAMS CALLED:
    !   GDSWZD        GRID DESCRIPTION SECTION WIZARD
    !   MOVECT        MOVE A VECTOR ALONG A GREAT CIRCLE
    !   POLFIXV       MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
    !   CHECK_GRIDS0V DETERMINE IF INPUT OR OUTPUT GRIDS HAVE CHANGED
    !                 BETWEEN CALLS TO THIS ROUTINE.
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,            INTENT(IN   ) :: IPOPT(20),IBI(KM),MI,MO,KM
    INTEGER,            INTENT(INOUT) :: NO
    INTEGER,            INTENT(  OUT) :: IRET, IBO(KM)
    !
    LOGICAL*1,          INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,          INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,               INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,               INTENT(INOUT) :: RLAT(MO),RLON(MO),CROT(MO),SROT(MO)
    REAL,               INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    !
    REAL,               PARAMETER     :: FILL=-9999.
    !
    INTEGER                           :: IJX(2),IJY(2)
    INTEGER                           :: MP,N,I,J,K,NK,NV
    !
    LOGICAL                           :: SAME_GRIDI, SAME_GRIDO
    !
    REAL                              :: CM,SM,UROT,VROT
    REAL                              :: PMP,XIJ,YIJ,XF,YF,U,V,W
    REAL                              :: XPTS(MO),YPTS(MO)
    REAL                              :: WX(2),WY(2)
    REAL                              :: XPTI(MI),YPTI(MI)
    REAL                              :: RLOI(MI),RLAI(MI)
    REAL                              :: CROI(MI),SROI(MI)

    logical :: to_station_points
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IRET=0
    MP=IPOPT(1)
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01

    if (.not. allocated(prev_grid_in) .or. .not. allocated(prev_grid_out)) then
       allocate(prev_grid_in, source = grid_in)
       allocate(prev_grid_out, source = grid_out)

       same_gridi = .false.
       same_grido = .false.
    else
       same_gridi = grid_in == prev_grid_in
       same_grido = grid_out == prev_grid_out

       if (.not. same_gridi .or. .not. same_grido) then
          deallocate(prev_grid_in)
          deallocate(prev_grid_out)

          allocate(prev_grid_in, source = grid_in)
          allocate(prev_grid_out, source = grid_out)
       end if
    end if

    select type(grid_out)
    type is(ip_station_points_grid)
       to_station_points = .true.
       class default
       to_station_points = .false.
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SAVE OR SKIP WEIGHT COMPUTATION
    IF(IRET.EQ.0.AND.(to_station_points.OR..NOT.SAME_GRIDI.OR..NOT.SAME_GRIDO))THEN
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
       IF(.not. to_station_points) THEN
          CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT, &
               NO,CROT,SROT)
          IF(NO.EQ.0) IRET=3
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       CALL GDSWZD(grid_in, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI,NV,&
            CROI,SROI)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  ALLOCATE AND SAVE GRID DATA
       IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,CROTX,SROTX,NXY,WXY,CXY,SXY)
          ALLOCATE(RLATX(NO),RLONX(NO),CROTX(NO),SROTX(NO), &
               NXY(2,2,NO),WXY(2,2,NO),CXY(2,2,NO),SXY(2,2,NO))
          NOX=NO
       ENDIF
       IRETX=IRET
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE WEIGHTS
       IF(IRET.EQ.0) THEN
          !$OMP PARALLEL DO PRIVATE(N,XIJ,YIJ,IJX,IJY,XF,YF,J,I,WX,WY,CM,SM) SCHEDULE(STATIC)
          DO N=1,NO
             RLONX(N)=RLON(N)
             RLATX(N)=RLAT(N)
             CROTX(N)=CROT(N)
             SROTX(N)=SROT(N)
             XIJ=XPTS(N)
             YIJ=YPTS(N)
             IF(XIJ.NE.FILL.AND.YIJ.NE.FILL) THEN
                IJX(1:2)=FLOOR(XIJ)+(/0,1/)
                IJY(1:2)=FLOOR(YIJ)+(/0,1/)
                XF=XIJ-IJX(1)
                YF=YIJ-IJY(1)
                WX(1)=(1-XF)
                WX(2)=XF
                WY(1)=(1-YF)
                WY(2)=YF
                DO J=1,2
                   DO I=1,2
                      nxy(i, j, n) = grid_in%field_pos(ijx(i), ijy(j))
                      WXY(I,J,N)=WX(I)*WY(J)
                      IF(NXY(I,J,N).GT.0) THEN
                         CALL MOVECT(RLAI(NXY(I,J,N)),RLOI(NXY(I,J,N)), &
                              RLAT(N),RLON(N),CM,SM)
                         CXY(I,J,N)=CM*CROI(NXY(I,J,N))+SM*SROI(NXY(I,J,N))
                         SXY(I,J,N)=SM*CROI(NXY(I,J,N))-CM*SROI(NXY(I,J,N))
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                NXY(:,:,N)=0
             ENDIF
          ENDDO
       ENDIF  ! IS IRET 0?
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE OVER ALL FIELDS
    IF(IRET.EQ.0.AND.IRETX.EQ.0) THEN
       IF(.not. to_station_points) THEN
          NO=NOX
          DO N=1,NO
             RLON(N)=RLONX(N)
             RLAT(N)=RLATX(N)
             CROT(N)=CROTX(N)
             SROT(N)=SROTX(N)
          ENDDO
       ENDIF
       !$OMP PARALLEL DO PRIVATE(NK,K,N,U,V,W,UROT,VROT,J,I) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          U=0
          V=0
          W=0
          DO J=1,2
             DO I=1,2
                IF(NXY(I,J,N).GT.0) THEN
                   IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
                      UROT=CXY(I,J,N)*UI(NXY(I,J,N),K)-SXY(I,J,N)*VI(NXY(I,J,N),K)
                      VROT=SXY(I,J,N)*UI(NXY(I,J,N),K)+CXY(I,J,N)*VI(NXY(I,J,N),K)
                      U=U+WXY(I,J,N)*UROT
                      V=V+WXY(I,J,N)*VROT
                      W=W+WXY(I,J,N)
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          LO(N,K)=W.GE.PMP
          IF(LO(N,K)) THEN
             UROT=CROT(N)*U-SROT(N)*V
             VROT=SROT(N)*U+CROT(N)*V
             UO(N,K)=UROT/W
             VO(N,K)=VROT/W
          ELSE
             UO(N,K)=0.
             VO(N,K)=0.
          ENDIF
       ENDDO  ! NK LOOP
       DO K=1,KM
          IBO(K)=IBI(K)
          IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
       ENDDO

       select type(grid_out)
       type is(ip_equid_cylind_grid)
          CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
       end select

    ELSE
       IF(IRET.EQ.0) IRET=IRETX
       IF(.not. to_station_points) NO=0
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE INTERPOLATE_BILINEAR_VECTOR

end module bilinear_interp_mod
