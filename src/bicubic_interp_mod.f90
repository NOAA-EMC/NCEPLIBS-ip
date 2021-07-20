!> @file
!! @brief Bicubic interpolation routines for scalars and vectors
!! @author Mark Iredell, Kyle Gerheiser

!> Bicubic interpolation routines for scalars and vectors
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
module bicubic_interp_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  implicit none

  private
  public :: interpolate_bicubic

  interface interpolate_bicubic
     module procedure interpolate_bicubic_scalar
     module procedure interpolate_bicubic_vector
  end interface interpolate_bicubic

contains

  !> @brief This subprogram performs bicubic interpolation
  !! from any grid to any grid for scalar fields.
  !!
  !! @details Bitmaps are now allowed even when invalid points are within
  !! the bicubic template provided the minimum weight is reached.
  !!
  !! Options allow choices between straight bicubic (ipopt(1)=0)
  !! and constrained bicubic (ipopt(1)=1) where the value is
  !! confined within the range of the surrounding 16 points.
  !!
  !! Another option is the minimum percentage for mask,
  !! i.e. percent valid input data required to make output data,
  !! (ipopt(2)) which defaults to 50 (if ipopt(2)=-1).
  !!
  !! Bilinear used within one grid length of boundaries.
  !! Only horizontal interpolation is performed.
  !!
  !! The code recognizes the following projections, where
  !! for the input and output grids, respectively:
  !! as an added bonus the number of output grid points
  !! and their latitudes and longitudes are also returned.
  !! On the other hand, the output can be a set of station points
  !! if igdtnumo<0, in which case the number of points
  !! and their latitudes and longitudes must be input.
  !! output bitmaps will only be created when the output grid
  !! extends outside of the domain of the input grid.
  !! the output field is set to 0 where the output bitmap is off.
  !!
  !! @param[in] ipopt Interpolation options.
  !! - ipopt(1)=0 For straight bicubic;
  !! - ipopt(1)=1 For constrained bicubic where value is confined within the range of the surrounding 4 points.
  !! - ipopt(2) Is minimum percentage for mask (defaults to 50 if ipopt(2)=-1)
  !!
  !! @param[in] grid_in Input grid.
  !! @param[in] grid_out Output grid.
  !! @param[in]  mi Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !! @param[out] mo Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !! @param[in]  km Number of fields to interpolate.
  !! @param[in]  ibi Input bitmap flags.
  !! @param[in]  li Input bitmaps (if some ibi(k)=1).
  !! @param[in]  gi Input fields to interpolate.
  !! @param[in,out] no  Number of output points (only if igdtnumo<0).
  !! @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0).
  !! @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0).
  !! @param[out] ibo Output bitmap flags.
  !! @param[out] lo Output bitmaps (always output).
  !! @param[out] go Output fields interpolated.
  !! @param[out] iret Return code.
  !! - 0 successful interpolation,
  !! - 2 unrecognized input grid or no grid overlap
  !! - 3 unrecognized output grid
  SUBROUTINE interpolate_bicubic_scalar(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,                INTENT(IN   ) :: IPOPT(20)
    INTEGER,                INTENT(IN   ) :: MI,MO,KM
    INTEGER,                INTENT(IN   ) :: IBI(KM)
    INTEGER,                INTENT(INOUT) :: NO
    INTEGER,                INTENT(  OUT) :: IRET, IBO(KM)
    !
    LOGICAL*1,              INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,              INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,                   INTENT(IN   ) :: GI(MI,KM)
    REAL,                   INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,                   INTENT(  OUT) :: GO(MO,KM)
    !
    REAL,                   PARAMETER     :: FILL=-9999.
    !
    INTEGER                               :: IJX(4),IJY(4)
    INTEGER                               :: MCON,MP,N,I,J,K
    INTEGER                               :: NK,NV
    LOGICAL                               :: SAME_GRIDI, SAME_GRIDO
    !
    REAL                                  :: PMP,XIJ,YIJ,XF,YF
    REAL                                  :: G,W,GMIN,GMAX
    REAL                                  :: WX(4),WY(4)
    REAL                                  :: XPTS(MO),YPTS(MO)
    logical :: to_station_points

    ! Save coeffecients between calls and only compute if grids have changed
    REAL,           ALLOCATABLE,SAVE  :: RLATX(:),RLONX(:)
    REAL,           ALLOCATABLE,SAVE  :: WXY(:,:,:)
    INTEGER,                    SAVE  :: NOX=-1,IRETX=-1
    INTEGER,        ALLOCATABLE,SAVE  :: NXY(:,:,:),NC(:)
    class(ip_grid), allocatable,save :: prev_grid_in, prev_grid_out

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IRET=0
    MCON=IPOPT(1)
    MP=IPOPT(2)
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
          CALL GDSWZD(grid_out,0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
          IF(NO.EQ.0) IRET=3
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  ALLOCATE AND SAVE GRID DATA
       IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,NC,NXY,WXY)
          ALLOCATE(RLATX(NO),RLONX(NO),NC(NO),NXY(4,4,NO),WXY(4,4,NO))
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
                IJX(1:4)=FLOOR(XIJ-1)+(/0,1,2,3/)
                IJY(1:4)=FLOOR(YIJ-1)+(/0,1,2,3/)
                XF=XIJ-IJX(2)
                YF=YIJ-IJY(2)
                DO J=1,4
                   DO I=1,4
                      NXY(I,J,N) = grid_in%field_pos(ijx(i), ijy(j))
                   ENDDO
                ENDDO
                IF(MINVAL(NXY(1:4,1:4,N)).GT.0) THEN
                   !  BICUBIC WHERE 16-POINT STENCIL IS AVAILABLE
                   NC(N)=1
                   WX(1)=XF*(1-XF)*(2-XF)/(-6.)
                   WX(2)=(XF+1)*(1-XF)*(2-XF)/2.
                   WX(3)=(XF+1)*XF*(2-XF)/2.
                   WX(4)=(XF+1)*XF*(1-XF)/(-6.)
                   WY(1)=YF*(1-YF)*(2-YF)/(-6.)
                   WY(2)=(YF+1)*(1-YF)*(2-YF)/2.
                   WY(3)=(YF+1)*YF*(2-YF)/2.
                   WY(4)=(YF+1)*YF*(1-YF)/(-6.)
                ELSE
                   !  BILINEAR ELSEWHERE NEAR THE EDGE OF THE GRID
                   NC(N)=2
                   WX(1)=0
                   WX(2)=(1-XF)
                   WX(3)=XF
                   WX(4)=0
                   WY(1)=0
                   WY(2)=(1-YF)
                   WY(3)=YF
                   WY(4)=0
                ENDIF
                DO J=1,4
                   DO I=1,4
                      WXY(I,J,N)=WX(I)*WY(J)
                   ENDDO
                ENDDO
             ELSE
                NC(N)=0
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
       !$OMP PARALLEL DO PRIVATE(NK,K,N,G,W,GMIN,GMAX,J,I) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          IF(NC(N).GT.0) THEN
             G=0
             W=0
             IF(MCON.GT.0) GMIN=HUGE(GMIN)
             IF(MCON.GT.0) GMAX=-HUGE(GMAX)
             DO J=NC(N),5-NC(N)
                DO I=NC(N),5-NC(N)
                   IF(NXY(I,J,N).GT.0)THEN
                      IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K))THEN
                         G=G+WXY(I,J,N)*GI(NXY(I,J,N),K)
                         W=W+WXY(I,J,N)
                         IF(MCON.GT.0) GMIN=MIN(GMIN,GI(NXY(I,J,N),K))
                         IF(MCON.GT.0) GMAX=MAX(GMAX,GI(NXY(I,J,N),K))
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
             LO(N,K)=W.GE.PMP
             IF(LO(N,K)) THEN
                GO(N,K)=G/W
                IF(MCON.GT.0) GO(N,K)=MIN(MAX(GO(N,K),GMIN),GMAX)
             ELSE
                GO(N,K)=0.
             ENDIF
          ELSE
             LO(N,K)=.FALSE.
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
    ELSE
       IF(IRET.EQ.0) IRET=IRETX
       IF(.not. to_station_points) NO=0
    ENDIF
  END SUBROUTINE interpolate_bicubic_scalar

  !> @brief This subprogram performs bicubic interpolation
  !! from any grid to any grid for vector fields.
  !!
  !! @details Bbitmaps are now allowed even when invalid points are within
  !! the bicubic template provided the minimum weight is reached.
  !!
  !! Options allow choices between straight bicubic (ipopt(1)=0)
  !! and constrained bicubic (ipopt(1)=1) where the value is
  !! confined within the range of the surrounding 16 points.
  !!
  !! Another option is the minimum percentage for mask,
  !! i.e. percent valid input data required to make output data,
  !! (ipopt(2)) which defaults to 50 (if ipopt(2)=-1).
  !!
  !! Bilinear used within one grid length of boundaries.
  !! Only horizontal interpolation is performed.
  !!
  !! The code recognizes the following projections, where
  !! for the input and output grids, respectively:
  !! as an added bonus the number of output grid points
  !! and their latitudes and longitudes are also returned.
  !! On the other hand, the output can be a set of station points
  !! if igdtnumo<0, in which case the number of points
  !! and their latitudes and longitudes must be input.
  !!
  !! Output bitmaps will only be created when the output grid
  !! extends outside of the domain of the input grid.
  !! the output field is set to 0 where the output bitmap is off.
  !!
  !! @param[in] ipopt INTEGER (20) INTERPOLATION OPTIONS
  !! - IPOPT(1)=0 FOR STRAIGHT BICUBIC;
  !! - IPOPT(1)=1 FOR CONSTRAINED BICUBIC WHERE VALUE IS CONFINED WITHIN THE RANGE OF THE SURROUNDING 4 POINTS.
  !! - IPOPT(2) IS MINIMUM PERCENTAGE FOR MASK (DEFAULTS TO 50 IF IPOPT(2)=-1)
  !!
  !! @param[in] grid_in Input grid.
  !! @param[in] grid_out Output grid.
  !! @param[in]  mi Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !! @param[out] mo Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !! @param[in]  km Number of fields to interpolate.
  !! @param[in]  ibi Input bitmap flags.
  !! @param[in]  li Input bitmaps (if some ibi(k)=1).
  !! @param[in]  ui Input u-component fields to interpolate.
  !! @param[in]  vi Input v-component fields to interpolate.
  !! @param[in,out] no Number of output points (only if igdtnumo<0).
  !! @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0).
  !! @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0).
  !! @param[in,out] crot Vector rotation cosines (if igdtnumo<0) ugrid=crot*uearth-srot*vearth.
  !! @param[in,out] srot Vector rotation sines (if igdtnumo<0) vgrid=srot*uearth+crot*vearth).
  !! @param[out] ibo Output bitmap flags.
  !! @param[out] lo Output bitmaps (always output).
  !! @param[out] uo Output u-component fields interpolated.
  !! @param[out] vo Output v-component fields interpolated.
  !! @param[out] iret Return code.
  !! - 0 successful interpolation
  !! - 2 unrecognized input grid or no grid overlap
  !! - 3 unrecognized output grid
  subroutine interpolate_bicubic_vector(ipopt, grid_in, grid_out, &
       mi, mo, km, ibi, li, ui, vi, &
       no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,            INTENT(IN   ) :: IPOPT(20)
    INTEGER,            INTENT(IN   ) :: IBI(KM),MI,MO,KM
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
    INTEGER                           :: IJX(4),IJY(4)
    INTEGER                           :: MCON,MP,N,I,J,K,NK,NV
    !
    LOGICAL                           :: SAME_GRIDI,SAME_GRIDO
    !
    REAL                              :: CM,SM,UROT,VROT
    REAL                              :: PMP,XIJ,YIJ,XF,YF
    REAL                              :: U,V,W,UMIN,UMAX,VMIN,VMAX
    REAL                              :: XPTS(MO),YPTS(MO)
    REAL                              :: WX(4),WY(4)
    REAL                              :: XPTI(MI),YPTI(MI),RLOI(MI),RLAI(MI)
    REAL                              :: CROI(MI),SROI(MI)

    logical :: to_station_points
    
    ! Save coeffecients between calls and only compute if grids have changed
    REAL,           ALLOCATABLE, SAVE  :: RLATX(:),RLONX(:),CROTX(:),SROTX(:)
    REAL,           ALLOCATABLE, SAVE  :: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)
    INTEGER,                     SAVE  :: NOX=-1,IRETX=-1
    INTEGER,        ALLOCATABLE, SAVE  :: NXY(:,:,:),NC(:)
    class(ip_grid), allocatable, save :: prev_grid_in, prev_grid_out
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IRET=0
    MCON=IPOPT(1)
    MP=IPOPT(2)
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
       IF(.not. to_station_points) then
          CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT, &
               NO,CROT,SROT)
          IF(NO.EQ.0) IRET=3
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       CALL GDSWZD(grid_in, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI, &
            NV,CROI,SROI)
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  ALLOCATE AND SAVE GRID DATA
       IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,CROTX,SROTX,NC,NXY,WXY,CXY,SXY)
          ALLOCATE(RLATX(NO),RLONX(NO),CROTX(NO),SROTX(NO),NC(NO), &
               NXY(4,4,NO),WXY(4,4,NO),CXY(4,4,NO),SXY(4,4,NO))
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
                IJX(1:4)=FLOOR(XIJ-1)+(/0,1,2,3/)
                IJY(1:4)=FLOOR(YIJ-1)+(/0,1,2,3/)
                XF=XIJ-IJX(2)
                YF=YIJ-IJY(2)
                DO J=1,4
                   DO I=1,4
                      nxy(i,j,n) = grid_in%field_pos(ijx(i), ijy(j))
                   ENDDO
                ENDDO
                IF(MINVAL(NXY(1:4,1:4,N)).GT.0) THEN
                   !  BICUBIC WHERE 16-POINT STENCIL IS AVAILABLE
                   NC(N)=1
                   WX(1)=XF*(1-XF)*(2-XF)/(-6.)
                   WX(2)=(XF+1)*(1-XF)*(2-XF)/2.
                   WX(3)=(XF+1)*XF*(2-XF)/2.
                   WX(4)=(XF+1)*XF*(1-XF)/(-6.)
                   WY(1)=YF*(1-YF)*(2-YF)/(-6.)
                   WY(2)=(YF+1)*(1-YF)*(2-YF)/2.
                   WY(3)=(YF+1)*YF*(2-YF)/2.
                   WY(4)=(YF+1)*YF*(1-YF)/(-6.)
                ELSE
                   !  BILINEAR ELSEWHERE NEAR THE EDGE OF THE GRID
                   NC(N)=2
                   WX(1)=0
                   WX(2)=(1-XF)
                   WX(3)=XF
                   WX(4)=0
                   WY(1)=0
                   WY(2)=(1-YF)
                   WY(3)=YF
                   WY(4)=0
                ENDIF
                DO J=1,4
                   DO I=1,4
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
                NC(N)=0
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
             CROT(N)=CROTX(N)
             SROT(N)=SROTX(N)
          ENDDO
       ENDIF
       !$OMP PARALLEL DO PRIVATE(NK,K,N,U,V,W,UMIN,UMAX,VMIN,VMAX,UROT,VROT,J,I) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          IF(NC(N).GT.0) THEN
             U=0
             V=0
             W=0
             IF(MCON.GT.0) UMIN=HUGE(UMIN)
             IF(MCON.GT.0) UMAX=-HUGE(UMAX)
             IF(MCON.GT.0) VMIN=HUGE(VMIN)
             IF(MCON.GT.0) VMAX=-HUGE(VMAX)
             DO J=NC(N),5-NC(N)
                DO I=NC(N),5-NC(N)
                   IF(NXY(I,J,N).GT.0) THEN
                      IF(IBI(K).EQ.0.OR.LI(NXY(I,J,N),K)) THEN
                         UROT=CXY(I,J,N)*UI(NXY(I,J,N),K)-SXY(I,J,N)*VI(NXY(I,J,N),K)
                         VROT=SXY(I,J,N)*UI(NXY(I,J,N),K)+CXY(I,J,N)*VI(NXY(I,J,N),K)
                         U=U+WXY(I,J,N)*UROT
                         V=V+WXY(I,J,N)*VROT
                         W=W+WXY(I,J,N)
                         IF(MCON.GT.0) UMIN=MIN(UMIN,UROT)
                         IF(MCON.GT.0) UMAX=MAX(UMAX,UROT)
                         IF(MCON.GT.0) VMIN=MIN(VMIN,VROT)
                         IF(MCON.GT.0) VMAX=MAX(VMAX,VROT)
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
                IF(MCON.GT.0) UO(N,K)=MIN(MAX(UO(N,K),UMIN),UMAX)
                IF(MCON.GT.0) VO(N,K)=MIN(MAX(VO(N,K),VMIN),VMAX)
             ELSE
                UO(N,K)=0.
                VO(N,K)=0.
             ENDIF
          ELSE
             LO(N,K)=.FALSE.
             UO(N,K)=0.
             VO(N,K)=0.
          ENDIF
       ENDDO
       DO K=1,KM
          IBO(K)=IBI(K)
          IF(.NOT.ALL(LO(1:NO,K))) IBO(K)=1
       ENDDO
       select type(grid_out)
       type is(ip_equid_cylind_grid)
          CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
       end select
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ELSE
       IF(IRET.EQ.0) IRET=IRETX
       IF(.not. to_station_points) NO=0
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine interpolate_bicubic_vector

end module bicubic_interp_mod
