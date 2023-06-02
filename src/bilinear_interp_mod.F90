!> @file
!> @brief Bilinear interpolation routines for scalars and vectors.
!> @author Mark Iredell, Kyle Gerheiser, Eric Engle

!> @brief Bilinear interpolation routines for scalars and vectors.
!>
!> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
module bilinear_interp_mod
  use gdswzd_mod
  use ip_grids_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
  use polfix_mod
  implicit none

  private
  public :: interpolate_bilinear

  interface interpolate_bilinear
     module procedure interpolate_bilinear_scalar
     module procedure interpolate_bilinear_vector
  end interface interpolate_bilinear

contains

  !> This subprogram performs bilinear interpolation
  !> from any grid to any grid for scalar fields.
  !> 
  !> @details Options allow varying the minimum percentage for mask,
  !> i.e. percent valid input data required to make output data,
  !> (ipopt(1)) which defaults to 50 (if ipopt(1)=-1).
  !> only horizontal interpolation is performed.
  !> if no input data is found near the output point, a spiral
  !> search may be invoked by setting ipopt(2)> 0.
  !> no searching is done if output point is outside the input grid.
  !> as an added bonus the number of output grid points
  !> and their latitudes and longitudes are also returned.
  !> on the other hand, the output can be a set of station points
  !> if igdtnumo<0, in which case the number of points
  !> and their latitudes and longitudes must be input.
  !> input bitmaps will be interpolated to output bitmaps.
  !> output bitmaps will also be created when the output grid
  !> extends outside of the domain of the input grid.
  !>
  !> The output field is set to 0 where the output bitmap is off.
  !> @param[in] ipopt interpolation options
  !> - ipopt(1) is minimum percentage for mask (defaults to 50 if ipopt(1)=-1)
  !> - ipopt(2) is width of square to examine in spiral search (defaults to no search if ipopt(2)=-1)
  !> @param[in] grid_in input grid
  !> @param[in] grid_out output grid
  !> @param[in]  mi skip number between input grid fields if km>1 or dimension of input grid fields if km=1
  !> @param[out] mo skip number between output grid fields if km>1 or dimension of output grid fields if km=1
  !> @param[in]  km number of fields to interpolate
  !> @param[in]  ibi input bitmap flags
  !> @param[in]  li input bitmaps (if some ibi(k)=1)
  !> @param[in]  gi input fields to interpolate
  !> @param[in,out] no  number of output points (only if igdtnumo<0)
  !> @param[in,out] rlat output latitudes in degrees (if igdtnumo<0)
  !> @param[in,out] rlon output longitudes in degrees (if igdtnumo<0)
  !> @param[out] ibo output bitmap flags
  !> @param[out] lo output bitmaps (always output)
  !> @param[out] go output fields interpolated
  !> @param[out] iret return code
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !>
  !> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
  subroutine interpolate_bilinear_scalar(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
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

    ! Save coeffecients between calls and only compute if grids have changed
    INTEGER,                    SAVE  :: NOX=-1,IRETX=-1
    INTEGER,        ALLOCATABLE,SAVE  :: NXY(:,:,:)
    REAL,           ALLOCATABLE,SAVE  :: RLATX(:),RLONX(:)
    REAL,           ALLOCATABLE,SAVE  :: WXY(:,:,:)
    class(ip_grid), allocatable,save :: prev_grid_in, prev_grid_out

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
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
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
    ELSE
       IF(IRET.EQ.0) IRET=IRETX
       IF(.not. to_station_points) NO=0
    ENDIF

  end subroutine interpolate_bilinear_scalar

  !> This subprogram performs bilinear interpolation from any grid to
  !> any grid for vector fields.
  !>
  !> Options allow varying the minimum percentage for mask,
  !> i.e. percent valid input data required to make output data,
  !> (ipopt(1)) which defaults to 50 (if ipopt(1)=-1).
  !>
  !> Only horizontal interpolation is performed.
  !> the input and output vectors are rotated so that they are
  !> either resolved relative to the defined grid
  !> in the direction of increasing x and y coordinates
  !> or resolved relative to easterly and northerly directions,
  !> as designated by their respective grid description sections.
  !>
  !> As an added bonus the number of output grid points
  !> and their latitudes and longitudes are also returned
  !> along with their vector rotation parameters.
  !> on the other hand, the data may be interpolated to a set of
  !> station points if igdtnumo < 0, in which case the number
  !> of points and their latitudes and longitudes must be
  !> input along with their vector rotation parameters.
  !> input bitmaps will be interpolated to output bitmaps.
  !> output bitmaps will also be created when the output grid
  !> extends outside of the domain of the input grid.
  !> the output field is set to 0 where the output bitmap is off.
  !>
  !> @param[in] ipopt interpolation options
  !> - ipopt(1) is minimum percentage for mask (defaults to 50 if ipopt(1)=-1)
  !> @param[in] grid_in Input grid
  !> @param[in] grid_out Output grid
  !> @param[in]  mi skip number between input grid fields if km>1 or dimension of input grid fields if km=1
  !> @param[out] mo skip number between output grid fields if km>1 or dimension of output grid fields if km=1
  !> @param[in]  km number of fields to interpolate
  !> @param[in]  ibi input bitmap flags
  !> @param[in]  li input bitmaps (if some ibi(k)=1)
  !> @param[in]  ui input u-component fields to interpolate
  !> @param[in]  vi input v-component fields to interpolate
  !> @param[in,out] no  number of output points (only if igdtnumo<0)
  !> @param[in,out] rlat output latitudes in degrees (if igdtnumo<0)
  !> @param[in,out] rlon output longitudes in degrees (if igdtnumo<0)
  !> @param[in,out] crot vector rotation cosines (if igdtnumo<0) ugrid=crot*uearth-srot*vearth;
  !> @param[in,out] srot vector rotation sines (if igdtnumo<0) vgrid=srot*uearth+crot*vearth)
  !> @param[out] ibo output bitmap flags
  !> @param[out] lo output bitmaps (always output)
  !> @param[out] uo output u-component fields interpolated
  !> @param[out] vo output v-component fields interpolated
  !> @param[out] iret return code
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !>
  !> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
  SUBROUTINE interpolate_bilinear_vector(ipopt,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
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
    
    ! Save coeffecients between calls and only compute if grids have changed
    INTEGER,                    SAVE  :: NOX=-1,IRETX=-1
    INTEGER,        ALLOCATABLE,SAVE  :: NXY(:,:,:)
    REAL,           ALLOCATABLE,SAVE  :: RLATX(:),RLONX(:)
    REAL,           ALLOCATABLE,SAVE  :: CROTX(:),SROTX(:)
    REAL,           ALLOCATABLE,SAVE  :: WXY(:,:,:),CXY(:,:,:),SXY(:,:,:)
    class(ip_grid), allocatable,save :: prev_grid_in, prev_grid_out

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
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       CALL GDSWZD(grid_in, 0,MI,FILL,XPTI,YPTI,RLOI,RLAI,NV,CROI,SROI)
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
