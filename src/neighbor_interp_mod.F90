!> @file
!> @brief Interpolate scalar and vector fields with neighbor interpolation.
!> @author Mark Iredell @date 96-04-10

!> @brief Interpolate scalar fields (neighbor).
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | Iredell | initial
!> 1999-04-08 | Iredell | split ijkgds into two pieces
!> 2001-06-18 | Iredell | include spiral search option
!> 2006-01-04 | Gayno | minor bug fix
!> 2007-10-30 | Iredell | save weights and thread for performance
!> 2012-06-26 | Gayno | fix out-of-bounds error. see nceplibs ticket #9.
!> 2015-01-27 | Gayno | replace calls to gdswiz() with new merged version of gdswzd().
!> 2015-07-13 | Gayno | replace grib 1 kgds arrays with grib 2 grid definition template arrays.
!>
!> @author Mark Iredell @date 96-04-10
module neighbor_interp_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  implicit none

  private
  public :: interpolate_neighbor

  interface interpolate_neighbor
     module procedure interpolate_neighbor_scalar
     module procedure interpolate_neighbor_vector
  end interface interpolate_neighbor

  ! Smallest positive real value (use for equality comparisons)
  REAL :: TINYREAL=TINY(1.0)

contains

  !> Interpolate scalar fields (neighbor).
  !>
  !> This subprogram performs neighbor interpolation from any grid to
  !> any grid for scalar fields.
  !>
  !> Options allow choosing the width of the grid square (ipopt(1)) to
  !> search for valid data, which defaults to 1 (if ipopt(1)=-1). Odd
  !> width squares are centered on the nearest input grid point; even
  !> width squares are centered on the nearest four input grid points.
  !> Squares are searched for valid data in a spiral pattern starting
  !> from the center. No searching is done where the output grid is
  !> outside the input grid. Only horizontal interpolation is
  !> performed.
  !>
  !> The code recognizes the following projections, where "igdtnumi/o"
  !> is the grib 2 grid defintion template number for the input and
  !> output grids, respectively:
  !> - (igdtnumi/o=00) equidistant cylindrical
  !> - (igdtnumi/o=01) rotated equidistant cylindrical. "e" and non-"e" staggered
  !> - (igdtnumi/o=10) mercator cylindrical
  !> - (igdtnumi/o=20) polar stereographic azimuthal
  !> - (igdtnumi/o=30) lambert conformal conical
  !> - (igdtnumi/o=40) gaussian cylindrical
  !>
  !> As an added bonus the number of output grid points and their
  !> latitudes and longitudes are also returned. On the other hand,
  !> the output can be a set of station points if igdtnumo<0, in which
  !> case the number of points and their latitudes and longitudes must
  !> be input.
  !>
  !> Input bitmaps will be interpolated to output bitmaps.
  !>
  !> Output bitmaps will also be created when the output grid extends
  !> outside of the domain of the input grid. The output field is set
  !> to 0 where the output bitmap is off.
  !>        
  !> @param[in] ipopt (20) interpolation options ipopt(1) is width of
  !> square to examine in spiral search (defaults to 1 if ipopt(1)=-1)
  !> @param[in] grid_in The input grid.
  !> @param[in] grid_out The output grid.
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags
  !> @param[in] li (mi,km) input bitmaps (if some ibi(k)=1)
  !> @param[in] gi (mi,km) input fields to interpolate
  !> @param[inout] no number of output points (only if igdtnumo<0)
  !> @param[inout] rlat (no) output latitudes in degrees (if igdtnumo<0)
  !> @param[inout] rlon (no) output longitudes in degrees (if igdtnumo<0)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] go (mo,km) output fields interpolated
  !> @param[out] iret return code
  !> - 0    successful interpolation
  !> - 2    unrecognized input grid or no grid overlap
  !> - 3    unrecognized output grid
  !>
  !> @author Mark Iredell @date 96-04-10
  !> @author Eric Engle @date 23-05-04
  SUBROUTINE interpolate_neighbor_scalar(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,GI,  &
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
    INTEGER                              :: I1,J1,IXS,JXS
    INTEGER                              :: MSPIRAL,N,K,NK
    INTEGER                              :: NV
    INTEGER                              :: MX,KXS,KXT,IX,JX,NX
    !
    LOGICAL                              :: SAME_GRIDI, SAME_GRIDO
    !
    REAL                                 :: XPTS(MO),YPTS(MO)
    logical :: to_station_points

    INTEGER,                     SAVE :: NOX=-1,IRETX=-1
    INTEGER,        ALLOCATABLE, SAVE :: NXY(:)
    REAL,           ALLOCATABLE, SAVE :: RLATX(:),RLONX(:),XPTSX(:),YPTSX(:)
    class(ip_grid), allocatable, save :: prev_grid_in, prev_grid_out
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IRET=0
    MSPIRAL=MAX(IPOPT(1),1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  LOCATE INPUT POINTS
       CALL GDSWZD(grid_in,-1,NO,FILL,XPTS,YPTS,RLON,RLAT,NV)
       IF(IRET.EQ.0.AND.NV.EQ.0) IRET=2
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  ALLOCATE AND SAVE GRID DATA
       IF(NOX.NE.NO) THEN
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,XPTSX,YPTSX,NXY)
          ALLOCATE(RLATX(NO),RLONX(NO),XPTSX(NO),YPTSX(NO),NXY(NO))
          NOX=NO
       ENDIF
       IRETX=IRET
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE WEIGHTS
       IF(IRET.EQ.0) THEN
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NO
             RLONX(N)=RLON(N)
             RLATX(N)=RLAT(N)
             XPTSX(N)=XPTS(N)
             YPTSX(N)=YPTS(N)
             IF(ABS(XPTS(N)-FILL).GT.TINYREAL.AND.ABS(YPTS(N)-FILL).GT.TINYREAL) THEN
                nxy(n) = grid_in%field_pos(NINT(XPTS(N)), NINT(YPTS(N)))
             ELSE
                NXY(N)=0
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
       DO N=1,NO
          XPTS(N)=XPTSX(N)
          YPTS(N)=YPTSX(N)
       ENDDO
       !$OMP PARALLEL DO PRIVATE(NK,K,N,I1,J1,IXS,JXS,MX,KXS,KXT,IX,JX,NX) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          GO(N,K)=0
          LO(N,K)=.FALSE.
          IF(NXY(N).GT.0) THEN
             IF(IBI(K).EQ.0.OR.LI(NXY(N),K)) THEN
                GO(N,K)=GI(NXY(N),K)
                LO(N,K)=.TRUE.
                ! SPIRAL AROUND UNTIL VALID DATA IS FOUND.
             ELSEIF(MSPIRAL.GT.1) THEN
                I1=NINT(XPTS(N))
                J1=NINT(YPTS(N))
                IXS=INT(SIGN(1.,XPTS(N)-I1))
                JXS=INT(SIGN(1.,YPTS(N)-J1))
                DO MX=2,MSPIRAL**2
                   KXS=INT(SQRT(4*MX-2.5))
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
                   nx = grid_in%field_pos(ix, jx)
                   IF(NX.GT.0) THEN
                      IF(LI(NX,K)) THEN
                         GO(N,K)=GI(NX,K)
                         LO(N,K)=.TRUE.
                         EXIT
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
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
  END SUBROUTINE interpolate_neighbor_scalar

  !> Interpolate vector fields (neighbor).
  !>
  !> This subprogram performs neighbor interpolation
  !> from any grid to any grid for vector fields.
  !>
  !> Options allow choosing the width of the grid square (ipopt(1)) to
  !> search for valid data, which defaults to 1 (IF IPOPT(1)=-1). Odd
  !> width squares are centered on the nearest input grid point; even
  !> width squares are centered on the nearest four input grid points.
  !>
  !> Squares are searched for valid data in a spiral pattern starting
  !> from the center. no searching is done where the output grid is
  !> outside the input grid. Only horizontal interpolation is
  !> performed.
  !>
  !> The input and output grids are defined by their grib 2 grid
  !> definition template as decoded by the ncep g2 library.
  !>
  !> The code recognizes the following projections, where "igdtnumi/o"
  !> is the grib 2 grid defintion template number for the input and
  !> output grids, respectively:
  !> - (igdtnumi/o=00) equidistant cylindrical
  !> - (igdtnumi/o=01) rotated equidistant cylindrical. "e" and non-"e" staggered
  !> - (igdtnumi/o=10) mercator cylindrical
  !> - (igdtnumi/o=20) polar stereographic azimuthal
  !> - (igdtnumi/o=30) lambert conformal conical
  !> - (igdtnumi/o=40) gaussian cylindrical
  !>
  !> The input and output vectors are rotated so that they are either
  !> resolved relative to the defined grid in the direction of
  !> increasing x and y coordinates or resolved relative to easterly
  !> and northerly directions, as designated by their respective grid
  !> definition sections.
  !>
  !> As an added bonus the number of output grid points and their
  !> latitudes and longitudes are also returned along with their
  !> vector rotation parameters. On the other hand, the output can be
  !> a set of station points if igdtnumo<0, in which case the number
  !> of points and their latitudes and longitudes must be input along
  !> with their vector rotation parameters.
  !>
  !> Input bitmaps will be interpolated to output bitmaps. output
  !> bitmaps will also be created when the output grid extends outside
  !> of the domain of the input grid. The output field is set to 0
  !> where the output bitmap is off.
  !>        
  !> @param[in] ipopt (20) interpolation options ipopt(1) is width of
  !> square to examine in spiral search (defaults to 1 if ipopt(1)=-1)
  !> @param[in] grid_in The input grid.
  !> @param[in] grid_out The output grid.
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags
  !> @param[in] li (mi,km) input bitmaps (if some ibi(k)=1)
  !> @param[in] ui (mi,km) input u-component fields to interpolate
  !> @param[in] vi (mi,km) input v-component fields to interpolate
  !> @param[inout] no number of output points (only if igdtnumo>=0)
  !> @param[inout] rlat (mo) output latitudes in degrees (if igdtnumo<0)
  !> @param[inout] rlon (mo) output longitudes in degrees (if igdtnumo<0)
  !> @param[inout] crot (mo) vector rotation cosines (if igdtnumo<0)
  !> @param[inout] srot (mo) vector rotation sines (if igdtnumo<0)
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] uo (mo,km) output u-component fields interpolated
  !> @param[out] vo (mo,km) output v-component fields interpolated
  !> @param[out] iret return code
  !> - 0    successful interpolation
  !> - 2    unrecognized input grid or no grid overlap
  !> - 3    unrecognized output grid
  !>
  !> @author Mark Iredell @date 96-04-10
  !> @author Eric Engle @date 23-05-04
  SUBROUTINE interpolate_neighbor_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)

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
    REAL,               INTENT(INOUT) :: CROT(MO),SROT(MO)
    REAL,               INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,               INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    !
    REAL,               PARAMETER     :: FILL=-9999.
    !
    INTEGER                           :: I1,J1,IXS,JXS,MX
    INTEGER                           :: KXS,KXT,IX,JX,NX
    INTEGER                           :: MSPIRAL,N,K,NK,NV
    !
    LOGICAL                           :: SAME_GRIDI, SAME_GRIDO
    !
    REAL                              :: CX,SX,CM,SM,UROT,VROT
    REAL                              :: XPTS(MO),YPTS(MO)
    REAL                              :: CROI(MI),SROI(MI)
    REAL                              :: XPTI(MI),YPTI(MI),RLOI(MI),RLAI(MI)

    logical :: to_station_points

    INTEGER,                     SAVE :: NOX=-1,IRETX=-1
    INTEGER,        ALLOCATABLE, SAVE :: NXY(:)

    REAL,           ALLOCATABLE, SAVE :: RLATX(:),RLONX(:),XPTSX(:),YPTSX(:)
    REAL,           ALLOCATABLE, SAVE :: CROTX(:),SROTX(:),CXY(:),SXY(:)
    class(ip_grid), allocatable, save :: prev_grid_in, prev_grid_out
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IRET=0
    MSPIRAL=MAX(IPOPT(1),1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
          IF(NOX.GE.0) DEALLOCATE(RLATX,RLONX,XPTSX,YPTSX,CROTX,SROTX,NXY,CXY,SXY)
          ALLOCATE(RLATX(NO),RLONX(NO),XPTSX(NO),YPTSX(NO), &
               CROTX(NO),SROTX(NO),NXY(NO),CXY(NO),SXY(NO))
          NOX=NO
       ENDIF
       IRETX=IRET
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  COMPUTE WEIGHTS
       IF(IRET.EQ.0) THEN
          !$OMP PARALLEL DO PRIVATE(N,CM,SM) SCHEDULE(STATIC)
          DO N=1,NO
             RLONX(N)=RLON(N)
             RLATX(N)=RLAT(N)
             XPTSX(N)=XPTS(N)
             YPTSX(N)=YPTS(N)
             CROTX(N)=CROT(N)
             SROTX(N)=SROT(N)
             IF(ABS(XPTS(N)-FILL).GT.TINYREAL.AND.ABS(YPTS(N)-FILL).GT.TINYREAL) THEN
                nxy(n) = grid_in%field_pos(NINT(XPTS(N)),NINT(YPTS(N)))
                IF(NXY(N).GT.0) THEN
                   CALL MOVECT(RLAI(NXY(N)),RLOI(NXY(N)),RLAT(N),RLON(N),CM,SM)
                   CXY(N)=CM*CROI(NXY(N))+SM*SROI(NXY(N))
                   SXY(N)=SM*CROI(NXY(N))-CM*SROI(NXY(N))
                ENDIF
             ELSE
                NXY(N)=0
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
       DO N=1,NO
          XPTS(N)=XPTSX(N)
          YPTS(N)=YPTSX(N)
       ENDDO
       !$OMP PARALLEL DO &
       !$OMP PRIVATE(NK,K,N,I1,J1,IXS,JXS,MX,KXS,KXT,IX,JX,NX) &
       !$OMP PRIVATE(CM,SM,CX,SX,UROT,VROT) SCHEDULE(STATIC)
       DO NK=1,NO*KM
          K=(NK-1)/NO+1
          N=NK-NO*(K-1)
          UO(N,K)=0
          VO(N,K)=0
          LO(N,K)=.FALSE.
          IF(NXY(N).GT.0) THEN
             IF(IBI(K).EQ.0.OR.LI(NXY(N),K)) THEN
                UROT=CXY(N)*UI(NXY(N),K)-SXY(N)*VI(NXY(N),K)
                VROT=SXY(N)*UI(NXY(N),K)+CXY(N)*VI(NXY(N),K)
                UO(N,K)=CROT(N)*UROT-SROT(N)*VROT
                VO(N,K)=SROT(N)*UROT+CROT(N)*VROT
                LO(N,K)=.TRUE.
                ! SPIRAL AROUND UNTIL VALID DATA IS FOUND.
             ELSEIF(MSPIRAL.GT.1) THEN
                I1=NINT(XPTS(N))
                J1=NINT(YPTS(N))
                IXS=INT(SIGN(1.,XPTS(N)-I1))
                JXS=INT(SIGN(1.,YPTS(N)-J1))
                DO MX=2,MSPIRAL**2
                   KXS=INT(SQRT(4*MX-2.5))
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
                   nx = grid_in%field_pos(ix, jx)
                   IF(NX.GT.0) THEN
                      IF(LI(NX,K)) THEN
                         CALL MOVECT(RLAI(NX),RLOI(NX),RLAT(N),RLON(N),CM,SM)
                         CX=CM*CROI(NX)+SM*SROI(NX)
                         SX=SM*CROI(NX)-CM*SROI(NX)
                         UROT=CX*UI(NX,K)-SX*VI(NX,K)
                         VROT=SX*UI(NX,K)+CX*VI(NX,K)
                         UO(N,K)=CROT(N)*UROT-SROT(N)*VROT
                         VO(N,K)=SROT(N)*UROT+CROT(N)*VROT
                         LO(N,K)=.TRUE.
                         EXIT
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
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
  END SUBROUTINE INTERPOLATE_NEIGHBOR_VECTOR

end module neighbor_interp_mod
