!> @file
!> @brief Interpolate scalar and vector fields with neighbor budget interpolation.
!> @author Mark Iredell @date 96-04-10

!> @brief Interpolate scalar fields (neighbor).
!>
!> @author Mark Iredell @date 96-04-10
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

contains

  !> Interpolate scalar fields (budget).
  !>
  !> This subprogram performs budget interpolation from any grid to
  !> any grid for scalar fields.
  !>
  !> The algorithm simply computes (weighted) averages of neighbor
  !> points arranged in a square box centered around each output grid
  !> point and stretching nearly halfway to each of the neighboring
  !> grid points.
  !>
  !> Options allow choices of number of points in each radius from the
  !> center point (ipopt(1)) which defaults to 2 (if ipopt(1)=-1)
  !> meaning that 25 points will be averaged; further options are the
  !> respective weights for the radius points starting at the center
  !> point (ipopt(2:2+ipopt(1)) which defaults to all 1 (if
  !> ipopt(1)=-1 or ipopt(2)=-1).
  !>
  !> Another option is the minimum percentage for mask, i.e. percent
  !> valid input data required to make output data, (ipopt(3+ipopt(1))
  !> which defaults to 50 (if -1).
  !>
  !> Only horizontal interpolation is performed.
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
  !> latitudes and longitudes are also returned. Input bitmaps will
  !> be interpolated to output bitmaps. Output bitmaps will also be
  !> created when the output grid extends outside of the domain of the
  !> input grid. The output field is set to 0 where the output bitmap
  !> is off.
  !>        
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | Iredell | Initial
  !> 96-10-04 | Iredell | neighbor points not bilinear interpolation
  !> 1999-04-08 | Iredell | split ijkgds into two pieces
  !> 2001-06-18 | Iredell | include minimum mask percentage option
  !> 2015-01-27 | Gayno | replace calls to gdswiz with new merged version of gdswzd.
  !> 2015-07-13 | Gayno | replace grib 1 kgds arrays with grib 2 grid definition template arrays.
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1) is number of
  !> radius points (defaults to 2 if ipopt(1)=-1); ipopt(2:2+ipopt(1))
  !> are respective weights (defaults to all 1 if ipopt(1)=-1 or
  !> ipopt(2)=-1).  ipopt(3+ipopt(1)) is minimum percentage for mask
  !> (defaults to 50 if ipopt(3+ipopt(1)=-1)
  !> @param[in] grid_in The input grid.
  !> @param[in] grid_out The output grid.
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags
  !> @param[in] li (mi,km) input bitmaps (if some ibi(k)=1)
  !> @param[in] gi (mi,km) input fields to interpolate
  !> @param[out] no number of output points
  !> @param[out] rlat (mo) output latitudes in degrees
  !> @param[out] rlon (mo) output longitudes in degrees
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] go (mo,km) output fields interpolated
  !> @param[out] iret return code
  !> - 0    successful interpolation
  !> - 2    unrecognized input grid or no grid overlap
  !> - 3    unrecognized output grid
  !> - 31   invalid undefined output grid
  !> - 32   invalid budget method parameters
  !>
  !> @author Mark Iredell @date 96-04-10
  SUBROUTINE interpolate_neighbor_budget_scalar(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
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


  !> Interpolate vector fields (budget).
  !>
  !> This subprogram performs budget interpolation from any grid to
  !> any grid for vector fields.
  !>  
  !> The algorithm simply computes (weighted) averages of neighbor
  !> points arranged in a square box centered around each output grid
  !> point and stretching nearly halfway to each of the neighboring
  !> grid points.
  !>
  !> Options allow choices of number of points in each radius from the
  !> center point (ipopt(1)) which defaults to 2 (if ipopt(1)=-1)
  !> meaning that 25 points will be averaged; further options are the
  !> respective weights for the radius points starting at the center
  !> point (ipopt(2:2+ipopt(1)) which defaults to all 1 (if
  !> ipopt(1)=-1 or ipopt(2)=-1).
  !>
  !> Another option is the minimum percentage for mask, i.e. percent
  !> valid input data required to make output data, (ipopt(3+ipopt(1))
  !> which defaults to 50 (if -1).
  !>
  !> Only horizontal interpolation is performed.
  !>
  !> The input and output grids are defined by their grib 2 grid
  !> definition template as decoded by the ncep g2 library. the code
  !> recognizes the following projections, where "igdtnumi/o" is the
  !> grib 2 grid defintion template number for the input and output
  !> grids, respectively:
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
  !> description sections.
  !>
  !> As an added bonus the number of output grid points and their
  !> latitudes and longitudes are also returned along with their
  !> vector rotation parameters. Input bitmaps will be interpolated
  !> to output bitmaps.
  !>
  !> Output bitmaps will also be created when the output grid extends
  !> outside of the domain of the input grid. The output field is set
  !> to 0 where the output bitmap is off.
  !>        
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | Iredell | Initial
  !> 1999-04-08 | Iredell | split ijkgds into two pieces
  !> 2001-06-18 | Iredell | include minimum mask percentage option
  !> 2002-01-17 | Iredell | save data from last call for optimization
  !> 2015-01-27 | Gayno | replace calls to gdswiz with new merged routine gdswzd.
  !> 2015-07-13 | Gayno | replace grib 1 kgds arrays with grib 2 grid definition template arrays.
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1) is number of
  !> radius points (defaults to 2 if ipopt(1)=-1); ipopt(2:2+ipopt(1))
  !> are respective weights (defaults to all 1 if ipopt(1)=-1 or
  !> ipopt(2)=-1).  ipopt(3+ipopt(1)) is minimum percentage for mask
  !> (defaults to 50 if ipopt(3+ipopt(1)=-1)
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
  !> @param[out] no number of output points
  !> @param[out] rlat (mo) output latitudes in degrees
  !> @param[out] rlon (mo) output longitudes in degrees
  !> @param[out] crot (mo) vector rotation cosines
  !> @param[out] srot (mo) vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] uo (mo,km) output u-component fields interpolated
  !> @param[out] vo (mo,km) output v-component fields interpolated
  !> @param[out] iret return code
  !> - 0    successful interpolation
  !> - 2    unrecognized input grid or no grid overlap
  !> - 3    unrecognized output grid
  !> - 31   invalid undefined output grid
  !> - 32   invalid budget method parameters
  !>
  !> @author Mark Iredell @date 96-04-10  
  SUBROUTINE interpolate_neighbor_budget_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
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

    ! Save coeffecients between runs and only compute if grid has changed
    INTEGER,                 SAVE  :: MIX=-1
    REAL,        ALLOCATABLE,SAVE  :: CROI(:),SROI(:)
    REAL,        ALLOCATABLE,SAVE  :: XPTI(:),YPTI(:)
    REAL,        ALLOCATABLE,SAVE  :: RLOI(:),RLAI(:)
    class(ip_grid), allocatable, save :: prev_grid_in

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
