!> @file
!! @brief Budget interpolation routines for scalars and vectors
!! @author Mark Iredell, Kyle Gerheiser
!! @date July 2021

!> Budget interpolation routines for scalars and vectors.
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
module budget_interp_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
  implicit none

  private
  public :: interpolate_budget

  interface interpolate_budget
     module procedure interpolate_budget_scalar
     module procedure interpolate_budget_vector
  end interface interpolate_budget

contains

  !> Performs budget interpolation
  !! from any grid to any grid (or to random station
  !! points) for scalar fields.
  !!
  !! The algorithm simply computes (weighted) averages of
  !! bilinearly interpolated points arranged in a square box
  !! centered around each output grid point and stretching
  !! nearly halfway to each of the neighboring grid points.
  !!
  !! Options allow choices of number of points in each radius
  !! from the center point (ipopt(1)) which defaults to 2
  !! (if ipopt(1)=-1) meaning that 25 points will be averaged;
  !! further options are the respective weights for the radius
  !! points starting at the center point (ipopt(2:2+ipopt(1))
  !! which defaults to all 1 (if ipopt(1)=-1 or ipopt(2)=-1).
  !!
  !! A special interpolation is done if ipopt(2)=-2.
  !! in this case, the boxes stretch nearly all the way to
  !! each of the neighboring grid points and the weights
  !! are the adjoint of the bilinear interpolation weights.
  !! This case gives quasi-second-order budget interpolation.
  !!
  !! Another option is the minimum percentage for mask,
  !! i.e. percent valid input data required to make output data,
  !! (ipopt(3+ipopt(1)) which defaults to 50 (if -1).
  !!
  !! In cases where there is no or insufficient valid input data,
  !! the user may choose to search for the nearest valid data. 
  !! this is invoked by setting ipopt(20) to the width of 
  !! the search square. The default is 1 (no search). Squares are
  !! searched for valid data in a spiral pattern
  !! starting from the center. No searching is done where
  !! the output grid is outside the input grid.
  !!
  !! Only horizontal interpolation is performed.
  !!
  !! @param[in] ipopt Interpolation options
  !! - ipopt(1) is number of radius points (defaults to 2 if ipopt(1)=-1).
  !! - ipopt(2:2+ipopt(1)) are respective weights (defaults to all 1 if ipopt(1)=-1 or ipopt(2)=-1).
  !! - ipopt(3+ipopt(1)) is minimum percentage for mask (defaults to 50 if ipopt(3+ipopt(1)=-1).
  !! @param[in] grid_in Input grid
  !! @param[in] grid_out Output grid
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
  !! - 0 Successful interpolation.
  !! - 2 Unrecognized input grid or no grid overlap.
  !! - 3 Unrecognized output grid.
  !! - 32 Invalid budget method parameters.
  !!
  !! @author Marke Iredell, George Gayno, Kyle Gerheiser
  !! @date July 2021
  SUBROUTINE interpolate_budget_scalar(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,    INTENT(IN   )     :: IBI(KM), IPOPT(20)
    INTEGER,    INTENT(IN   )     :: KM, MI, MO
    INTEGER,    INTENT(  OUT)     :: IBO(KM), IRET, NO
    !
    LOGICAL*1,  INTENT(IN   )     :: LI(MI,KM)
    LOGICAL*1,  INTENT(  OUT)     :: LO(MO,KM)
    !
    REAL,       INTENT(IN   )     :: GI(MI,KM)
    REAL,       INTENT(INOUT)     :: RLAT(MO), RLON(MO)
    REAL,       INTENT(  OUT)     :: GO(MO,KM)
    !
    REAL,       PARAMETER         :: FILL=-9999.
    !
    INTEGER                       :: IJKGDS1, I1, J1, I2, J2, IB, JB
    INTEGER                       :: IJKGDSA(20), IX, JX, IXS, JXS
    INTEGER                       :: K, KXS, KXT, IGDTNUMO2
    INTEGER                       :: LB, LSW, MP, MSPIRAL, MX
    INTEGER                       :: N, NB, NB1, NB2, NB3, NB4, NV, NX
    INTEGER                       :: N11(MO),N21(MO),N12(MO),N22(MO)
    !
    REAL                          :: GB, LAT(1), LON(1)
    REAL                          :: PMP, RB2, RLOB(MO), RLAB(MO), WB
    REAL                          :: W11(MO), W21(MO), W12(MO), W22(MO)
    REAL                          :: WO(MO,KM), XF, YF, XI, YI, XX, YY
    REAL                          :: XPTS(MO),YPTS(MO),XPTB(MO),YPTB(MO)
    REAL                          :: XXX(1), YYY(1)
    class(ip_grid), allocatable :: grid_out2
    class(ip_grid_descriptor), allocatable :: grid_desc_out2
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    !  DO SUBSECTION OF GRID IF KGDSO(1) IS SUBTRACTED FROM 255.
    IRET=0

    select type(grid_out)
    type is(ip_station_points_grid)
       grid_desc_out2 = grid_out%descriptor
       grid_desc_out2%grid_num = 255 + grid_out%descriptor%grid_num
       grid_out2 = init_grid(grid_desc_out2)

       CALL GDSWZD(grid_out2,-1,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) then
          IRET=3
       end if
    class default
       allocate(grid_out2, source = grid_out)
       CALL GDSWZD(grid_out2, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
    end select

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IF(IPOPT(1).GT.16) IRET=32  
    MSPIRAL=MAX(IPOPT(20),1)
    NB1=IPOPT(1)
    IF(NB1.EQ.-1) NB1=2
    IF(IRET.EQ.0.AND.NB1.LT.0) IRET=32
    LSW=1
    IF(IPOPT(2).EQ.-2) LSW=2
    IF(IPOPT(1).EQ.-1.OR.IPOPT(2).EQ.-1) LSW=0
    IF(IRET.EQ.0.AND.LSW.EQ.1.AND.NB1.GT.15) IRET=32
    MP=IPOPT(3+IPOPT(1))
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01
    IF(IRET.EQ.0) THEN
       NB2=2*NB1+1
       RB2=1./NB2
       NB3=NB2*NB2
       NB4=NB3
       IF(LSW.EQ.2) THEN
          RB2=1./(NB1+1)
          NB4=(NB1+1)**4
       ELSEIF(LSW.EQ.1) THEN
          NB4=IPOPT(2)
          DO IB=1,NB1
             NB4=NB4+8*IB*IPOPT(2+IB)
          ENDDO
       ENDIF
    ELSE
       NB3=0
       NB4=1
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
       IF(LSW.EQ.2) THEN
          WB=(NB1+1-ABS(IB))*(NB1+1-ABS(JB))
       ELSEIF(LSW.EQ.1) THEN
          WB=IPOPT(2+LB)
       ENDIF
       IF(WB.NE.0) THEN
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NO
             XPTB(N)=XPTS(N)+IB*RB2
             YPTB(N)=YPTS(N)+JB*RB2
          ENDDO
          !$OMP END PARALLEL DO
          CALL GDSWZD(grid_out2, 1,NO,FILL,XPTB,YPTB,RLOB,RLAB,NV)
          CALL GDSWZD(grid_in,-1,NO,FILL,XPTB,YPTB,RLOB,RLAB,NV)
          IF(IRET.EQ.0.AND.NV.EQ.0.AND.LB.EQ.0) IRET=2
          !$OMP PARALLEL DO PRIVATE(N,XI,YI,I1,I2,J1,J2,XF,YF) SCHEDULE(STATIC)
          DO N=1,NO
             XI=XPTB(N)
             YI=YPTB(N)
             IF(XI.NE.FILL.AND.YI.NE.FILL) THEN
                I1=XI
                I2=I1+1
                J1=YI
                J2=J1+1
                XF=XI-I1
                YF=YI-J1
                N11(N)=grid_in%field_pos(i1, j1)                
                N21(N)=grid_in%field_pos(i2, j1)
                N12(N)=grid_in%field_pos(i1, j2)
                N22(N)=grid_in%field_pos(i2, j2)
                IF(MIN(N11(N),N21(N),N12(N),N22(N)).GT.0) THEN
                   W11(N)=(1-XF)*(1-YF)
                   W21(N)=XF*(1-YF)
                   W12(N)=(1-XF)*YF
                   W22(N)=XF*YF
                ELSE
                   N11(N)=0
                   N21(N)=0
                   N12(N)=0
                   N22(N)=0
                ENDIF
             ELSE
                N11(N)=0
                N21(N)=0
                N12(N)=0
                N22(N)=0
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !  INTERPOLATE WITH OR WITHOUT BITMAPS
          !$OMP PARALLEL DO PRIVATE(K,N,GB) SCHEDULE(STATIC)
          DO K=1,KM
             DO N=1,NO
                IF(N11(N).GT.0) THEN
                   IF(IBI(K).EQ.0) THEN
                      GB=W11(N)*GI(N11(N),K)+W21(N)*GI(N21(N),K) &
                           +W12(N)*GI(N12(N),K)+W22(N)*GI(N22(N),K)
                      GO(N,K)=GO(N,K)+WB*GB
                      WO(N,K)=WO(N,K)+WB
                   ELSE
                      IF(LI(N11(N),K)) THEN
                         GO(N,K)=GO(N,K)+WB*W11(N)*GI(N11(N),K)
                         WO(N,K)=WO(N,K)+WB*W11(N)
                      ENDIF
                      IF(LI(N21(N),K)) THEN
                         GO(N,K)=GO(N,K)+WB*W21(N)*GI(N21(N),K)
                         WO(N,K)=WO(N,K)+WB*W21(N)
                      ENDIF
                      IF(LI(N12(N),K)) THEN
                         GO(N,K)=GO(N,K)+WB*W12(N)*GI(N12(N),K)
                         WO(N,K)=WO(N,K)+WB*W12(N)
                      ENDIF
                      IF(LI(N22(N),K)) THEN
                         GO(N,K)=GO(N,K)+WB*W22(N)*GI(N22(N),K)
                         WO(N,K)=WO(N,K)+WB*W22(N)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
    ENDDO   ! sub-grid points
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE OUTPUT BITMAPS AND FIELDS
    ! KM is often 1 .. do not do OMP PARALLEL DO here
    KM_LOOP : DO K=1,KM
       IBO(K)=IBI(K)
       !$OMP PARALLEL DO PRIVATE(N,LAT,LON,XXX,YYY,NV,XX,YY,IXS,JXS,MX,KXS,KXT,IX,JX,NX) SCHEDULE(STATIC)
       N_LOOP : DO N=1,NO
          LO(N,K)=WO(N,K).GE.PMP*NB4
          IF(LO(N,K)) THEN
             GO(N,K)=GO(N,K)/WO(N,K)
          ELSEIF (MSPIRAL.GT.1) THEN
             LAT(1)=RLAT(N)
             LON(1)=RLON(N)
             CALL GDSWZD(grid_in,-1,1,FILL,XXX,YYY,LON,LAT,NV)
             XX=XXX(1)
             YY=YYY(1)
             IF(NV.EQ.1)THEN
                I1=NINT(XX)
                J1=NINT(YY)
                IXS=SIGN(1.,XX-I1)
                JXS=SIGN(1.,YY-J1)
                SPIRAL_LOOP : DO MX=2,MSPIRAL**2
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
                      IF(LI(NX,K).OR.IBI(K).EQ.0) THEN
                         GO(N,K)=GI(NX,K)
                         LO(N,K)=.TRUE.
                         CYCLE N_LOOP           
                      ENDIF
                   ENDIF
                ENDDO SPIRAL_LOOP
                IBO(K)=1
                GO(N,K)=0.
             ELSE
                IBO(K)=1
                GO(N,K)=0.
             ENDIF
          ELSE  ! no spiral search option
             IBO(K)=1
             GO(N,K)=0.
          ENDIF
       ENDDO N_LOOP
       !$OMP END PARALLEL DO
    ENDDO KM_LOOP
    
    select type(grid_out2)
    type is(ip_equid_cylind_grid)
       CALL POLFIXS(NO,MO,KM,RLAT,IBO,LO,GO)
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE interpolate_budget_scalar

  
  !> This subprogram performs budget interpolation
  !! from any grid to any grid (or to random station
  !! points) for vector fields.
  !!
  !! The algorithm simply computes (weighted) averages of
  !! bilinearly interpolated points arranged in a square box
  !! centered around each output grid point and stretching
  !! nearly halfway to each of the neighboring grid points.
  !!
  !! Options allow choices of number of points in each radius
  !! from the center point (ipopt(1)) which defaults to 2
  !! (if ipopt(1)=-1) meaning that 25 points will be averaged;
  !! further options are the respective weights for the radius
  !! points starting at the center point (ipopt(2:2+ipopt(1))
  !! which defaults to all 1 (if ipopt(1)=-1 or ipopt(2)=-1).
  !!
  !! A special interpolation is done if ipopt(2)=-2.
  !! in this case, the boxes stretch nearly all the way to
  !! each of the neighboring grid points and the weights
  !! are the adjoint of the bilinear interpolation weights.
  !! This case gives quasi-second-order budget interpolation.
  !!
  !! Another option is the minimum percentage for mask,
  !! i.e. percent valid input data required to make output data,
  !! (ipopt(3+ipopt(1)) which defaults to 50 (if -1).
  !!
  !! In cases where there is no or insufficient valid input data,
  !! the user may choose to search for the nearest valid data. 
  !! this is invoked by setting ipopt(20) to the width of 
  !! the search square. The default is 1 (no search). Squares are
  !! searched for valid data in a spiral pattern
  !! starting from the center. No searching is done where
  !! the output grid is outside the input grid.
  !!
  !! Only horizontal interpolation is performed.
  !!
  !! param[in] ipopt interpolation options
  !! ipopt(1) Number of radius points (defaults to 2 if ipopt(1)=-1);
  !! ipopt(2:2+ipopt(1)) Respective weights (defaults to all 1 if ipopt(1)=-1 or ipopt(2)=-1).
  !! ipopt(3+ipopt(1)) Minimum percentage for mask (defaults to 50 if ipopt(3+ipopt(1)=-1)
  !! @param[in] grid_in Input grid.
  !! @param[in] grid_out Output grid.
  !! @param[in]  mi skip Number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !! @param[out] mo skip Number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !! @param[in]  km Number of fields to interpolate.
  !! @param[in]  ibi Input bitmap flags.
  !! @param[in]  li Input bitmaps (if some ibi(k)=1).
  !! @param[in]  ui Input u-component fields to interpolate.
  !! @param[in]  vi Input v-component fields to interpolate.
  !! @param[in,out] no  Number of output points (only if igdtnumo<0)
  !! @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0)
  !! @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0)
  !! @param[in,out] crot Vector rotation cosines. If interpolating subgrid ugrid=crot * uearth - srot * vearth.
  !! @param[in,out] srot Vector rotation sines. If interpolating subgrid vgrid = srot * uearth + crot * vearth.
  !! @param[out] ibo Output bitmap flags.
  !! @param[out] lo Output bitmaps (always output).
  !! @param[out] uo Output u-component fields interpolated.
  !! @param[out] vo Output v-component fields interpolated.
  !! @param[out] iret Return code.
  !! - 0 Successful interpolation.
  !! - 2 Unrecognized input grid or no grid overlap.
  !! - 3 Unrecognized output grid.
  !! - 32 Invalid budget method parameters.
  !! 
  !! @author Marke Iredell, George Gayno, Kyle Gerheiser
  !! @date July 2021
  SUBROUTINE interpolate_budget_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,          INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,          INTENT(IN   ) :: KM, MI, MO
    INTEGER,          INTENT(  OUT) :: IRET, NO, IBO(KM)
    !
    LOGICAL*1,        INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,             INTENT(  OUT) :: CROT(MO),SROT(MO)
    !
    REAL,             PARAMETER     :: FILL=-9999.
    !
    INTEGER                         :: IGDTNUMO2
    INTEGER                         :: I1,I2,J1,J2,IB,JB,LSW,MP
    INTEGER                         :: K,LB,N,NB,NB1,NB2,NB3,NB4,NV
    INTEGER                         :: N11(MO),N21(MO),N12(MO),N22(MO)
    !
    LOGICAL                         :: SAME_GRID
    !
    REAL                            :: CM11,SM11,CM12,SM12
    REAL                            :: CM21,SM21,CM22,SM22
    REAL                            :: PMP,RB2
    REAL                            :: C11(MO),C21(MO),C12(MO),C22(MO)
    REAL                            :: S11(MO),S21(MO),S12(MO),S22(MO)
    REAL                            :: W11(MO),W21(MO),W12(MO),W22(MO)
    REAL                            :: UB,VB,WB,UROT,VROT
    REAL                            :: U11,V11,U21,V21,U12,V12,U22,V22
    REAL                            :: WI1,WJ1,WI2,WJ2
    REAL                            :: WO(MO,KM),XI,YI
    REAL                            :: XPTS(MO),YPTS(MO)
    REAL                            :: XPTB(MO),YPTB(MO),RLOB(MO),RLAB(MO)

    class(ip_grid_descriptor), allocatable :: desc_out_subgrid
    class(ip_grid), allocatable :: grid_out2
    
    ! Save coeffecients between calls and only compute if grids have changed
    INTEGER,          SAVE          :: MIX=-1
    REAL,         ALLOCATABLE, SAVE :: CROI(:),SROI(:)
    REAL,         ALLOCATABLE, SAVE :: XPTI(:),YPTI(:),RLOI(:),RLAI(:)

    class(ip_grid), allocatable, save :: prev_grid_in

    IRET=0

    ! Negative grid number means interpolate to subgrid
    ! The type of the subgrid is calculated by 255 + 
    select type(grid_out)
    type is(ip_station_points_grid)
       desc_out_subgrid = grid_out%descriptor
       desc_out_subgrid%grid_num = 255 + grid_out%descriptor%grid_num

       grid_out2 = init_grid(desc_out_subgrid)
       CALL GDSWZD(grid_out2,-1,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
    class default
       allocate(grid_out2, source = grid_out)
       CALL GDSWZD(grid_out2, 0,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
    end select

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
       CALL GDSWZD(grid_in, 0,MI,FILL,XPTI,YPTI, &
            RLOI,RLAI,NV,CROI,SROI)
    ENDIF

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    NB1=IPOPT(1)
    IF(NB1.EQ.-1) NB1=2
    IF(IRET.EQ.0.AND.NB1.LT.0) IRET=32
    LSW=1
    IF(IPOPT(2).EQ.-2) LSW=2
    IF(IPOPT(1).EQ.-1.OR.IPOPT(2).EQ.-1) LSW=0
    IF(IRET.EQ.0.AND.LSW.EQ.1.AND.NB1.GT.15) IRET=32
    MP=IPOPT(3+IPOPT(1))
    IF(MP.EQ.-1.OR.MP.EQ.0) MP=50
    IF(MP.LT.0.OR.MP.GT.100) IRET=32
    PMP=MP*0.01
    IF(IRET.EQ.0) THEN
       NB2=2*NB1+1
       RB2=1./NB2
       NB3=NB2*NB2
       NB4=NB3
       IF(LSW.EQ.2) THEN
          RB2=1./(NB1+1)
          NB4=(NB1+1)**4
       ELSEIF(LSW.EQ.1) THEN
          NB4=IPOPT(2)
          DO IB=1,NB1
             NB4=NB4+8*IB*IPOPT(2+IB)
          ENDDO
       ENDIF
    ELSE
       NB3=0
       NB4=1
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
       IF(IPOPT(2).EQ.-2) THEN
          WB=(NB1+1-ABS(IB))*(NB1+1-ABS(JB))
       ELSEIF(IPOPT(2).NE.-1) THEN
          WB=IPOPT(2+LB)
       ENDIF
       IF(WB.NE.0) THEN
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NO
             XPTB(N)=XPTS(N)+IB*RB2
             YPTB(N)=YPTS(N)+JB*RB2
          ENDDO
          !$OMP END PARALLEL DO
          CALL GDSWZD(grid_out2, 1,NO,FILL,XPTB,YPTB, &
               RLOB,RLAB,NV)
          CALL GDSWZD(grid_in,-1,NO,FILL,XPTB,YPTB, &
               RLOB,RLAB,NV)
          IF(IRET.EQ.0.AND.NV.EQ.0.AND.LB.EQ.0) IRET=2
          !$OMP PARALLEL DO PRIVATE(N,XI,YI,I1,I2,WI1,WI2,J1,J2,WJ1,WJ2,CM11,CM21,CM12,CM22,SM11,SM21,SM12,SM22) &
          !$OMP SCHEDULE(STATIC)
          DO N=1,NO
             XI=XPTB(N)
             YI=YPTB(N)
             IF(XI.NE.FILL.AND.YI.NE.FILL) THEN
                I1=XI
                I2=I1+1
                WI2=XI-I1
                WI1=1-WI2
                J1=YI
                J2=J1+1
                WJ2=YI-J1
                WJ1=1-WJ2
                N11(N) = grid_in%field_pos(i1,j1)
                N21(N) = grid_in%field_pos(i2, j1)
                N12(N) = grid_in%field_pos(i1, j2)
                N22(N) = grid_in%field_pos(i2, j2)
                IF(MIN(N11(N),N21(N),N12(N),N22(N)).GT.0) THEN
                   W11(N)=WI1*WJ1
                   W21(N)=WI2*WJ1
                   W12(N)=WI1*WJ2
                   W22(N)=WI2*WJ2
                   CALL MOVECT(RLAI(N11(N)),RLOI(N11(N)),RLAT(N),RLON(N),CM11,SM11)
                   CALL MOVECT(RLAI(N21(N)),RLOI(N21(N)),RLAT(N),RLON(N),CM21,SM21)
                   CALL MOVECT(RLAI(N12(N)),RLOI(N12(N)),RLAT(N),RLON(N),CM12,SM12)
                   CALL MOVECT(RLAI(N22(N)),RLOI(N22(N)),RLAT(N),RLON(N),CM22,SM22)
                   C11(N)=CM11*CROI(N11(N))+SM11*SROI(N11(N))
                   S11(N)=SM11*CROI(N11(N))-CM11*SROI(N11(N))
                   C21(N)=CM21*CROI(N21(N))+SM21*SROI(N21(N))
                   S21(N)=SM21*CROI(N21(N))-CM21*SROI(N21(N))
                   C12(N)=CM12*CROI(N12(N))+SM12*SROI(N12(N))
                   S12(N)=SM12*CROI(N12(N))-CM12*SROI(N12(N))
                   C22(N)=CM22*CROI(N22(N))+SM22*SROI(N22(N))
                   S22(N)=SM22*CROI(N22(N))-CM22*SROI(N22(N))
                ELSE
                   N11(N)=0
                   N21(N)=0
                   N12(N)=0
                   N22(N)=0
                ENDIF
             ELSE
                N11(N)=0
                N21(N)=0
                N12(N)=0
                N22(N)=0
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO 
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          !  INTERPOLATE WITH OR WITHOUT BITMAPS
          !  KM IS OFTEN 1 .. DO NO PUT OMP PARALLEL DO HERE
          DO K=1,KM
             !$OMP PARALLEL DO PRIVATE(N,U11,U12,U21,U22,UB,V11,V12,V21,V22,VB) SCHEDULE(STATIC)
             DO N=1,NO
                IF(N11(N).GT.0) THEN
                   IF(IBI(K).EQ.0) THEN
                      U11=C11(N)*UI(N11(N),K)-S11(N)*VI(N11(N),K)
                      V11=S11(N)*UI(N11(N),K)+C11(N)*VI(N11(N),K)
                      U21=C21(N)*UI(N21(N),K)-S21(N)*VI(N21(N),K)
                      V21=S21(N)*UI(N21(N),K)+C21(N)*VI(N21(N),K)
                      U12=C12(N)*UI(N12(N),K)-S12(N)*VI(N12(N),K)
                      V12=S12(N)*UI(N12(N),K)+C12(N)*VI(N12(N),K)
                      U22=C22(N)*UI(N22(N),K)-S22(N)*VI(N22(N),K)
                      V22=S22(N)*UI(N22(N),K)+C22(N)*VI(N22(N),K)
                      UB=W11(N)*U11+W21(N)*U21+W12(N)*U12+W22(N)*U22
                      VB=W11(N)*V11+W21(N)*V21+W12(N)*V12+W22(N)*V22
                      UO(N,K)=UO(N,K)+WB*UB
                      VO(N,K)=VO(N,K)+WB*VB
                      WO(N,K)=WO(N,K)+WB
                   ELSE
                      IF(LI(N11(N),K)) THEN
                         U11=C11(N)*UI(N11(N),K)-S11(N)*VI(N11(N),K)
                         V11=S11(N)*UI(N11(N),K)+C11(N)*VI(N11(N),K)
                         UO(N,K)=UO(N,K)+WB*W11(N)*U11
                         VO(N,K)=VO(N,K)+WB*W11(N)*V11
                         WO(N,K)=WO(N,K)+WB*W11(N)
                      ENDIF
                      IF(LI(N21(N),K)) THEN
                         U21=C21(N)*UI(N21(N),K)-S21(N)*VI(N21(N),K)
                         V21=S21(N)*UI(N21(N),K)+C21(N)*VI(N21(N),K)
                         UO(N,K)=UO(N,K)+WB*W21(N)*U21
                         VO(N,K)=VO(N,K)+WB*W21(N)*V21
                         WO(N,K)=WO(N,K)+WB*W21(N)
                      ENDIF
                      IF(LI(N12(N),K)) THEN
                         U12=C12(N)*UI(N12(N),K)-S12(N)*VI(N12(N),K)
                         V12=S12(N)*UI(N12(N),K)+C12(N)*VI(N12(N),K)
                         UO(N,K)=UO(N,K)+WB*W12(N)*U12
                         VO(N,K)=VO(N,K)+WB*W12(N)*V12
                         WO(N,K)=WO(N,K)+WB*W12(N)
                      ENDIF
                      IF(LI(N22(N),K)) THEN
                         U22=C22(N)*UI(N22(N),K)-S22(N)*VI(N22(N),K)
                         V22=S22(N)*UI(N22(N),K)+C22(N)*VI(N22(N),K)
                         UO(N,K)=UO(N,K)+WB*W22(N)*U22
                         VO(N,K)=VO(N,K)+WB*W22(N)*V22
                         WO(N,K)=WO(N,K)+WB*W22(N)
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
             !$OMP END PARALLEL DO
          ENDDO
       ENDIF
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE OUTPUT BITMAPS AND FIELDS
    ! KM is often 1, do not put OMP PARALLEL here
    DO K=1,KM
       IBO(K)=IBI(K)
       !$OMP PARALLEL DO PRIVATE(N,UROT,VROT) SCHEDULE(STATIC)
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
       !$OMP END PARALLEL DO
    ENDDO

    select type(grid_out2)
    type is(ip_equid_cylind_grid)
       CALL POLFIXV(NO,MO,KM,RLAT,RLON,IBO,LO,UO,VO)
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE INTERPOLATE_BUDGET_VECTOR

end module budget_interp_mod
