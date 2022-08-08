!> @file
!! @brief Interpolate spectral.
!! @author Mark Iredell @date 96-04-10

!> @brief Interpolate spectral.
!!
!! @author Mark Iredell @date 96-04-10
module spectral_interp_mod
  use gdswzd_mod
  use ip_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
  use earth_radius_mod
  implicit none

  private
  public :: interpolate_spectral

  interface interpolate_spectral
     module procedure interpolate_spectral_scalar
     module procedure interpolate_spectral_vector
  end interface interpolate_spectral

  interface polates4
     module procedure polates4_grib1
     module procedure polates4_grib2
  end interface polates4

  interface polatev4
     module procedure polatev4_grib1
     module procedure polatev4_grib2
  end interface polatev4

contains

  !> Interpolate spectral scalar.
  !>
  !> @param ipopt ???
  !> @param grid_in ???
  !> @param grid_out ???
  !> @param MI ???
  !> @param MO ???
  !> @param KM ???
  !> @param IBI ???
  !> @param GI ???
  !> @param NO ???
  !> @param RLAT ???
  !> @param RLON ???
  !> @param IBO ???
  !> @param LO ???
  !> @param GO ???
  !> @param IRET ???
  !>
  !! @author Mark Iredell @date 96-04-10
  subroutine interpolate_spectral_scalar(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    INTEGER,          INTENT(IN   ) :: IPOPT(20)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,          INTENT(IN   ) :: MI, MO
    INTEGER,          INTENT(IN   ) :: IBI(KM), KM
    INTEGER,          INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: GI(MI,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: GO(MO,KM)


    select type(desc_in => grid_in%descriptor)
    type is(grib1_descriptor)
       select type(desc_out => grid_out%descriptor)
       type is(grib1_descriptor)
          CALL POLATES4(IPOPT,desc_in%gds,desc_out%gds,MI,MO,KM,IBI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
       end select

    type is(grib2_descriptor)
       select type(desc_out => grid_out%descriptor)
       type is(grib2_descriptor)
          CALL POLATES4(IPOPT,desc_in%gdt_num,desc_in%gdt_tmpl,desc_in%gdt_len, &
               desc_out%gdt_num,desc_out%gdt_tmpl,desc_out%gdt_len, &
               MI,MO,KM,IBI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
       end select
    end select
  end subroutine interpolate_spectral_scalar

  !> Interpolate spectral vector.
  !>
  !> @param ipopt ???
  !> @param grid_in ???
  !> @param grid_out ???
  !> @param MI ???
  !> @param MO ???
  !> @param KM ???
  !> @param IBI ???
  !> @param UI ???
  !> @param VI ???
  !> @param NO ???
  !> @param RLAT ???
  !> @param RLON ???
  !> @param CROT ???
  !> @param SROT ???
  !> @param IBO ???
  !> @param LO ???
  !> @param UO ???
  !> @param VO ???
  !> @param IRET ???
  !>
  !! @author Mark Iredell @date 96-04-10
  subroutine interpolate_spectral_vector(IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,          INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,          INTENT(IN   ) :: KM, MI, MO
    INTEGER,          INTENT(  OUT) :: IRET, IBO(KM), NO
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,             INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: CROT(MO),SROT(MO)


    select type(desc_in => grid_in%descriptor)
    type is(grib1_descriptor)
       select type(desc_out => grid_out%descriptor)
       type is(grib1_descriptor)
         CALL polatev4_grib1(IPOPT,desc_in%gds,desc_out%gds,MI,MO,KM,IBI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       end select

    type is(grib2_descriptor)
       select type(desc_out => grid_out%descriptor)
       type is(grib2_descriptor)
           CALL POLATEV4(IPOPT,desc_in%gdt_num,desc_in%gdt_tmpl,desc_in%gdt_len, &
            desc_out%gdt_num,desc_out%gdt_tmpl,desc_out%gdt_len, &
            MI,MO,KM,IBI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
       end select
    end select


  end subroutine interpolate_spectral_vector

  !> Interpolate scalar fields (spectral).
  !>
  !> This subprogram performs spectral interpolation from any grid to
  !> any grid for scalar fields. It requires that the input fields be
  !> uniformly global. Options allow choices between triangular shape
  !> (ipopt(1)=0) and rhomboidal shape (ipopt(1)=1) which has no
  !> default; a second option is the truncation (ipopt(2)) which
  !> defaults to a sensible truncation for the input grid (if
  !> opt(2)=-1).
  !>
  !> @note If the output grid is not found in a special list, then the
  !> transform back to grid is not very fast. This special list
  !> contains global cylindrical grids, polar stereographic grids
  !> centered at the pole and mercator grids.
  !>
  !> Only horizontal interpolation is performed.
  !>
  !> The code recognizes the following projections, where "igdtnumi/o"
  !> is the GRIB2 grid defintion template number for the input and
  !> onutput grids, respectively:
  !> - igdtnumi/o = 00 equidistant cylindrical
  !> - igdtnumo = 01 rotated equidistant cylindrical. "e" and non-"e" staggered
  !> - igdtnumo = 10 mercator cylindrical
  !> - igdtnumo = 20 polar stereographic azimuthal
  !> - igdtnumo = 30 lambert conformal conical
  !> - igdtnumi/o = 40 gaussian cylindrical
  !>
  !> As an added bonus the number of output grid points and their
  !> latitudes and longitudes are also returned. On the other hand,
  !> the output can be a set of station points if igdtnumo < 0, in which
  !> case the number of points and their latitudes and longitudes must
  !> be input. Output bitmaps will not be created.
  !>        
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !>   96-04-10 | Iredell | initial
  !> 2001-06-18 | Iredell | improve detection of special fast transform
  !> 2015-01-27 | Gayno | replace calls to gdswiz with new merged version of gdswzd.
  !> 2015-07-13 | Gayno | convert to grib 2. replace grib 1 kgds arrays with grib 2 grid definition template arrays.
  !>
  !> @param[in] ipopt (20) interpolation options; ipopt(1)=0 for
  !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
  !> number (defaults to sensible if ipopt(2)=-1).
  !> @param[in] igdtnumi grid definition template number - input
  !> grid. Corresponds to the gfld%igdtnum component of the
  !> [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) library
  !> gridmod data structure.
  !> - 00 - equidistant cylindrical
  !> - 01 - rotated equidistant cylindrical.  "e" and non-"e" staggered
  !> - 10 - mercator cyclindrical
  !> - 20 - polar stereographic azimuthal
  !> - 30 - lambert conformal conical
  !> - 40 - gaussian equidistant cyclindrical
  !> @param[in] igdtmpli (igdtleni) grid definition template array -
  !> input grid. corresponds to the gfld%igdtmpl component of the ncep
  !> g2 library gridmod data structure: (section 3 info). See comments
  !> in routine ipolates() for complete definition.
  !> @param[in] igdtleni number of elements of the grid definition
  !> template array - input grid.  corresponds to the gfld%igdtlen
  !> component of the ncep g2 library gridmod data structure.
  !> @param[in] igdtnumo grid definition template number - output
  !> grid. Corresponds to the gfld%igdtnum component of the
  !> [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) library
  !> gridmod data structure. igdtnumo<0 means interpolate to random
  !> station points. Otherwise, same definition as igdtnumi.
  !> @param[in] igdtmplo (igdtleno) grid definition template array -
  !> output grid. Corresponds to the gfld%igdtmpl component of the
  !> ncep g2 library gridmod data structure (section 3 info). See
  !> comments in routine ipolates() for complete definition.
  !> @param[in] igdtleno number of elements of the grid definition
  !> template array - output grid. Corresponds to the gfld%igdtlen
  !> component of the
  !> [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) library
  !> gridmod data structure.
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1
  !> @param[out] km number of fields to interpolate
  !> @param[out] ibi (km) input bitmap flags (must be all 0)
  !> @param[out] gi (mi,km) input fields to interpolate
  !> @param[out] no number of output points (only if igdtnumo>=0)
  !> @param[out] rlat (mo) output latitudes in degrees (if igdtnumo<0)
  !> @param[out] rlon (mo) output longitudes in degrees (if igdtnumo<0)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] go (mo,km) output fields interpolated
  !> @param[out] iret return code
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !> - 41 invalid nonglobal input grid
  !> - 42 invalid spectral method parameters
  !>
  !>! @author Mark Iredell @date 96-04-10
  SUBROUTINE POLATES4_grib2(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    INTEGER,          INTENT(IN   ) :: IGDTNUMI, IGDTLENI
    INTEGER,          INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
    INTEGER,          INTENT(IN   ) :: IGDTNUMO, IGDTLENO
    INTEGER,          INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
    INTEGER,          INTENT(IN   ) :: IPOPT(20)
    INTEGER,          INTENT(IN   ) :: MI, MO
    INTEGER,          INTENT(IN   ) :: IBI(KM), KM
    INTEGER,          INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: GI(MI,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: GO(MO,KM)
    !
    REAL,             PARAMETER     :: FILL=-9999.
    REAL,             PARAMETER     :: PI=3.14159265358979
    REAL,             PARAMETER     :: DPR=180./PI
    !
    INTEGER                         :: IDRTI, IDRTO, IG, JG, IM, JM
    INTEGER                         :: IGO, JGO, IMO, JMO
    INTEGER                         :: ISCAN, JSCAN, NSCAN
    INTEGER                         :: ISCANO, JSCANO, NSCANO
    INTEGER                         :: ISKIPI, JSKIPI, ISCALE
    INTEGER                         :: IMAXI, JMAXI, ISPEC
    INTEGER                         :: IP, IPRIME, IPROJ, IROMB, K
    INTEGER                         :: MAXWV, N, NI, NJ, NPS
    !
    REAL                            :: DE, DR, DY
    REAL                            :: DLAT, DLON, DLATO, DLONO
    REAL                            :: GO2(MO,KM), H, HI, HJ
    REAL                            :: ORIENT, SLAT, RERTH, E2
    REAL                            :: RLAT1, RLON1, RLAT2, RLON2, RLATI
    REAL                            :: XMESH, XP, YP
    REAL                            :: XPTS(MO), YPTS(MO)

    type(grib2_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(igdtnumi, igdtleni, igdtmpli)
    desc_out = init_descriptor(igdtnumo, igdtleno, igdtmplo)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(IGDTNUMO.GE.0) THEN
       !CALL GDSWZD(IGDTNUMO,IGDTMPLO,IGDTLENO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  AFFIRM APPROPRIATE INPUT GRID
    !    LAT/LON OR GAUSSIAN
    !    NO BITMAPS
    !    FULL ZONAL COVERAGE
    !    FULL MERIDIONAL COVERAGE
    IDRTI=IGDTNUMI
    IF(IDRTI==40) IDRTI=4
    IF(IDRTI==0.OR.IDRTI==4)THEN
       IM=IGDTMPLI(8)
       JM=IGDTMPLI(9)
       ISCALE=IGDTMPLI(10)*IGDTMPLI(11)
       IF(ISCALE==0) ISCALE=10**6
       RLON1=FLOAT(IGDTMPLI(13))/FLOAT(ISCALE)
       RLON2=FLOAT(IGDTMPLI(16))/FLOAT(ISCALE)
       ISCAN=MOD(IGDTMPLI(19)/128,2)
       JSCAN=MOD(IGDTMPLI(19)/64,2)
       NSCAN=MOD(IGDTMPLI(19)/32,2)
    ELSE
       IRET=41
    ENDIF
    DO K=1,KM
       IF(IBI(K).NE.0) IRET=41
    ENDDO
    IF(IRET.EQ.0) THEN
       IF(ISCAN.EQ.0) THEN
          DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
       ELSE
          DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
       ENDIF
       IG=NINT(360/ABS(DLON))
       IPRIME=1+MOD(-NINT(RLON1/DLON)+IG,IG)
       IMAXI=IG
       JMAXI=JM
       IF(MOD(IG,2).NE.0.OR.IM.LT.IG) IRET=41
    ENDIF
    IF(IRET.EQ.0.AND.IDRTI.EQ.0) THEN
       ISCALE=IGDTMPLI(10)*IGDTMPLI(11)
       IF(ISCALE==0) ISCALE=10**6
       RLAT1=FLOAT(IGDTMPLI(12))/FLOAT(ISCALE)
       RLAT2=FLOAT(IGDTMPLI(15))/FLOAT(ISCALE)
       DLAT=(RLAT2-RLAT1)/(JM-1)
       JG=NINT(180/ABS(DLAT))
       IF(JM.EQ.JG) IDRTI=256
       IF(JM.NE.JG.AND.JM.NE.JG+1) IRET=41
    ELSEIF(IRET.EQ.0.AND.IDRTI.EQ.4) THEN
       JG=IGDTMPLI(18)*2
       IF(JM.NE.JG) IRET=41
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IF(IRET.EQ.0) THEN
       IROMB=IPOPT(1)
       MAXWV=IPOPT(2)
       IF(MAXWV.EQ.-1) THEN
          IF(IROMB.EQ.0.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)
          IF(IROMB.EQ.1.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.0.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/4
          IF(IROMB.EQ.0.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/4
       ENDIF
       IF((IROMB.NE.0.AND.IROMB.NE.1).OR.MAXWV.LT.0) IRET=42
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE
    IF(IRET.EQ.0) THEN
       IF(NSCAN.EQ.0) THEN
          ISKIPI=1
          JSKIPI=IM
       ELSE
          ISKIPI=JM
          JSKIPI=1
       ENDIF
       IF(ISCAN.EQ.1) ISKIPI=-ISKIPI
       IF(JSCAN.EQ.0) JSKIPI=-JSKIPI
       ISPEC=0
       !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
       IF((IGDTNUMO.EQ.0.OR.IGDTNUMO.EQ.40).AND. &
            MOD(IGDTMPLO(8),2).EQ.0.AND.IGDTMPLO(13).EQ.0.AND.IGDTMPLO(19).EQ.0) THEN
          IDRTO=IGDTNUMO
          IF(IDRTO==40)IDRTO=4
          IMO=IGDTMPLO(8)
          JMO=IGDTMPLO(9)
          ISCALE=IGDTMPLO(10)*IGDTMPLO(11)
          IF(ISCALE==0) ISCALE=10**6
          RLON2=FLOAT(IGDTMPLO(16))/FLOAT(ISCALE)
          DLONO=(MOD(RLON2-1+3600,360.)+1)/(IMO-1)
          IGO=NINT(360/ABS(DLONO))
          IF(IMO.EQ.IGO.AND.IDRTO.EQ.0) THEN
             RLAT1=FLOAT(IGDTMPLO(12))/FLOAT(ISCALE)
             RLAT2=FLOAT(IGDTMPLO(15))/FLOAT(ISCALE)
             DLAT=(RLAT2-RLAT1)/(JMO-1)
             JGO=NINT(180/ABS(DLAT))
             IF(JMO.EQ.JGO) IDRTO=256
             IF(JMO.EQ.JGO.OR.JMO.EQ.JGO+1) ISPEC=1
          ELSEIF(IMO.EQ.IGO.AND.IDRTO.EQ.4) THEN
             JGO=IGDTMPLO(18)*2
             IF(JMO.EQ.JGO) ISPEC=1
          ENDIF
          IF(ISPEC.EQ.1) THEN
             CALL SPTRUN(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMO,JMO, &
                  KM,IPRIME,ISKIPI,JSKIPI,MI,0,0,MO,0,GI,GO)
          ENDIF
          !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
       ELSEIF(IGDTNUMO.EQ.20.AND. &
            IGDTMPLO(8).EQ.IGDTMPLO(9).AND.MOD(IGDTMPLO(8),2).EQ.1.AND. &
            IGDTMPLO(15).EQ.IGDTMPLO(16).AND.IGDTMPLO(18).EQ.64) THEN
          NPS=IGDTMPLO(8)
          RLAT1=FLOAT(IGDTMPLO(10))*1.E-6
          RLON1=FLOAT(IGDTMPLO(11))*1.E-6
          ORIENT=FLOAT(IGDTMPLO(14))*1.E-6
          XMESH=FLOAT(IGDTMPLO(15))*1.E-3
          IPROJ=MOD(IGDTMPLO(17)/128,2)
          IP=(NPS+1)/2
          H=(-1.)**IPROJ
          SLAT=FLOAT(ABS(IGDTMPLO(13)))*1.E-6
          CALL EARTH_RADIUS(IGDTMPLO,IGDTLENO,RERTH,E2)
          DE=(1.+SIN(SLAT/DPR))*RERTH
          DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
          XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/XMESH
          YP=1+COS((RLON1-ORIENT)/DPR)*DR/XMESH
          IF(NINT(XP).EQ.IP.AND.NINT(YP).EQ.IP) THEN
             IF(IPROJ.EQ.0) THEN
                CALL SPTRUNS(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     SLAT,XMESH,ORIENT,GI,GO,GO2)
             ELSE
                CALL SPTRUNS(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     SLAT,XMESH,ORIENT,GI,GO2,GO)
             ENDIF
             ISPEC=1
          ENDIF
          !  SPECIAL CASE OF MERCATOR GRID
       ELSEIF(IGDTNUMO.EQ.10) THEN
          NI=IGDTMPLO(8)
          NJ=IGDTMPLO(9)
          RLAT1=FLOAT(IGDTMPLO(10))*1.0E-6
          RLON1=FLOAT(IGDTMPLO(11))*1.0E-6
          RLON2=FLOAT(IGDTMPLO(15))*1.0E-6
          RLATI=FLOAT(IGDTMPLO(13))*1.0E-6
          ISCANO=MOD(IGDTMPLO(16)/128,2)
          JSCANO=MOD(IGDTMPLO(16)/64,2)
          NSCANO=MOD(IGDTMPLO(16)/32,2)
          DY=FLOAT(IGDTMPLO(19))*1.0E-3
          HI=(-1.)**ISCANO
          HJ=(-1.)**(1-JSCANO)
          CALL EARTH_RADIUS(IGDTMPLO,IGDTLENO,RERTH,E2)
          DLONO=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(NI-1)
          DLATO=HJ*DY/(RERTH*COS(RLATI/DPR))*DPR
          IF(NSCANO.EQ.0) THEN
             CALL SPTRUNM(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NI,NJ, &
                  IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                  RLAT1,RLON1,DLATO,DLONO,GI,GO)
             ISPEC=1
          ENDIF
       ENDIF
       !  GENERAL SLOW CASE
       IF(ISPEC.EQ.0) THEN
          CALL SPTRUNG(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
               IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON,GI,GO)
       ENDIF
       DO K=1,KM
          IBO(K)=0
          DO N=1,NO
             LO(N,K)=.TRUE.
          ENDDO
       ENDDO
    ELSE
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             GO(N,K)=0.
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE POLATES4_GRIB2

  !> Interpolate scalar fields (spectral).
  !>
  !> This subprogram performs spectral interpolation from any grid to
  !> any grid for scalar fields.  It requires that the input fields be
  !> uniformly global.
  !>
  !> Options allow choices between triangular shape (ipopt(1)=0) and
  !> rhomboidal shape (ipopt(1)=1) which has no default; a second
  !> option is the truncation (ipopt(2)) which defaults to a sensible
  !> truncation for the input grid (if opt(2)=-1).
  !>
  !> @note If the output grid is not found in a special list, then the
  !> transform back to grid is not very fast. This special list
  !> contains global cylindrical grids, polar stereographic grids
  !> centered at the pole and mercator grids.
  !>
  !> Only horizontal interpolation is performed. The grids are defined
  !> by their grid description sections (passed in integer form as
  !> decoded by subprogram w3fi63()).
  !>
  !> The current code recognizes the following projections:
  !> - kgds(1) = 000 equidistant cylindrical
  !> - kgds(1) = 001 mercator cylindrical
  !> - kgds(1) = 003 lambert conformal conical
  !> - kgds(1) = 004 gaussian cylindrical (spectral native)
  !> - kgds(1) = 005 polar stereographic azimuthal
  !> - kgds(1) = 203 rotated equidistant cylindrical (e-stagger)
  !> - kgds(1) = 205 rotated equidistant cylindrical (b-stagger)
  !>
  !> Where kgds could be either input kgdsi or output kgdso. As an
  !> added bonus the number of output grid points and their latitudes
  !> and longitudes are also returned. On the other hand, the output
  !> can be a set of station points if kgdso(1)<0, in which case the
  !> number of points and their latitudes and longitudes must be
  !> input. Output bitmaps will not be created.
  !>        
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | Iredell | Initial
  !> 2001-06-18 | Iredell | improve detection of special fast transform
  !> 2015-01-27 | Gayno | replace calls to gdswiz() with new merged version of gdswzd().
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1)=0 for
  !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
  !> number (defaults to sensible if ipopt(2)=-1).
  !> @param[in] kgdsi (200) input gds parameters as decoded by w3fi63
  !> @param[in] kgdso (200) output gds parameters (kgdso(1)<0 implies random station points)
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags (must be all 0)
  !> @param[in] gi (mi,km) input fields to interpolate
  !> @param[out] no number of output points (only if kgdso(1)<0)
  !> @param[out] rlat (no) output latitudes in degrees (if kgdso(1)<0)
  !> @param[out] rlon (no) output longitudes in degrees (if kgdso(1)<0)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] go (mo,km) output fields interpolated
  !> @param[out] iret return code
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !> - 41 invalid nonglobal input grid
  !> - 42 invalid spectral method parameters
  !>
  !> @author Iredell @date 96-04-10
  suBROUTINE POLATES4_grib1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    INTEGER,          INTENT(IN   ) :: IPOPT(20), KGDSI(200)
    INTEGER,          INTENT(IN   ) :: KGDSO(200), MI, MO
    INTEGER,          INTENT(IN   ) :: IBI(KM), KM
    INTEGER,          INTENT(  OUT) :: IBO(KM), IRET
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: GI(MI,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: GO(MO,KM)
    !
    REAL,             PARAMETER     :: FILL=-9999.
    REAL,             PARAMETER     :: RERTH=6.3712E6
    REAL,             PARAMETER     :: PI=3.14159265358979
    REAL,             PARAMETER     :: DPR=180./PI
    !
    INTEGER                         :: IDRTI, IDRTO, IG, JG, IM, JM
    INTEGER                         :: IGO, JGO, IMO, JMO
    INTEGER                         :: ISCAN, JSCAN, NSCAN
    INTEGER                         :: ISCANO, JSCANO, NSCANO
    INTEGER                         :: ISKIPI, JSKIPI
    INTEGER                         :: IMAXI, JMAXI, ISPEC
    INTEGER                         :: IP, IPRIME, IPROJ, IROMB, K
    INTEGER                         :: MAXWV, N, NI, NJ, NPS, NO
    !
    REAL                            :: DE, DR, DY
    REAL                            :: DLAT, DLON, DLATO, DLONO
    REAL                            :: GO2(MO,KM), H, HI, HJ
    REAL                            :: ORIENT
    REAL                            :: RLAT1, RLON1, RLAT2, RLON2, RLATI
    REAL                            :: XMESH, XP, YP
    REAL                            :: XPTS(MO), YPTS(MO)

    type(grib1_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(kgdsi)
    desc_out = init_descriptor(kgdso)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(KGDSO(1).GE.0) THEN
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
       IF(NO.EQ.0) IRET=3
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  AFFIRM APPROPRIATE INPUT GRID
    !    LAT/LON OR GAUSSIAN
    !    NO BITMAPS
    !    FULL ZONAL COVERAGE
    !    FULL MERIDIONAL COVERAGE
    IDRTI=KGDSI(1)
    IM=KGDSI(2)
    JM=KGDSI(3)
    RLON1=KGDSI(5)*1.E-3
    RLON2=KGDSI(8)*1.E-3
    ISCAN=MOD(KGDSI(11)/128,2)
    JSCAN=MOD(KGDSI(11)/64,2)
    NSCAN=MOD(KGDSI(11)/32,2)
    IF(IDRTI.NE.0.AND.IDRTI.NE.4) IRET=41
    DO K=1,KM
       IF(IBI(K).NE.0) IRET=41
    ENDDO
    IF(IRET.EQ.0) THEN
       IF(ISCAN.EQ.0) THEN
          DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
       ELSE
          DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
       ENDIF
       IG=NINT(360/ABS(DLON))
       IPRIME=1+MOD(-NINT(RLON1/DLON)+IG,IG)
       IMAXI=IG
       JMAXI=JM
       IF(MOD(IG,2).NE.0.OR.IM.LT.IG) IRET=41
    ENDIF
    IF(IRET.EQ.0.AND.IDRTI.EQ.0) THEN
       RLAT1=KGDSI(4)*1.E-3
       RLAT2=KGDSI(7)*1.E-3
       DLAT=(RLAT2-RLAT1)/(JM-1)
       JG=NINT(180/ABS(DLAT))
       IF(JM.EQ.JG) IDRTI=256
       IF(JM.NE.JG.AND.JM.NE.JG+1) IRET=41
    ELSEIF(IRET.EQ.0.AND.IDRTI.EQ.4) THEN
       JG=KGDSI(10)*2
       IF(JM.NE.JG) IRET=41
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IF(IRET.EQ.0) THEN
       IROMB=IPOPT(1)
       MAXWV=IPOPT(2)
       IF(MAXWV.EQ.-1) THEN
          IF(IROMB.EQ.0.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)
          IF(IROMB.EQ.1.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.0.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/4
          IF(IROMB.EQ.0.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/4
       ENDIF
       IF((IROMB.NE.0.AND.IROMB.NE.1).OR.MAXWV.LT.0) IRET=42
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE
    IF(IRET.EQ.0) THEN
       IF(NSCAN.EQ.0) THEN
          ISKIPI=1
          JSKIPI=IM
       ELSE
          ISKIPI=JM
          JSKIPI=1
       ENDIF
       IF(ISCAN.EQ.1) ISKIPI=-ISKIPI
       IF(JSCAN.EQ.0) JSKIPI=-JSKIPI
       ISPEC=0
       !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
       IF((KGDSO(1).EQ.0.OR.KGDSO(1).EQ.4).AND. &
            MOD(KGDSO(2),2).EQ.0.AND.KGDSO(5).EQ.0.AND.KGDSO(11).EQ.0) THEN
          IDRTO=KGDSO(1)
          IMO=KGDSO(2)
          JMO=KGDSO(3)
          RLON2=KGDSO(8)*1.E-3
          DLONO=(MOD(RLON2-1+3600,360.)+1)/(IMO-1)
          IGO=NINT(360/ABS(DLONO))
          IF(IMO.EQ.IGO.AND.IDRTO.EQ.0) THEN
             RLAT1=KGDSO(4)*1.E-3
             RLAT2=KGDSO(7)*1.E-3
             DLAT=(RLAT2-RLAT1)/(JMO-1)
             JGO=NINT(180/ABS(DLAT))
             IF(JMO.EQ.JGO) IDRTO=256
             IF(JMO.EQ.JGO.OR.JMO.EQ.JGO+1) ISPEC=1
          ELSEIF(IMO.EQ.IGO.AND.IDRTO.EQ.4) THEN
             JGO=KGDSO(10)*2
             IF(JMO.EQ.JGO) ISPEC=1
          ENDIF
          IF(ISPEC.EQ.1) THEN
             CALL SPTRUN(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMO,JMO, &
                  KM,IPRIME,ISKIPI,JSKIPI,MI,0,0,MO,0,GI,GO)
          ENDIF
          !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
       ELSEIF(KGDSO(1).EQ.5.AND. &
            KGDSO(2).EQ.KGDSO(3).AND.MOD(KGDSO(2),2).EQ.1.AND. &
            KGDSO(8).EQ.KGDSO(9).AND.KGDSO(11).EQ.64) THEN
          NPS=KGDSO(2)
          RLAT1=KGDSO(4)*1.E-3
          RLON1=KGDSO(5)*1.E-3
          ORIENT=KGDSO(7)*1.E-3
          XMESH=KGDSO(8)
          IPROJ=MOD(KGDSO(10)/128,2)
          IP=(NPS+1)/2
          H=(-1.)**IPROJ
          DE=(1.+SIN(60./DPR))*RERTH
          DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
          XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/XMESH
          YP=1+COS((RLON1-ORIENT)/DPR)*DR/XMESH
          IF(NINT(XP).EQ.IP.AND.NINT(YP).EQ.IP) THEN
             IF(IPROJ.EQ.0) THEN
                CALL SPTRUNS(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     60.,XMESH,ORIENT,GI,GO,GO2)
             ELSE
                CALL SPTRUNS(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     60.,XMESH,ORIENT,GI,GO2,GO)
             ENDIF
             ISPEC=1
          ENDIF
          !  SPECIAL CASE OF MERCATOR GRID
       ELSEIF(KGDSO(1).EQ.1) THEN
          NI=KGDSO(2)
          NJ=KGDSO(3)
          RLAT1=KGDSO(4)*1.E-3
          RLON1=KGDSO(5)*1.E-3
          RLON2=KGDSO(8)*1.E-3
          RLATI=KGDSO(9)*1.E-3
          ISCANO=MOD(KGDSO(11)/128,2)
          JSCANO=MOD(KGDSO(11)/64,2)
          NSCANO=MOD(KGDSO(11)/32,2)
          DY=KGDSO(13)
          HI=(-1.)**ISCANO
          HJ=(-1.)**(1-JSCANO)
          DLONO=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(NI-1)
          DLATO=HJ*DY/(RERTH*COS(RLATI/DPR))*DPR
          IF(NSCANO.EQ.0) THEN
             CALL SPTRUNM(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NI,NJ, &
                  IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                  RLAT1,RLON1,DLATO,DLONO,GI,GO)
             ISPEC=1
          ENDIF
       ENDIF
       !  GENERAL SLOW CASE
       IF(ISPEC.EQ.0) THEN
          CALL SPTRUNG(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
               IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON,GI,GO)
       ENDIF
       DO K=1,KM
          IBO(K)=0
          DO N=1,NO
             LO(N,K)=.TRUE.
          ENDDO
       ENDDO
    ELSE
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             GO(N,K)=0.
          ENDDO
       ENDDO
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLATES4_GRIB1


  !> Interpolate vector fields (spectral).
  !>
  !> This subprogram performs spectral interpolation from any grid to
  !> any grid for vector fields. It requires that the input fields be
  !> uniformly global. Options allow choices between triangular shape
  !> (ipopt(1)=0) and rhomboidal shape (ipopt(1)=1) which has no
  !> default; a second option is the truncation (ipopt(2)) which
  !> defaults to a sensible truncation for the input grid (if
  !> opt(2)=-1).
  !>
  !> @note If the output grid is not found in a special list, then the
  !> transform back to grid is not very fast.  This special list
  !> contains global cylindrical grids, polar stereographic grids
  !> centered at the pole and mercator grids. Only horizontal
  !> interpolation is performed.
  !>
  !> The input and output grids are defined by their grib 2 grid
  !> definition template as decoded by the ncep g2 library. The code
  !> recognizes the following projections, where "igdtnumi/o" is the
  !> grib 2 grid defintion template number for the input and output
  !> grids, respectively:
  !> - igdtnumi/o=00 equidistant cylindrical
  !> - igdtnumo  =01 rotated equidistant cylindrical. "e" and non-"e" staggered
  !> - igdtnumo  =10 mercator cylindrical
  !> - igdtnumo  =20 polar stereographic azimuthal
  !> - igdtnumo  =30 lambert conformal conical
  !> - igdtnumi/o=40 gaussian cylindrical
  !>
  !> The input and output vectors are rotated so that they are either
  !> resolved relative to the defined grid in the direction of
  !> increasing x and y coordinates or resolved relative to easterly
  !> and northerly directions, as designated by their respective grid
  !> description sections.
  !>
  !> As an added bonus the number of output grid points and their
  !> latitudes and longitudes are also returned along with their
  !> vector rotation parameters.  On the other hand, the output can be
  !> a set of station points if igdtnumo<0, in which case the number
  !> of points and their latitudes and longitudes must be input along
  !> with their vector rotation parameters.
  !>
  !> Output bitmaps will only be created when the output grid extends
  !> outside of the domain of the input grid.  the output field is set
  !> to 0 where the output bitmap is off.
  !>        
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | iredell | initial
  !> 2001-06-18 | iredell | improve detection of special fast transform
  !> 2015-01-27 | gayno | replace calls to gdswiz() with new merged routine gdswzd().
  !> 2015-07-13 | gayno | convert to grib 2. replace grib 1 kgds arrays with grib 2 grid definition template arrays.
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1)=0 for
  !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
  !> number (defaults to sensible if ipopt(2)=-1).
  !> @param[in] igdtnumi grid definition template number - input
  !> grid. Corresponds to the gfld%igdtnum component of the ncep g2
  !> library gridmod data structure:
  !> - 00 equidistant cylindrical
  !> - 01 rotated equidistant cylindrical.  "e" and non-"e" staggered
  !> - 10 mercator cyclindrical
  !> - 20 polar stereographic azimuthal
  !> - 30 lambert conformal conical
  !> - 40 gaussian equidistant cyclindrical
  !> @param[in] igdtmpli (igdtleni) grid definition template array - input
  !> grid. corresponds to the gfld%igdtmpl component of the ncep g2
  !> library gridmod data structure (section 3 info).  see comments in
  !> routine ipolatev for complete definition.
  !> @param[in] igdtleni number of elements of the grid definition
  !> template array - input grid.  corresponds to the gfld%igdtlen
  !> component of the ncep g2 library gridmod data structure.
  !> @param[in] igdtnumo grid definition template number - output
  !> grid. Corresponds to the gfld%igdtnum component of the ncep g2
  !> library gridmod data structure. igdtnumo<0 means interpolate to
  !> random station points. Otherwise, same definition as "igdtnumi".
  !> @param[in] igdtmplo (igdtleno) grid definition template array -
  !> output grid. corresponds to the gfld%igdtmpl component of the
  !> ncep g2 library gridmod data structure (section 3 info).  see
  !> comments in routine ipolatev() for complete definition.
  !> @param[in] igdtleno number of elements of the grid definition
  !> template array - output grid. Corresponds to the gfld%igdtlen
  !> component of the ncep g2 library gridmod data structure.
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1.
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags (must be all 0)
  !> @param[in] ui (mi,km) input u-component fields to interpolate
  !> @param[in] vi (mi,km) input v-component fields to interpolate
  !> @param[out] no number of output points (only if igdtnumo>=0)
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
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !> - 41 invalid nonglobal input grid
  !> - 42 invalid spectral method parameters
  !>
  !> @author IREDELL @date 96-04-10
  SUBROUTINE POLATEV4_grib2(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    INTEGER,          INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,          INTENT(IN   ) :: KM, MI, MO
    INTEGER,          INTENT(  OUT) :: IRET, IBO(KM), NO
    INTEGER,          INTENT(IN   ) :: IGDTNUMI, IGDTLENI
    INTEGER,          INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
    INTEGER,          INTENT(IN   ) :: IGDTNUMO, IGDTLENO
    INTEGER,          INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,             INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: CROT(MO),SROT(MO)
    !
    REAL,                 PARAMETER :: FILL=-9999.
    REAL,                 PARAMETER :: PI=3.14159265358979
    REAL,                 PARAMETER :: DPR=180./PI
    !
    INTEGER                         :: IDRTO, IROMB, ISKIPI, ISPEC
    INTEGER                         :: IDRTI, IMAXI, JMAXI, IM, JM
    INTEGER                         :: IPRIME, IG, IMO, JMO, IGO, JGO
    INTEGER                         :: ISCAN, JSCAN, NSCAN
    INTEGER                         :: ISCANO, JSCANO, NSCANO
    INTEGER                         :: ISCALE, IP, IPROJ, JSKIPI, JG
    INTEGER                         :: K, MAXWV, N, NI, NJ, NPS
    !
    REAL                            :: DLAT, DLON, DLATO, DLONO, DE, DR, DY
    REAL                            :: DUM, E2, H, HI, HJ
    REAL                            :: ORIENT, RERTH, SLAT
    REAL                            :: RLAT1, RLON1, RLAT2, RLON2, RLATI
    REAL                            :: UROT, VROT, UO2(MO,KM),VO2(MO,KM)
    REAL                            :: XMESH, X, XP, YP, XPTS(MO),YPTS(MO)

    type(grib2_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(igdtnumi, igdtleni, igdtmpli)
    desc_out = init_descriptor(igdtnumo, igdtleno, igdtmplo)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(IGDTNUMO.GE.0) THEN
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS, &
            RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  AFFIRM APPROPRIATE INPUT GRID
    !    LAT/LON OR GAUSSIAN
    !    NO BITMAPS
    !    FULL ZONAL COVERAGE
    !    FULL MERIDIONAL COVERAGE
    IDRTI=IGDTNUMI
    IF(IDRTI==40) IDRTI=4
    IF(IDRTI==0.OR.IDRTI==4)THEN
       IM=IGDTMPLI(8)
       JM=IGDTMPLI(9)
       ISCALE=IGDTMPLI(10)*IGDTMPLI(11)
       IF(ISCALE==0) ISCALE=10**6
       RLON1=FLOAT(IGDTMPLI(13))/FLOAT(ISCALE)
       RLON2=FLOAT(IGDTMPLI(16))/FLOAT(ISCALE)
       ISCAN=MOD(IGDTMPLI(19)/128,2)
       JSCAN=MOD(IGDTMPLI(19)/64,2)
       NSCAN=MOD(IGDTMPLI(19)/32,2)
    ELSE
       IRET=41
    ENDIF
    DO K=1,KM
       IF(IBI(K).NE.0) IRET=41
    ENDDO
    IF(IRET.EQ.0) THEN
       IF(ISCAN.EQ.0) THEN
          DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
       ELSE
          DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
       ENDIF
       IG=NINT(360/ABS(DLON))
       IPRIME=1+MOD(-NINT(RLON1/DLON)+IG,IG)
       IMAXI=IG
       JMAXI=JM
       IF(MOD(IG,2).NE.0.OR.IM.LT.IG) IRET=41
    ENDIF
    IF(IRET.EQ.0.AND.IDRTI.EQ.0) THEN
       ISCALE=IGDTMPLI(10)*IGDTMPLI(11)
       IF(ISCALE==0) ISCALE=10**6
       RLAT1=FLOAT(IGDTMPLI(12))/FLOAT(ISCALE)
       RLAT2=FLOAT(IGDTMPLI(15))/FLOAT(ISCALE)
       DLAT=(RLAT2-RLAT1)/(JM-1)
       JG=NINT(180/ABS(DLAT))
       IF(JM.EQ.JG) IDRTI=256
       IF(JM.NE.JG.AND.JM.NE.JG+1) IRET=41
    ELSEIF(IRET.EQ.0.AND.IDRTI.EQ.4) THEN
       JG=IGDTMPLI(18)*2
       IF(JM.NE.JG) IRET=41
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IF(IRET.EQ.0) THEN
       IROMB=IPOPT(1)
       MAXWV=IPOPT(2)
       IF(MAXWV.EQ.-1) THEN
          IF(IROMB.EQ.0.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)
          IF(IROMB.EQ.1.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.0.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/4
          IF(IROMB.EQ.0.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/4
       ENDIF
       IF((IROMB.NE.0.AND.IROMB.NE.1).OR.MAXWV.LT.0) IRET=42
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE
    IF(IRET.EQ.0) THEN
       IF(NSCAN.EQ.0) THEN
          ISKIPI=1
          JSKIPI=IM
       ELSE
          ISKIPI=JM
          JSKIPI=1
       ENDIF
       IF(ISCAN.EQ.1) ISKIPI=-ISKIPI
       IF(JSCAN.EQ.0) JSKIPI=-JSKIPI
       ISPEC=0
       !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
       IF((IGDTNUMO.EQ.0.OR.IGDTNUMO.EQ.40).AND. &
            MOD(IGDTMPLO(8),2).EQ.0.AND.IGDTMPLO(13).EQ.0.AND. &
            IGDTMPLO(19).EQ.0) THEN
          IDRTO=IGDTNUMO
          IF(IDRTO==40)IDRTO=4
          IMO=IGDTMPLO(8)
          JMO=IGDTMPLO(9)
          ISCALE=IGDTMPLO(10)*IGDTMPLO(11)
          IF(ISCALE==0) ISCALE=10**6
          RLON2=FLOAT(IGDTMPLO(16))/FLOAT(ISCALE)
          DLONO=(MOD(RLON2-1+3600,360.)+1)/(IMO-1)
          IGO=NINT(360/ABS(DLONO))
          IF(IMO.EQ.IGO.AND.IDRTO.EQ.0) THEN
             RLAT1=FLOAT(IGDTMPLO(12))/FLOAT(ISCALE)
             RLAT2=FLOAT(IGDTMPLO(15))/FLOAT(ISCALE)
             DLAT=(RLAT2-RLAT1)/(JMO-1)
             JGO=NINT(180/ABS(DLAT))
             IF(JMO.EQ.JGO) IDRTO=256
             IF(JMO.EQ.JGO.OR.JMO.EQ.JGO+1) ISPEC=1
          ELSEIF(IMO.EQ.IGO.AND.IDRTO.EQ.4) THEN
             JGO=IGDTMPLO(18)*2
             IF(JMO.EQ.JGO) ISPEC=1
          ENDIF
          IF(ISPEC.EQ.1) THEN
             CALL SPTRUNV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMO,JMO, &
                  KM,IPRIME,ISKIPI,JSKIPI,MI,0,0,MO,0,UI,VI, &
                  .TRUE.,UO,VO,.FALSE.,DUM,DUM,.FALSE.,DUM,DUM)
          ENDIF
          !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
       ELSEIF(IGDTNUMO.EQ.20.AND. &
            IGDTMPLO(8).EQ.IGDTMPLO(9).AND.MOD(IGDTMPLO(8),2).EQ.1.AND. &
            IGDTMPLO(15).EQ.IGDTMPLO(16).AND.IGDTMPLO(18).EQ.64.AND. &
            MOD(IGDTMPLO(12)/8,2).EQ.1) THEN
          NPS=IGDTMPLO(8)
          RLAT1=FLOAT(IGDTMPLO(10))*1.E-6
          RLON1=FLOAT(IGDTMPLO(11))*1.E-6
          ORIENT=FLOAT(IGDTMPLO(14))*1.E-6
          XMESH=FLOAT(IGDTMPLO(15))*1.E-3
          IPROJ=MOD(IGDTMPLO(17)/128,2)
          IP=(NPS+1)/2
          H=(-1.)**IPROJ
          SLAT=FLOAT(ABS(IGDTMPLO(13)))*1.E-6
          CALL EARTH_RADIUS(IGDTMPLO,IGDTLENO,RERTH,E2)
          DE=(1.+SIN(SLAT/DPR))*RERTH
          DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
          XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/XMESH
          YP=1+COS((RLON1-ORIENT)/DPR)*DR/XMESH
          IF(NINT(XP).EQ.IP.AND.NINT(YP).EQ.IP) THEN
             IF(IPROJ.EQ.0) THEN
                CALL SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     SLAT,XMESH,ORIENT,UI,VI,.TRUE.,UO,VO,UO2,VO2, &
                     .FALSE.,DUM,DUM,DUM,DUM, &
                     .FALSE.,DUM,DUM,DUM,DUM)
             ELSE
                CALL SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     SLAT,XMESH,ORIENT,UI,VI,.TRUE.,UO2,VO2,UO,VO, &
                     .FALSE.,DUM,DUM,DUM,DUM, &
                     .FALSE.,DUM,DUM,DUM,DUM)
             ENDIF
             ISPEC=1
          ENDIF
          !  SPECIAL CASE OF MERCATOR GRID
       ELSEIF(IGDTNUMO.EQ.10) THEN
          NI=IGDTMPLO(8)
          NJ=IGDTMPLO(9)
          RLAT1=FLOAT(IGDTMPLO(10))*1.0E-6
          RLON1=FLOAT(IGDTMPLO(11))*1.0E-6
          RLON2=FLOAT(IGDTMPLO(15))*1.0E-6
          RLATI=FLOAT(IGDTMPLO(13))*1.0E-6
          ISCANO=MOD(IGDTMPLO(16)/128,2)
          JSCANO=MOD(IGDTMPLO(16)/64,2)
          NSCANO=MOD(IGDTMPLO(16)/32,2)
          DY=FLOAT(IGDTMPLO(19))*1.0E-3
          HI=(-1.)**ISCANO
          HJ=(-1.)**(1-JSCANO)
          CALL EARTH_RADIUS(IGDTMPLO,IGDTLENO,RERTH,E2)
          DLONO=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(NI-1)
          DLATO=HJ*DY/(RERTH*COS(RLATI/DPR))*DPR
          IF(NSCANO.EQ.0) THEN
             CALL SPTRUNMV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NI,NJ, &
                  IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                  RLAT1,RLON1,DLATO,DLONO,UI,VI, &
                  .TRUE.,UO,VO,.FALSE.,DUM,DUM,.FALSE.,DUM,DUM)
             ISPEC=1
          ENDIF
       ENDIF
       !  GENERAL SLOW CASE
       IF(ISPEC.EQ.0) THEN
          CALL SPTRUNGV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
               IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON, &
               UI,VI,.TRUE.,UO,VO,.FALSE.,X,X,.FALSE.,X,X)
          DO K=1,KM
             IBO(K)=0
             DO N=1,NO
                LO(N,K)=.TRUE.
                UROT=CROT(N)*UO(N,K)-SROT(N)*VO(N,K)
                VROT=SROT(N)*UO(N,K)+CROT(N)*VO(N,K)
                UO(N,K)=UROT
                VO(N,K)=VROT
             ENDDO
          ENDDO
       ENDIF
    ELSE
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             UO(N,K)=0.
             VO(N,K)=0.
          ENDDO
       ENDDO
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLATEV4_GRIB2

  !> Interpolate vector fields (spectral).
  !>
  !> This subprogram performs spectral interpolation from any grid to
  !> any grid for vector fields. It requires that the input fields be
  !> uniformly global. Options allow choices between triangular shape
  !> (ipopt(1)=0) and rhomboidal shape (ipopt(1)=1) which has no
  !> default; a second option is the truncation (ipopt(2)) which
  !> defaults to a sensible truncation for the input grid (if
  !> opt(2)=-1).
  !>
  !> @note If the output grid is not found in a special list, then the
  !> transform back to grid is not very fast.  This special list
  !> contains global cylindrical grids, polar stereographic grids
  !> centered at the pole and mercator grids.
  !>
  !> Only horizontal interpolation is performed. The grids are defined
  !> by their grid description sections (passed in integer form as
  !> decoded by subprogram w3fi63).
  !>
  !> The current code recognizes the following projections:
  !> - (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
  !> - (KGDS(1)=001) MERCATOR CYLINDRICAL
  !> - (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
  !> - (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
  !> - (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
  !> - (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL (E-STAGGER)
  !> - (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL (B-STAGGER)
  !>
  !> Where kgds could be either input kgdsi or output kgdso.
  !>
  !> The input and output vectors are rotated so that they are either
  !> resolved relative to the defined grid in the direction of
  !> increasing x and y coordinates or resolved relative to easterly
  !> and northerly directions, as designated by their respective grid
  !> description sections. As an added bonus the number of output grid
  !> points and their latitudes and longitudes are also returned along
  !> with their vector rotation parameters. On the other hand, the
  !> output can be a set of station points if kgdso(1)<0, in which
  !> case the number of points and their latitudes and longitudes must
  !> be input along with their vector rotation parameters.
  !>
  !> Output bitmaps will only be created when the output grid extends
  !> outside of the domain of the input grid. The output field is set
  !> to 0 where the output bitmap is off.
  !>        
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | iredell | initial.
  !> 2001-06-18 | iredell |  improve detection of special fast transform
  !> 2015-01-27 | gayno | replace calls to gdswiz() with new merged routine gdswzd().
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1)=0 for
  !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
  !> number (defaults to sensible if ipopt(2)=-1).
  !> @param[in] kgdsi (200) input gds parameters as decoded by w3fi63.
  !> @param[in] kgdso (200) output gds parameters (kgdso(1)<0 implies
  !> random station points).
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1.
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags (must be all 0)
  !> @param[in] ui (mi,km) input u-component fields to interpolate
  !> @param[in] vi (mi,km) input v-component fields to interpolate
  !> @param[out] no number of output points (only if kgdso(1)<0)
  !> @param[out] rlat (no) output latitudes in degrees (if kgdso(1)<0)
  !> @param[out] rlon (no) output longitudes in degrees (if kgdso(1)<0)
  !> @param[out] crot (no) vector rotation cosines (if kgdso(1)<0)
  !> @param[out] srot (no) vector rotation sines (if kgdso(1)<0)
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] uo (mo,km) output u-component fields interpolated
  !> @param[out] vo (mo,km) output v-component fields interpolated
  !> @param[out] iret return code
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !> - 41 invalid nonglobal input grid
  !> - 42 invalid spectral method parameters
  !>
  !> @author IREDELL @date 96-04-10
  SUBROUTINE POLATEV4_grib1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    INTEGER,          INTENT(IN   ) :: IPOPT(20), IBI(KM)
    INTEGER,          INTENT(IN   ) :: KM, MI, MO
    INTEGER,          INTENT(  OUT) :: IRET, IBO(KM)
    INTEGER,          INTENT(IN) :: KGDSI(200),KGDSO(200)
    !
    LOGICAL*1,        INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,             INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,             INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    REAL,             INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,             INTENT(  OUT) :: CROT(MO),SROT(MO)
    !
    REAL,                 PARAMETER :: FILL=-9999.
    REAL,                 PARAMETER :: RERTH=6.3712E6
    REAL,                 PARAMETER :: PI=3.14159265358979
    REAL,                 PARAMETER :: DPR=180./PI
    !
    INTEGER                         :: IDRTO, IROMB, ISKIPI, ISPEC
    INTEGER                         :: IDRTI, IMAXI, JMAXI, IM, JM
    INTEGER                         :: IPRIME, IG, IMO, JMO, IGO, JGO
    INTEGER                         :: ISCAN, JSCAN, NSCAN
    INTEGER                         :: ISCANO, JSCANO, NSCANO
    INTEGER                         :: IP, IPROJ, JSKIPI, JG
    INTEGER                         :: K, MAXWV, N, NI, NJ, NO, NPS
    !
    REAL                            :: DLAT, DLON, DLATO, DLONO, DE, DR, DY
    REAL                            :: DUM, H, HI, HJ
    REAL                            :: ORIENT
    REAL                            :: RLAT1, RLON1, RLAT2, RLON2, RLATI
    REAL                            :: UROT, VROT, UO2(MO,KM),VO2(MO,KM)
    REAL                            :: XMESH, X, XP, YP, XPTS(MO),YPTS(MO)

    type(grib1_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(kgdsi)
    desc_out = init_descriptor(kgdso)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
    IRET=0
    IF(KGDSO(1).GE.0) THEN
       CALL GDSWZD(grid_out, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO,CROT,SROT)
       IF(NO.EQ.0) IRET=3
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  AFFIRM APPROPRIATE INPUT GRID
    !    LAT/LON OR GAUSSIAN
    !    NO BITMAPS
    !    FULL ZONAL COVERAGE
    !    FULL MERIDIONAL COVERAGE
    IDRTI=KGDSI(1)
    IM=KGDSI(2)
    JM=KGDSI(3)
    RLON1=KGDSI(5)*1.E-3
    RLON2=KGDSI(8)*1.E-3
    ISCAN=MOD(KGDSI(11)/128,2)
    JSCAN=MOD(KGDSI(11)/64,2)
    NSCAN=MOD(KGDSI(11)/32,2)
    IF(IDRTI.NE.0.AND.IDRTI.NE.4) IRET=41
    DO K=1,KM
       IF(IBI(K).NE.0) IRET=41
    ENDDO
    IF(IRET.EQ.0) THEN
       IF(ISCAN.EQ.0) THEN
          DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
       ELSE
          DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
       ENDIF
       IG=NINT(360/ABS(DLON))
       IPRIME=1+MOD(-NINT(RLON1/DLON)+IG,IG)
       IMAXI=IG
       JMAXI=JM
       IF(MOD(IG,2).NE.0.OR.IM.LT.IG) IRET=41
    ENDIF
    IF(IRET.EQ.0.AND.IDRTI.EQ.0) THEN
       RLAT1=KGDSI(4)*1.E-3
       RLAT2=KGDSI(7)*1.E-3
       DLAT=(RLAT2-RLAT1)/(JM-1)
       JG=NINT(180/ABS(DLAT))
       IF(JM.EQ.JG) IDRTI=256
       IF(JM.NE.JG.AND.JM.NE.JG+1) IRET=41
    ELSEIF(IRET.EQ.0.AND.IDRTI.EQ.4) THEN
       JG=KGDSI(10)*2
       IF(JM.NE.JG) IRET=41
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    IF(IRET.EQ.0) THEN
       IROMB=IPOPT(1)
       MAXWV=IPOPT(2)
       IF(MAXWV.EQ.-1) THEN
          IF(IROMB.EQ.0.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)
          IF(IROMB.EQ.1.AND.IDRTI.EQ.4) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.0.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.0) MAXWV=(JMAXI-3)/4
          IF(IROMB.EQ.0.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/2
          IF(IROMB.EQ.1.AND.IDRTI.EQ.256) MAXWV=(JMAXI-1)/4
       ENDIF
       IF((IROMB.NE.0.AND.IROMB.NE.1).OR.MAXWV.LT.0) IRET=42
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE
    IF(IRET.EQ.0) THEN
       IF(NSCAN.EQ.0) THEN
          ISKIPI=1
          JSKIPI=IM
       ELSE
          ISKIPI=JM
          JSKIPI=1
       ENDIF
       IF(ISCAN.EQ.1) ISKIPI=-ISKIPI
       IF(JSCAN.EQ.0) JSKIPI=-JSKIPI
       ISPEC=0
       !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
       IF((KGDSO(1).EQ.0.OR.KGDSO(1).EQ.4).AND. &
            MOD(KGDSO(2),2).EQ.0.AND.KGDSO(5).EQ.0.AND. &
            KGDSO(11).EQ.0) THEN
          IDRTO=KGDSO(1)
          IMO=KGDSO(2)
          JMO=KGDSO(3)
          RLON2=KGDSO(8)*1.E-3
          DLONO=(MOD(RLON2-1+3600,360.)+1)/(IMO-1)
          IGO=NINT(360/ABS(DLONO))
          IF(IMO.EQ.IGO.AND.IDRTO.EQ.0) THEN
             RLAT1=KGDSO(4)*1.E-3
             RLAT2=KGDSO(7)*1.E-3
             DLAT=(RLAT2-RLAT1)/(JMO-1)
             JGO=NINT(180/ABS(DLAT))
             IF(JMO.EQ.JGO) IDRTO=256
             IF(JMO.EQ.JGO.OR.JMO.EQ.JGO+1) ISPEC=1
          ELSEIF(IMO.EQ.IGO.AND.IDRTO.EQ.4) THEN
             JGO=KGDSO(10)*2
             IF(JMO.EQ.JGO) ISPEC=1
          ENDIF
          IF(ISPEC.EQ.1) THEN
             CALL SPTRUNV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMO,JMO, &
                  KM,IPRIME,ISKIPI,JSKIPI,MI,0,0,MO,0,UI,VI, &
                  .TRUE.,UO,VO,.FALSE.,DUM,DUM,.FALSE.,DUM,DUM)
          ENDIF
          !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
       ELSEIF(KGDSO(1).EQ.5.AND. &
            KGDSO(2).EQ.KGDSO(3).AND.MOD(KGDSO(2),2).EQ.1.AND. &
            KGDSO(8).EQ.KGDSO(9).AND.KGDSO(11).EQ.64.AND. &
            MOD(KGDSO(6)/8,2).EQ.1) THEN
          NPS=KGDSO(2)
          RLAT1=KGDSO(4)*1.E-3
          RLON1=KGDSO(5)*1.E-3
          ORIENT=KGDSO(7)*1.E-3
          XMESH=KGDSO(8)
          IPROJ=MOD(KGDSO(10)/128,2)
          IP=(NPS+1)/2
          H=(-1.)**IPROJ
          DE=(1.+SIN(60./DPR))*RERTH
          DR=DE*COS(RLAT1/DPR)/(1+H*SIN(RLAT1/DPR))
          XP=1-H*SIN((RLON1-ORIENT)/DPR)*DR/XMESH
          YP=1+COS((RLON1-ORIENT)/DPR)*DR/XMESH
          IF(NINT(XP).EQ.IP.AND.NINT(YP).EQ.IP) THEN
             IF(IPROJ.EQ.0) THEN
                CALL SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     60.,XMESH,ORIENT,UI,VI,.TRUE.,UO,VO,UO2,VO2, &
                     .FALSE.,DUM,DUM,DUM,DUM, &
                     .FALSE.,DUM,DUM,DUM,DUM)
             ELSE
                CALL SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NPS, &
                     IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                     60.,XMESH,ORIENT,UI,VI,.TRUE.,UO2,VO2,UO,VO, &
                     .FALSE.,DUM,DUM,DUM,DUM, &
                     .FALSE.,DUM,DUM,DUM,DUM)
             ENDIF
             ISPEC=1
          ENDIF
          !  SPECIAL CASE OF MERCATOR GRID
       ELSEIF(KGDSO(1).EQ.1) THEN
          NI=KGDSO(2)
          NJ=KGDSO(3)
          RLAT1=KGDSO(4)*1.E-3
          RLON1=KGDSO(5)*1.E-3
          RLON2=KGDSO(8)*1.E-3
          RLATI=KGDSO(9)*1.E-3
          ISCANO=MOD(KGDSO(11)/128,2)
          JSCANO=MOD(KGDSO(11)/64,2)
          NSCANO=MOD(KGDSO(11)/32,2)
          DY=KGDSO(13)
          HI=(-1.)**ISCANO
          HJ=(-1.)**(1-JSCANO)
          DLONO=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(NI-1)
          DLATO=HJ*DY/(RERTH*COS(RLATI/DPR))*DPR
          IF(NSCANO.EQ.0) THEN
             CALL SPTRUNMV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NI,NJ, &
                  IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0, &
                  RLAT1,RLON1,DLATO,DLONO,UI,VI, &
                  .TRUE.,UO,VO,.FALSE.,DUM,DUM,.FALSE.,DUM,DUM)
             ISPEC=1
          ENDIF
       ENDIF
       !  GENERAL SLOW CASE
       IF(ISPEC.EQ.0) THEN
          CALL SPTRUNGV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
               IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON, &
               UI,VI,.TRUE.,UO,VO,.FALSE.,X,X,.FALSE.,X,X)
          DO K=1,KM
             IBO(K)=0
             DO N=1,NO
                LO(N,K)=.TRUE.
                UROT=CROT(N)*UO(N,K)-SROT(N)*VO(N,K)
                VROT=SROT(N)*UO(N,K)+CROT(N)*VO(N,K)
                UO(N,K)=UROT
                VO(N,K)=VROT
             ENDDO
          ENDDO
       ENDIF
    ELSE
       DO K=1,KM
          IBO(K)=1
          DO N=1,NO
             LO(N,K)=.FALSE.
             UO(N,K)=0.
             VO(N,K)=0.
          ENDDO
       ENDDO
    ENDIF
  END SUBROUTINE POLATEV4_GRIB1
end module spectral_interp_mod
