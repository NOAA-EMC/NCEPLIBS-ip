module spectral_interp_mod
  use gdswzd_mod
  use ip_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
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


  SUBROUTINE POLATES4_grib2(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATES4   INTERPOLATE SCALAR FIELDS (SPECTRAL)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS SPECTRAL INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
    !           IT REQUIRES THAT THE INPUT FIELDS BE UNIFORMLY GLOBAL.
    !           OPTIONS ALLOW CHOICES BETWEEN TRIANGULAR SHAPE (IPOPT(1)=0)
    !           AND RHOMBOIDAL SHAPE (IPOPT(1)=1) WHICH HAS NO DEFAULT;
    !           A SECOND OPTION IS THE TRUNCATION (IPOPT(2)) WHICH DEFAULTS 
    !           TO A SENSIBLE TRUNCATION FOR THE INPUT GRID (IF OPT(2)=-1).
    !           NOTE THAT IF THE OUTPUT GRID IS NOT FOUND IN A SPECIAL LIST,
    !           THEN THE TRANSFORM BACK TO GRID IS NOT VERY FAST.
    !           THIS SPECIAL LIST CONTAINS GLOBAL CYLINDRICAL GRIDS,
    !           POLAR STEREOGRAPHIC GRIDS CENTERED AT THE POLE
    !           AND MERCATOR GRIDS. ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
    !
    !           THE CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OnUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMO  =01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMO  =10) MERCATOR CYLINDRICAL
    !             (IGDTNUMO  =20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMO  =30) LAMBERT CONFORMAL CONICAL
    !             (IGDTNUMI/O=40) GAUSSIAN CYLINDRICAL
    !           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
    !           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
    !           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
    !           IF IGDTNUMO<0, IN WHICH CASE THE NUMBER OF POINTS
    !           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
    !           OUTPUT BITMAPS WILL NOT BE CREATED.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 2001-06-18  IREDELL  IMPROVE DETECTION OF SPECIAL FAST TRANSFORM
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      VERSION OF GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATES4(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,GI, &
    !                         NO,RLAT,RLON,IBO,LO,GO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1)=0 FOR TRIANGULAR, IPOPT(1)=1 FOR RHOMBOIDAL;
    !                IPOPT(2) IS TRUNCATION NUMBER
    !                (DEFAULTS TO SENSIBLE IF IPOPT(2)=-1).
    !     IGDTNUMI - INTEGER GRID DEFINITION TEMPLATE NUMBER - INPUT GRID.
    !                CORRESPONDS TO THE GFLD%IGDTNUM COMPONENT OF THE
    !                NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !                  00 - EQUIDISTANT CYLINDRICAL
    !                  01 - ROTATED EQUIDISTANT CYLINDRICAL.  "E"
    !                       AND NON-"E" STAGGERED
    !                  10 - MERCATOR CYCLINDRICAL
    !                  20 - POLAR STEREOGRAPHIC AZIMUTHAL
    !                  30 - LAMBERT CONFORMAL CONICAL
    !                  40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
    !     IGDTMPLI - INTEGER (IGDTLENI) GRID DEFINITION TEMPLATE ARRAY -
    !                INPUT GRID. CORRESPONDS TO THE GFLD%IGDTMPL COMPONENT
    !                OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE:
    !                (SECTION 3 INFO).  SEE COMMENTS IN ROUTINE
    !                IPOLATES FOR COMPLETE DEFINITION.
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
    !                IPOLATES FOR COMPLETE DEFINITION.
    !     IGDTLENO - INTEGER NUMBER OF ELEMENTS OF THE GRID DEFINITION
    !                TEMPLATE ARRAY - OUTPUT GRID.  CORRESPONDS TO THE GFLD%IGDTLEN
    !                COMPONENT OF THE NCEP G2 LIBRARY GRIDMOD DATA STRUCTURE.
    !     MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
    !     MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
    !                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
    !     KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS (MUST BE ALL 0)
    !     GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO<0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO<0)
    !
    !   OUTPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF IGDTNUMO>=0)
    !     RLAT     - REAL (MO) OUTPUT LATITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     RLON     - REAL (MO) OUTPUT LONGITUDES IN DEGREES (IF IGDTNUMO>=0)
    !     IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
    !     LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
    !     GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
    !     IRET     - INTEGER RETURN CODE
    !                0    SUCCESSFUL INTERPOLATION
    !                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
    !                3    UNRECOGNIZED OUTPUT GRID
    !                41   INVALID NONGLOBAL INPUT GRID
    !                42   INVALID SPECTRAL METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   EARTH_RADIUS DETERMINE SIZE/SHAPE OF EARTH
    !   GDSWZD       GRID DESCRIPTION SECTION WIZARD
    !   SPTRUN       SPECTRALLY TRUNCATE GRIDDED SCALAR FIELDS
    !   SPTRUNS      SPECTRALLY INTERPOLATE SCALARS TO POLAR STEREO.
    !   SPTRUNM      SPECTRALLY INTERPOLATE SCALARS TO MERCATOR
    !   SPTRUNG      SPECTRALLY INTERPOLATE SCALARS TO STATIONS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
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

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)
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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLATES4_GRIB2


  !> @file
  !! INTERPOLATE SCALAR FIELDS (SPECTRAL)
  !! @author IREDELL @date 96-04-10
  !
  !> THIS SUBPROGRAM PERFORMS SPECTRAL INTERPOLATION
  !!           FROM ANY GRID TO ANY GRID FOR SCALAR FIELDS.
  !!           IT REQUIRES THAT THE INPUT FIELDS BE UNIFORMLY GLOBAL.
  !!           OPTIONS ALLOW CHOICES BETWEEN TRIANGULAR SHAPE (IPOPT(1)=0)
  !!           AND RHOMBOIDAL SHAPE (IPOPT(1)=1) WHICH HAS NO DEFAULT;
  !!           A SECOND OPTION IS THE TRUNCATION (IPOPT(2)) WHICH DEFAULTS 
  !!           TO A SENSIBLE TRUNCATION FOR THE INPUT GRID (IF OPT(2)=-1).
  !!           NOTE THAT IF THE OUTPUT GRID IS NOT FOUND IN A SPECIAL LIST,
  !!           THEN THE TRANSFORM BACK TO GRID IS NOT VERY FAST.
  !!           THIS SPECIAL LIST CONTAINS GLOBAL CYLINDRICAL GRIDS,
  !!           POLAR STEREOGRAPHIC GRIDS CENTERED AT THE POLE
  !!           AND MERCATOR GRIDS.
  !!           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
  !!           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
  !!           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
  !!           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
  !!             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
  !!             (KGDS(1)=001) MERCATOR CYLINDRICAL
  !!             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
  !!             (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
  !!             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
  !!             (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL (E-STAGGER)
  !!             (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL (B-STAGGER)
  !!           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
  !!           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
  !!           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
  !!           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
  !!           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
  !!           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
  !!           OUTPUT BITMAPS WILL NOT BE CREATED.
  !!        
  !! PROGRAM HISTORY LOG:
  !! -  96-04-10  IREDELL
  !! - 2001-06-18  IREDELL  IMPROVE DETECTION OF SPECIAL FAST TRANSFORM
  !! - 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
  !!                      VERSION OF GDSWZD.
  !!
  !! @param IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
  !!                IPOPT(1)=0 FOR TRIANGULAR, IPOPT(1)=1 FOR RHOMBOIDAL;
  !!                IPOPT(2) IS TRUNCATION NUMBER
  !!                (DEFAULTS TO SENSIBLE IF IPOPT(2)=-1).
  !! @param KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
  !! @param KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
  !!                (KGDSO(1)<0 IMPLIES RANDOM STATION POINTS)
  !! @param MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
  !!                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
  !! @param MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
  !!                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
  !! @param KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
  !! @param IBI      - INTEGER (KM) INPUT BITMAP FLAGS (MUST BE ALL 0)
  !! @param LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
  !! @param GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
  !! @param[out] NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
  !! @param[out] RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
  !! @param[out] RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
  !! @param[out] IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
  !! @param[out] LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
  !! @param[out] GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
  !! @param[out] IRET     - INTEGER RETURN CODE
  !!                0    SUCCESSFUL INTERPOLATION
  !!                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
  !!                3    UNRECOGNIZED OUTPUT GRID
  !!                41   INVALID NONGLOBAL INPUT GRID
  !!                42   INVALID SPECTRAL METHOD PARAMETERS
  !!
  !! SUBPROGRAMS CALLED:
  !! -  GDSWZD       GRID DESCRIPTION SECTION WIZARD
  !! -  SPTRUN       SPECTRALLY TRUNCATE GRIDDED SCALAR FIELDS
  !! -  SPTRUNS      SPECTRALLY INTERPOLATE SCALARS TO POLAR STEREO.
  !! -  SPTRUNM      SPECTRALLY INTERPOLATE SCALARS TO MERCATOR
  !! -  SPTRUNG      SPECTRALLY INTERPOLATE SCALARS TO STATIONS
  !!
  SUBROUTINE POLATES4_grib1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,GI, &
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

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)
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


  SUBROUTINE POLATEV4_grib2(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLATEV4   INTERPOLATE VECTOR FIELDS (SPECTRAL)
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM PERFORMS SPECTRAL INTERPOLATION
    !           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
    !           IT REQUIRES THAT THE INPUT FIELDS BE UNIFORMLY GLOBAL.
    !           OPTIONS ALLOW CHOICES BETWEEN TRIANGULAR SHAPE (IPOPT(1)=0)
    !           AND RHOMBOIDAL SHAPE (IPOPT(1)=1) WHICH HAS NO DEFAULT;
    !           A SECOND OPTION IS THE TRUNCATION (IPOPT(2)) WHICH DEFAULTS 
    !           TO A SENSIBLE TRUNCATION FOR THE INPUT GRID (IF OPT(2)=-1).
    !           NOTE THAT IF THE OUTPUT GRID IS NOT FOUND IN A SPECIAL LIST,
    !           THEN THE TRANSFORM BACK TO GRID IS NOT VERY FAST.
    !           THIS SPECIAL LIST CONTAINS GLOBAL CYLINDRICAL GRIDS,
    !           POLAR STEREOGRAPHIC GRIDS CENTERED AT THE POLE
    !           AND MERCATOR GRIDS. ONLY HORIZONTAL INTERPOLATION
    !           IS PERFORMED.
    !
    !           THE INPUT AND OUTPUT GRIDS ARE DEFINED BY THEIR GRIB 2 GRID
    !           DEFINITION TEMPLATE AS DECODED BY THE NCEP G2 LIBRARY. THE
    !           CODE RECOGNIZES THE FOLLOWING PROJECTIONS, WHERE
    !           "IGDTNUMI/O" IS THE GRIB 2 GRID DEFINTION TEMPLATE NUMBER
    !           FOR THE INPUT AND OUTPUT GRIDS, RESPECTIVELY:
    !             (IGDTNUMI/O=00) EQUIDISTANT CYLINDRICAL
    !             (IGDTNUMO  =01) ROTATED EQUIDISTANT CYLINDRICAL. "E" AND
    !                             NON-"E" STAGGERED
    !             (IGDTNUMO  =10) MERCATOR CYLINDRICAL
    !             (IGDTNUMO  =20) POLAR STEREOGRAPHIC AZIMUTHAL
    !             (IGDTNUMO  =30) LAMBERT CONFORMAL CONICAL
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
    !           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
    !           IF IGDTNUMO<0, IN WHICH CASE THE NUMBER OF POINTS
    !           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT 
    !           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
    !
    !           OUTPUT BITMAPS WILL ONLY BE CREATED WHEN THE OUTPUT GRID
    !           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
    !           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    ! 2001-06-18  IREDELL  IMPROVE DETECTION OF SPECIAL FAST TRANSFORM
    ! 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
    !                      ROUTINE GDSWZD.
    ! 2015-07-13  GAYNO    CONVERT TO GRIB 2. REPLACE GRIB 1 KGDS ARRAYS
    !                      WITH GRIB 2 GRID DEFINITION TEMPLATE ARRAYS.
    !
    ! USAGE:    CALL POLATEV4(IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
    !                         IGDTNUMO,IGDTMPLO,IGDTLENO, &
    !                         MI,MO,KM,IBI,LI,UI,VI, &
    !                         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    !
    !   INPUT ARGUMENT LIST:
    !     IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
    !                IPOPT(1)=0 FOR TRIANGULAR, IPOPT(1)=1 FOR RHOMBOIDAL;
    !                IPOPT(2) IS TRUNCATION NUMBER
    !                (DEFAULTS TO SENSIBLE IF IPOPT(2)=-1).
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
    !     IBI      - INTEGER (KM) INPUT BITMAP FLAGS (MUST BE ALL 0)
    !     UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
    !     VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
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
    !                41   INVALID NONGLOBAL INPUT GRID
    !                42   INVALID SPECTRAL METHOD PARAMETERS
    !
    ! SUBPROGRAMS CALLED:
    !   EARTH_RADIUS DETERMINE SIZE/SHAPE OF EARTH
    !   GDSWZD       GRID DESCRIPTION SECTION WIZARD
    !   SPTRUNV      SPECTRALLY TRUNCATE GRIDDED VECTOR FIELDS
    !   SPTRUNSV     SPECTRALLY INTERPOLATE VECTORS TO POLAR STEREO.
    !   SPTRUNMV     SPECTRALLY INTERPOLATE VECTORS TO MERCATOR
    !   SPTRUNGV     SPECTRALLY INTERPOLATE VECTORS TO STATIONS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
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

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)

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

  !> @file
  !! INTERPOLATE VECTOR FIELDS (SPECTRAL)
  !! @author IREDELL @date 96-04-10
  !
  !> THIS SUBPROGRAM PERFORMS SPECTRAL INTERPOLATION
  !!           FROM ANY GRID TO ANY GRID FOR VECTOR FIELDS.
  !!           IT REQUIRES THAT THE INPUT FIELDS BE UNIFORMLY GLOBAL.
  !!           OPTIONS ALLOW CHOICES BETWEEN TRIANGULAR SHAPE (IPOPT(1)=0)
  !!           AND RHOMBOIDAL SHAPE (IPOPT(1)=1) WHICH HAS NO DEFAULT;
  !!           A SECOND OPTION IS THE TRUNCATION (IPOPT(2)) WHICH DEFAULTS 
  !!           TO A SENSIBLE TRUNCATION FOR THE INPUT GRID (IF OPT(2)=-1).
  !!           NOTE THAT IF THE OUTPUT GRID IS NOT FOUND IN A SPECIAL LIST,
  !!           THEN THE TRANSFORM BACK TO GRID IS NOT VERY FAST.
  !!           THIS SPECIAL LIST CONTAINS GLOBAL CYLINDRICAL GRIDS,
  !!           POLAR STEREOGRAPHIC GRIDS CENTERED AT THE POLE
  !!           AND MERCATOR GRIDS.
  !!           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
  !!           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
  !!           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
  !!           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
  !!             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
  !!             (KGDS(1)=001) MERCATOR CYLINDRICAL
  !!             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
  !!             (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
  !!             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
  !!             (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL (E-STAGGER)
  !!             (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL (B-STAGGER)
  !!           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
  !!           THE INPUT AND OUTPUT VECTORS ARE ROTATED SO THAT THEY ARE
  !!           EITHER RESOLVED RELATIVE TO THE DEFINED GRID
  !!           IN THE DIRECTION OF INCREASING X AND Y COORDINATES
  !!           OR RESOLVED RELATIVE TO EASTERLY AND NORTHERLY DIRECTIONS,
  !!           AS DESIGNATED BY THEIR RESPECTIVE GRID DESCRIPTION SECTIONS.
  !!           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
  !!           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED
  !!           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  !!           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
  !!           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
  !!           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT 
  !!           ALONG WITH THEIR VECTOR ROTATION PARAMETERS.
  !!           OUTPUT BITMAPS WILL ONLY BE CREATED WHEN THE OUTPUT GRID
  !!           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
  !!           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
  !!        
  !! PROGRAM HISTORY LOG:
  !! -  96-04-10  IREDELL
  !! - 2001-06-18  IREDELL  IMPROVE DETECTION OF SPECIAL FAST TRANSFORM
  !! - 2015-01-27  GAYNO    REPLACE CALLS TO GDSWIZ WITH NEW MERGED
  !!                      ROUTINE GDSWZD.
  !!
  !! @param IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
  !!                IPOPT(1)=0 FOR TRIANGULAR, IPOPT(1)=1 FOR RHOMBOIDAL;
  !!                IPOPT(2) IS TRUNCATION NUMBER
  !!                (DEFAULTS TO SENSIBLE IF IPOPT(2)=-1).
  !! @param KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
  !! @param KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
  !!                (KGDSO(1)<0 IMPLIES RANDOM STATION POINTS)
  !! @param MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
  !!                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
  !! @param MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
  !!                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
  !! @param KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
  !! @param IBI      - INTEGER (KM) INPUT BITMAP FLAGS (MUST BE ALL 0)
  !! @param LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF SOME IBI(K)=1)
  !! @param UI       - REAL (MI,KM) INPUT U-COMPONENT FIELDS TO INTERPOLATE
  !! @param VI       - REAL (MI,KM) INPUT V-COMPONENT FIELDS TO INTERPOLATE
  !! @param[out] NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
  !! @param[out] RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
  !! @param[out] RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
  !! @param[out] CROT     - REAL (NO) VECTOR ROTATION COSINES (IF KGDSO(1)<0)
  !! @param[out] SROT     - REAL (NO) VECTOR ROTATION SINES (IF KGDSO(1)<0)
  !!                (UGRID=CROT*UEARTH-SROT*VEARTH;
  !!                 VGRID=SROT*UEARTH+CROT*VEARTH)
  !! @param[out] IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
  !! @param[out] LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
  !! @param[out] UO       - REAL (MO,KM) OUTPUT U-COMPONENT FIELDS INTERPOLATED
  !! @param[out] VO       - REAL (MO,KM) OUTPUT V-COMPONENT FIELDS INTERPOLATED
  !! @param[out] IRET     - INTEGER RETURN CODE
  !!                0    SUCCESSFUL INTERPOLATION
  !!                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
  !!                3    UNRECOGNIZED OUTPUT GRID
  !!                41   INVALID NONGLOBAL INPUT GRID
  !!                42   INVALID SPECTRAL METHOD PARAMETERS
  !!
  !! SUBPROGRAMS CALLED:
  !! -  GDSWZD       GRID DESCRIPTION SECTION WIZARD
  !! -  SPTRUNV      SPECTRALLY TRUNCATE GRIDDED VECTOR FIELDS
  !! -  SPTRUNSV     SPECTRALLY INTERPOLATE VECTORS TO POLAR STEREO.
  !! -  SPTRUNMV     SPECTRALLY INTERPOLATE VECTORS TO MERCATOR
  !! -  SPTRUNGV     SPECTRALLY INTERPOLATE VECTORS TO STATIONS
  !!
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

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)
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
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLATEV4_GRIB1



end module spectral_interp_mod
