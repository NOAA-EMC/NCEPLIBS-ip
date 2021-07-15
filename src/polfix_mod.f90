module polfix_mod
  implicit none

  private
  public :: polfixs, polfixv

contains

  SUBROUTINE POLFIXS(NM,NX,KM,RLAT,IB,LO,GO)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLFIXS    MAKE MULTIPLE POLE SCALAR VALUES CONSISTENT
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM AVERAGES MULTIPLE POLE SCALAR VALUES
    !           ON A LATITUDE/LONGITUDE GRID.  BITMAPS MAY BE AVERAGED TOO.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !
    ! USAGE:    CALL POLFIXS(NM,NX,KM,RLAT,IB,LO,GO)
    !
    !   INPUT ARGUMENT LIST:
    !     NO       - INTEGER NUMBER OF GRID POINTS
    !     NX       - INTEGER LEADING DIMENSION OF FIELDS
    !     KM       - INTEGER NUMBER OF FIELDS
    !     RLAT     - REAL (NO) LATITUDES IN DEGREES
    !     IB       - INTEGER (KM) BITMAP FLAGS
    !     LO       - LOGICAL*1 (NX,KM) BITMAPS (IF SOME IB(K)=1)
    !     GO       - REAL (NX,KM) FIELDS
    !
    !   OUTPUT ARGUMENT LIST:
    !     LO       - LOGICAL*1 (NX,KM) BITMAPS (IF SOME IB(K)=1)
    !     GO       - REAL (NX,KM) FIELDS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER,    INTENT(IN   ) :: NM, NX, KM
    INTEGER,    INTENT(IN   ) :: IB(KM)
    !
    LOGICAL*1,  INTENT(INOUT) :: LO(NX,KM)
    !
    REAL,       INTENT(IN   ) :: RLAT(NM)
    REAL,       INTENT(INOUT) :: GO(NX,KM)
    !
    REAL,       PARAMETER     :: RLATNP=89.9995
    REAL,       PARAMETER     :: RLATSP=-RLATNP
    !
    INTEGER                   :: K, N
    !
    REAL                      :: WNP, GNP, TNP, WSP, GSP, TSP
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DO K=1,KM
       WNP=0.
       GNP=0.
       TNP=0.
       WSP=0.
       GSP=0.
       TSP=0.
       !  AVERAGE MULTIPLE POLE VALUES
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:WNP,GNP,TNP,WSP,GSP,TSP) SCHEDULE(STATIC)
       DO N=1,NM
          IF(RLAT(N).GE.RLATNP) THEN
             WNP=WNP+1
             IF(IB(K).EQ.0.OR.LO(N,K)) THEN
                GNP=GNP+GO(N,K)
                TNP=TNP+1
             ENDIF
          ELSEIF(RLAT(N).LE.RLATSP) THEN
             WSP=WSP+1
             IF(IB(K).EQ.0.OR.LO(N,K)) THEN
                GSP=GSP+GO(N,K)
                TSP=TSP+1
             ENDIF
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       !  DISTRIBUTE AVERAGE VALUES BACK TO MULTIPLE POLES
       IF(WNP.GT.1) THEN
          IF(TNP.GE.WNP/2) THEN
             GNP=GNP/TNP
          ELSE
             GNP=0.
          ENDIF
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NM
             IF(RLAT(N).GE.RLATNP) THEN
                IF(IB(K).NE.0) LO(N,K)=TNP.GE.WNP/2
                GO(N,K)=GNP
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
       IF(WSP.GT.1) THEN
          IF(TSP.GE.WSP/2) THEN
             GSP=GSP/TSP
          ELSE
             GSP=0.
          ENDIF
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NM
             IF(RLAT(N).LE.RLATSP) THEN
                IF(IB(K).NE.0) LO(N,K)=TSP.GE.WSP/2
                GO(N,K)=GSP
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLFIXS


  SUBROUTINE POLFIXV(NM,NX,KM,RLAT,RLON,IB,LO,UO,VO)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  POLFIXV    MAKE MULTIPLE POLE VECTOR VALUES CONSISTENT
    !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
    !
    ! ABSTRACT: THIS SUBPROGRAM AVERAGES MULTIPLE POLE VECTOR VALUES
    !           ON A LATITUDE/LONGITUDE GRID.  BITMAPS MAY BE AVERAGED TOO.
    !           VECTORS ARE ROTATED WITH RESPECT TO THEIR LONGITUDE.
    !        
    ! PROGRAM HISTORY LOG:
    !   96-04-10  IREDELL
    !
    ! USAGE:    CALL POLFIXV(NM,NX,KM,RLAT,RLON,IB,LO,UO,VO)
    !
    !   INPUT ARGUMENT LIST:
    !     NM       - INTEGER NUMBER OF GRID POINTS
    !     NX       - INTEGER LEADING DIMENSION OF FIELDS
    !     KM       - INTEGER NUMBER OF FIELDS
    !     RLAT     - REAL (NM) LATITUDES IN DEGREES
    !     RLON     - REAL (NM) LONGITUDES IN DEGREES
    !     IB       - INTEGER (KM) BITMAP FLAGS
    !     LO       - LOGICAL*1 (NX,KM) BITMAPS (IF SOME IB(K)=1)
    !     UO       - REAL (NX,KM) U-WINDS
    !     VO       - REAL (NX,KM) V-WINDS
    !
    !   OUTPUT ARGUMENT LIST:
    !     LO       - LOGICAL*1 (NX,KM) BITMAPS (IF SOME IB(K)=1)
    !     UO       - REAL (NX,KM) U-WINDS
    !     VO       - REAL (NX,KM) V-WINDS
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
    IMPLICIT NONE
    !
    INTEGER,      INTENT(IN   ) :: IB(KM), NM, NX, KM
    !
    LOGICAL*1,    INTENT(INOUT) :: LO(NX,KM)
    !
    REAL,         INTENT(IN   ) :: RLAT(NM), RLON(NM)
    REAL,         INTENT(INOUT) :: UO(NX,KM), VO(NX,KM)
    !
    REAL,         PARAMETER     :: RLATNP=89.9995
    REAL,         PARAMETER     :: RLATSP=-RLATNP
    REAL,         PARAMETER     :: PI=3.14159265358979
    REAL,         PARAMETER     :: DPR=180./PI
    !
    INTEGER                     :: K, N
    !
    REAL                        :: CLON(NM),SLON(NM)
    REAL                        :: TNP, UNP, VNP, WNP
    REAL                        :: TSP, USP, VSP, WSP
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
    DO N=1,NM
       CLON(N)=COS(RLON(N)/DPR)
       SLON(N)=SIN(RLON(N)/DPR)
    ENDDO
    !$OMP END PARALLEL DO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    DO K=1,KM
       WNP=0.
       UNP=0.
       VNP=0.
       TNP=0.
       WSP=0.
       USP=0.
       VSP=0.
       TSP=0.
       !  AVERAGE MULTIPLE POLE VALUES
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:WNP,UNP,VNP,TNP,WSP,USP,VSP,TSP) SCHEDULE(STATIC)
       DO N=1,NM
          IF(RLAT(N).GE.RLATNP) THEN
             WNP=WNP+1
             IF(IB(K).EQ.0.OR.LO(N,K)) THEN
                UNP=UNP+(CLON(N)*UO(N,K)-SLON(N)*VO(N,K))
                VNP=VNP+(SLON(N)*UO(N,K)+CLON(N)*VO(N,K))
                TNP=TNP+1
             ENDIF
          ELSEIF(RLAT(N).LE.RLATSP) THEN
             WSP=WSP+1
             IF(IB(K).EQ.0.OR.LO(N,K)) THEN
                USP=USP+(CLON(N)*UO(N,K)+SLON(N)*VO(N,K))
                VSP=VSP+(-SLON(N)*UO(N,K)+CLON(N)*VO(N,K))
                TSP=TSP+1
             ENDIF
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       !  DISTRIBUTE AVERAGE VALUES BACK TO MULTIPLE POLES
       IF(WNP.GT.1) THEN
          IF(TNP.GE.WNP/2) THEN
             UNP=UNP/TNP
             VNP=VNP/TNP
          ELSE
             UNP=0.
             VNP=0.
          ENDIF
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NM
             IF(RLAT(N).GE.RLATNP) THEN
                IF(IB(K).NE.0) LO(N,K)=TNP.GE.WNP/2
                UO(N,K)=CLON(N)*UNP+SLON(N)*VNP
                VO(N,K)=-SLON(N)*UNP+CLON(N)*VNP
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
       IF(WSP.GT.1) THEN
          IF(TSP.GE.WSP/2) THEN
             USP=USP/WSP
             VSP=VSP/WSP
          ELSE
             USP=0.
             VSP=0.
          ENDIF
          !$OMP PARALLEL DO PRIVATE(N) SCHEDULE(STATIC)
          DO N=1,NM
             IF(RLAT(N).LE.RLATSP) THEN
                IF(IB(K).NE.0) LO(N,K)=TSP.GE.WSP/2
                UO(N,K)=CLON(N)*USP-SLON(N)*VSP
                VO(N,K)=SLON(N)*USP+CLON(N)*VSP
             ENDIF
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE POLFIXV
end module polfix_mod
