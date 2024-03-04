! This is a test from the NCEPLIBS-sp project.
!
! This test tests the sptrungv() subrroutine.
!
! Alex Richert, Oct 2023
program test_sptrungv
  use sp_mod
  implicit none

  !
  !
  INTEGER                         :: I
  INTEGER                         :: IROMB=0, MAXWV=89
  INTEGER                         :: IDRTI=0, IMAXI=360, JMAXI=181
  INTEGER,parameter               :: KM=1
  INTEGER                         :: NO=26553, IPRIME=1, ISKIPI=1, JSKIPI=-360
  INTEGER,parameter               :: MI=65160, MO=26553
  REAL                            :: RLAT(MO),RLON(MO)
  REAL                            :: UI(MI,KM),VI(MI,KM)
  REAL                            :: UO(MO,KM),VO(MO,KM)
  REAL                            :: UOREF(MO,KM),VOREF(MO,KM)
  REAL*4                          :: RDRLAT(MO),RDRLON(MO)
  REAL*4                          :: RDUI(MI,KM),RDVI(MI,KM)
  REAL*4                          :: RDUOREF(MO,KM),RDVOREF(MO,KM)
  REAL                            :: X(1)=0.0
  REAL                            :: TOL=1e-2

  OPEN (12, file="data/sptrungv.uv.in", access='direct', recl=MI*KM*4, convert='little_endian')
  READ (12, rec=1) RDUI
  READ (12, rec=2) RDVI
  CLOSE (12)

  OPEN (13, file="data/sptrungv.ll.in", access="direct", recl=MO*4, convert='little_endian')
  READ (13, rec=1) RDRLAT
  READ (13, rec=2) RDRLON
  CLOSE (13)

  OPEN (14, file="data/sptrungv.uv.out", access="direct", recl=MO*KM*4, convert='little_endian')
  READ (14, rec=1) RDUOREF
  READ (14, rec=2) RDVOREF
  CLOSE (14)

  UI = REAL(RDUI)
  VI = REAL(RDVI)
  RLAT = REAL(RDRLAT)
  RLON = REAL(RDRLON)
  UOREF = REAL(RDUOREF)
  VOREF = REAL(RDVOREF)

  CALL SPTRUNGV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
       IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON, &
       UI,VI,.TRUE.,UO,VO,.FALSE.,X,X,.FALSE.,X,X)

  DO I=1, MO
    IF (ABS(UO(I,KM)-RDUOREF(I,KM)) .GT. TOL) STOP 1
    IF (ABS(VO(I,KM)-RDVOREF(I,KM)) .GT. TOL) STOP 2
  ENDDO

end program test_sptrungv
