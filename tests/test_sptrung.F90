! This is a test from the NCEPLIBS-sp project.
!
! This test tests the sptrung() subrroutine.
!
! Alex Richert, Oct 2023
program test_sptrung
  use sp_mod
  implicit none

  INTEGER                         :: I
  INTEGER                         :: IROMB=0, MAXWV=89
  INTEGER                         :: IDRTI=256, IMAXI=360, JMAXI=180
  INTEGER,parameter               :: KM=1
  INTEGER                         :: NO=4, IPRIME=181, ISKIPI=1, JSKIPI=360
  INTEGER,parameter               :: MI=64800, MO=4
  REAL                            :: RLAT(MO),RLON(MO)
  REAL                            :: GI(MI,KM)
  REAL                            :: GO(MO,KM)
  REAL                            :: GOREF(MO,KM)
  REAL*4                          :: RDGI(MI,KM)
  REAL                            :: TINI=4e-3


  OPEN (12, file="data/sptrung.gi.in", form='unformatted', recl=MI*KM*4, convert='little_endian')
  READ (12) RDGI
  CLOSE (12)

  GI = REAL(RDGI)

  GOREF(:,1) = (/77.3174667,70.4562683,72.8935242,51.0591698/)
  RLAT = (/45.0,35.0,40.0,35.0/)
  RLON = (/-100.0,-100.0,-90.0,-120.0/)

  CALL SPTRUNG(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KM,NO, &
       IPRIME,ISKIPI,JSKIPI,MI,MO,0,0,0,RLAT,RLON, &
       GI,GO)

  DO I=1,MO
    IF (ABS(GO(I,KM)-GOREF(I,KM)) .GT. TINI) STOP 1
  ENDDO

end program test_sptrung
