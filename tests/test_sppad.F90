! This is a test from the NCEPLIBS-sp project.
!
! This test tests the sppad() subrroutine.
!
! Alex Richert, Oct 2023
PROGRAM TEST_SPPAD
  IMPLICIT NONE

  INTEGER WHICH, I, IMAX1, IMAX2
  INTEGER, DIMENSION(6) :: I1=(/0,1,0,0,1,1/), I2=(/0,1,0,1,0,1/)
  INTEGER, DIMENSION(6) :: M1=(/384,384,1,1,1,1/), M2=(/384,384,2,2,2,2/)
  REAL, ALLOCATABLE :: Q1(:), Q2(:)
  REAL :: TINI=TINY(1.0)
  REAL, DIMENSION(18) :: W4REF, W6REF

  W4REF=(/1.0/6.0,1.0/3.0,0.5,2.0/3.0,0.0,0.0, &
    5.0/6.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
  W6REF=(/0.1250,0.25,0.375,0.5,0.0,0.0,0.625,0.75,0.875,1.0, &
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

  DO WHICH=1,6
    IMAX1=(M1(WHICH)+1)*((I1(WHICH)+1)*M1(WHICH)+2)
    IMAX2=(M2(WHICH)+1)*((I2(WHICH)+1)*M2(WHICH)+2)
    IF ((WHICH.LE.2) .AND. (IMAX1.NE.IMAX2)) STOP 1
    ALLOCATE(Q1(1:IMAX1))
    ALLOCATE(Q2(1:IMAX2))
    ! Make all test values positive to distinguish from padding zeros
    DO I=1,IMAX1
      Q1(I)=REAL(I)/REAL(IMAX1)
    ENDDO
    CALL SPPAD(I1(WHICH),M1(WHICH),Q1,I2(WHICH),M2(WHICH),Q2)
    ! When I1==I2, the arrays should be unchanged
    IF (WHICH.EQ.1) THEN
      IF (.NOT.ALL(ABS(Q1-Q2).LT.TINI)) STOP 2
    ENDIF
    IF (WHICH.EQ.2) THEN
      IF (.NOT.ALL(ABS(Q1-Q2).LT.TINI)) STOP 3
    ENDIF
    ! Non-pad values (i.e., non-zeros) should be unchanged
    IF (.NOT. ALL(ABS(Q1-PACK(Q2,Q2>TINI)).LT.TINI)) STOP 4
    IF (WHICH.EQ.4) THEN
      IF (.NOT.ALL(ABS(Q2-W4REF).LT.TINI)) STOP 5
    ENDIF
    IF (WHICH.EQ.6) THEN
      IF (.NOT.ALL(ABS(Q2-W6REF).LT.TINI)) STOP 6
    ENDIF
    DEALLOCATE(Q1)
    DEALLOCATE(Q2)
  ENDDO ! WHICH=1,6

END PROGRAM TEST_SPPAD
