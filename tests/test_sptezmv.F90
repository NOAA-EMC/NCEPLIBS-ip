! This is a test for NCEPLIBS-ip library.
!
! This tests the sptezmv() subroutine.
! 
! Alyson Stahl 8/2026
program test_sptezmv
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-5
  integer, parameter :: IROMB=0, MAXWV=7
  integer, parameter :: IMAX=16, JMAX=8, KMAX=8
  integer, parameter :: MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
  !
  integer :: IDRT=4
  integer :: IDIR
  ! SPTRANFV() inputs for expected values
  integer :: IP=1, IS=1, JB=1
  integer :: JN=IMAX, JS=-IMAX, KW=2*MX, KG=IMAX*JMAX, JE=JMAX/2
  integer :: JC
  !
  integer :: I, J, K, res
  real :: WAVED((MAXWV+1)*((IROMB+1)*MAXWV+2),KMAX)
  real :: WAVEZ((MAXWV+1)*((IROMB+1)*MAXWV+2),KMAX)
  real :: GRIDU(IMAX,JMAX,KMAX)
  real :: GRIDV(IMAX,JMAX,KMAX)
  ! 
  real :: EXP_WAVED((MAXWV+1)*((IROMB+1)*MAXWV+2),KMAX)
  real :: EXP_WAVEZ((MAXWV+1)*((IROMB+1)*MAXWV+2),KMAX)
  real :: EXP_GRIDU(IMAX,JMAX,KMAX)
  real :: EXP_GRIDV(IMAX,JMAX,KMAX)

  JC=NCPUS()

  ! Testing WAVE => GRID (IDIR > 0)
  IDIR = 1
  WAVED = 0.0
  WAVEZ = 0.0
  do K = 1, KMAX
    WAVED(1,K) = (2.0 + real(K)) * sqrt(2.0)
    WAVEZ(1,K) = (0.5 + 0.25 * real(K)) * sqrt(2.0)
  end do
  GRIDU = 0.0
  GRIDV = 0.0

  EXP_WAVED = WAVED
  EXP_WAVEZ = WAVEZ

  ! Caclulate expected values
  call SPTRANFV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
                  IP,IS,JN,JS,KW,KG,JB,JE,JC, &
                  EXP_WAVED,EXP_WAVEZ,        &
                  EXP_GRIDU,EXP_GRIDU(1,JMAX,1),EXP_GRIDV, &
                  EXP_GRIDV(1,JMAX,1),IDIR)

  call SPTEZMV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
               WAVED,WAVEZ,GRIDU,GRIDV,IDIR)
  

  ! WAVED & WAVEZ not expected to change, but we will check anyways
  res = 0
  do K = 1, KMAX
    do I = 1, MX
      if (abs(WAVED(I,K) - EXP_WAVED(I,K)) > TOL) then
        print *, 'Mismatch at WAVED(', I, ',', K, '): expected ', EXP_WAVED(I,K), ' but got ', WAVED(I,K)
        res = 1
      end if
      if (abs(WAVEZ(I,K) - EXP_WAVEZ(I,K)) > TOL) then
        print *, 'Mismatch at WAVEZ(', I, ',', K, '): expected ', EXP_WAVEZ(I,K), ' but got ', WAVEZ(I,K)
        res = 1
      end if
    end do
  end do

  do K = 1, KMAX
    do I = 1, IMAX
      do J = 1, JMAX
        if (abs(GRIDU(I,J,K) - EXP_GRIDU(I,J,K)) > TOL) then
          print *, 'Mismatch at GRIDU(', I, ',', J, ',', K, '): expected ', EXP_GRIDU(I,J,K), ' but got ', GRIDU(I,J,K)
          res = 1
        end if
        if (abs(GRIDV(I,J,K) - EXP_GRIDV(I,J,K)) > TOL) then
          print *, 'Mismatch at GRIDV(', I, ',', J, ',', K, '): expected ', EXP_GRIDV(I,J,K), ' but got ', GRIDV(I,J,K)
          res = 1
        end if
      end do
    end do
  end do

  if (res .ne. 0) stop 1

  ! Testing GRID => WAVE (IDIR < 0)
  IDIR = -1

  EXP_WAVED = 0.0
  EXP_WAVEZ = 0.0
  
  ! Caclulate expected values
  call SPTRANFV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
                  IP,IS,JN,JS,KW,KG,JB,JE,JC, &
                  EXP_WAVED,EXP_WAVEZ,        &
                  EXP_GRIDU,EXP_GRIDU(1,JMAX,1),EXP_GRIDV, &
                  EXP_GRIDV(1,JMAX,1),IDIR)

  call SPTEZMV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
               WAVED,WAVEZ,GRIDU,GRIDV,IDIR)

  res = 0
  do K = 1, KMAX
    do I = 1, MX
      if (abs(WAVED(I,K) - EXP_WAVED(I,K)) > TOL) then
        print *, 'Mismatch at WAVED(', I, ',', K, '): expected ', EXP_WAVED(I,K), ' but got ', WAVED(I,K)
        res = 1
      end if
      if (abs(WAVEZ(I,K) - EXP_WAVEZ(I,K)) > TOL) then
        print *, 'Mismatch at WAVEZ(', I, ',', K, '): expected ', EXP_WAVEZ(I,K), ' but got ', WAVEZ(I,K)
        res = 1
      end if
    end do
  end do

  ! GRIDU and GRIDV not expected to change, but we will check anyways
  do K = 1, KMAX
    do I = 1, IMAX
      do J = 1, JMAX
        if (abs(GRIDU(I,J,K) - EXP_GRIDU(I,J,K)) > TOL) then
          print *, 'Mismatch at GRIDU(', I, ',', J, ',', K, '): expected ', EXP_GRIDU(I,J,K), ' but got ', GRIDU(I,J,K)
          res = 1
        end if
        if (abs(GRIDV(I,J,K) - EXP_GRIDV(I,J,K)) > TOL) then
          print *, 'Mismatch at GRIDV(', I, ',', J, ',', K, '): expected ', EXP_GRIDV(I,J,K), ' but got ', GRIDV(I,J,K)
          res = 1
        end if
      end do
    end do
  end do

  if (res .ne. 0) stop 2

  print *, "Success!"
end program test_sptezmv