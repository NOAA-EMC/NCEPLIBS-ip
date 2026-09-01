! This is a test for NCEPLIBS-ip library.
!
! This tests the sptezd() subroutine.
! 
! Alyson Stahl 8/2026
program test_sptezd
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-5
  integer, parameter :: IROMB=0, MAXWV=7
  integer, parameter :: IMAX=16, JMAX=8
  integer, parameter :: MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)
  !
  integer :: IDRT=4
  ! SPTRAND() inputs for expected values
  integer, parameter :: KMAX=1
  integer :: IPRIME=0, JBEG=0, JEND=0, IDIR=1
  integer :: ISKIP=0, JNSKIP=0, JSSKIP=0, KWSKIP=0, KGSKIP=0
  integer :: JCPU
  !
  integer :: I, J, K, res
  real :: WAVE(MX)
  real :: GRIDMN(KMAX), GRIDX(IMAX,JMAX), GRIDY(IMAX,JMAX)
  !
  real :: EXP_WAVE(MX)
  real :: EXP_GRIDMN(KMAX)
  real :: EXP_GRIDX(IMAX,JMAX)
  real :: EXP_GRIDY(IMAX,JMAX)

  JCPU = NCPUS()

  GRIDMN = 0.0
  GRIDX = 0.0
  GRIDY = 0.0
  WAVE = 0.0
  WAVE(1) = 12.0 * sqrt(2.0)

  ! Caclulate expected values
  EXP_WAVE = WAVE
  EXP_GRIDMN = 0.0
  EXP_GRIDX = 0.0
  EXP_GRIDY = 0.0

  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP, &
              JBEG,JEND,JCPU,EXP_WAVE,EXP_GRIDMN, &
              EXP_GRIDX,EXP_GRIDX(1,JMAX),&
              EXP_GRIDY,EXP_GRIDY(1,JMAX),IDIR)

  ! NOTE: SPTEZD() takes IDIR as input, but it does not use it. 
  ! Output for IDIR > 0 is automatically given.
  call SPTEZD(IROMB,MAXWV,IDRT,IMAX,JMAX, &
              WAVE,GRIDMN,GRIDX,GRIDY,IDIR)
              
  ! WAVE not expected to change, but we will check anyways
  res = 0
  do I = 1, KMAX*MX
    if (abs(WAVE(I) - EXP_WAVE(I)) > TOL) then
      print *, 'Mismatch at WAVE(', I, '): expected ', EXP_WAVE(I), ' but got ', WAVE(I)
      res = 1
    end if
  end do

  do K=1,KMAX
    if (abs(GRIDMN(K)-EXP_GRIDMN(K)) > TOL) then
      print *, 'Mismatch at GRIDMN(', K, '): expected ', EXP_GRIDMN(K), ' but got ', GRIDMN(K)
      res = 1
    end if
  end do

  do J=1,JMAX
    do I=1,IMAX
      if (abs(GRIDX(I,J)-EXP_GRIDX(I,J)) > TOL) then
        print *, 'Mismatch at GRIDX(', I, ',', J, '): expected ', EXP_GRIDX(I,J), ' but got ', GRIDX(I,J)
        res = 1
      end if
      if (abs(GRIDY(I,J)-EXP_GRIDY(I,J)) > TOL) then
        print *, 'Mismatch at GRIDY(', I, ',', J, '): expected ', EXP_GRIDY(I,J), ' but got ', GRIDY(I,J)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 1

  print *, "Success!"
end program test_sptezd