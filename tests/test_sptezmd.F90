! This is a test for NCEPLIBS-ip library.
!
! This tests the sptezmd() subroutine.
! 
! Alyson Stahl 8/2026
program test_sptezmd
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-5
  integer, parameter :: IROMB=0, MAXWV=7
  integer, parameter :: IMAX=16, JMAX=8, KMAX=2
  integer, parameter :: MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)
  !
  integer :: IDRT=4
  integer :: JCPU, IDIR
  !
  integer :: I, K, res
  real :: WAVE(MX,KMAX)
  real :: GRIDMN(KMAX), GRIDX(IMAX,JMAX,KMAX), GRIDY(IMAX,JMAX,KMAX)
  !
  real :: EXP_WAVE(MX,KMAX)

  JCPU = NCPUS()

  GRIDMN = 0.0
  GRIDX = 0.0
  GRIDY = 0.0

  WAVE = 0.0
  WAVE(1,1) = 12.0 * sqrt(2.0)
  WAVE(1,2) = 8.0 * sqrt(2.0)
  EXP_WAVE = WAVE
  
  IDIR = 1
  call SPTEZMD(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              WAVE,GRIDMN,GRIDX,GRIDY,IDIR)

  IDIR = -1
  call SPTEZMD(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              WAVE,GRIDMN,GRIDX,GRIDY,IDIR)
  
  res = 0
  do K = 1, KMAX
    do I = 1, MX
      if (abs(WAVE(I,K) - EXP_WAVE(I,K)) > TOL) then
        print *, 'Mismatch at WAVE(', I, ',', K, '): expected ', EXP_WAVE(I,K), ' but got ', WAVE(I,K)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 1

  print *, "Success!"
end program test_sptezmd