! This is a test for NCEPLIBS-ip library.
!
! This tests the sptrand() subroutine.
! 
! Alyson Stahl 8/2026
program test_sptrand
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-5
  integer, parameter :: IROMB=0, MAXWV=7
  integer, parameter :: IMAX=16, JMAX=8, KMAX=4
  integer, parameter :: MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)
  !
  integer :: IDRT=4, IPRIME=0, JEND=0
  integer :: ISKIP=0, JNSKIP=0, JSSKIP=0, KGSKIP=0
  integer :: KWSKIP, JBEG, IDIR, JCPU
  !
  real :: WAVE(MX,KMAX)
  real :: GRIDXN(IMAX,JMAX,KMAX), GRIDXS(IMAX,JMAX,KMAX)
  real :: GRIDYN(IMAX,JMAX,KMAX), GRIDYS(IMAX,JMAX,KMAX)
  !
  real :: EXP_WAVE(MX,KMAX)
  !
  integer :: I, J, res

  JCPU = NCPUS()

  EXP_WAVE = 0.0
  EXP_WAVE(1,1) = 12.0 * sqrt(2.0)

  ! Test Case 1: WAVE => GRID => WAVE with KWSKIP = 0 & JBEG = 0. Expect WAVE to remain unchanged.
  KWSKIP = 0
  JBEG = 0
  IDIR = 1
  WAVE = 0.0
  WAVE(1,1) = 12.0 * sqrt(2.0)
  GRIDMN = 0.0
  GRIDXN = 0.0
  GRIDXS = 0.0
  GRIDYN = 0.0
  GRIDYS = 0.0

  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  IDIR = -1
  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  res = 0
  do I = 1, MX
    do J = 1, KMAX
      if (abs(WAVE(I,J) - EXP_WAVE(I,J)) > TOL) then
        print *, 'Mismatch at WAVE(', I, ',', J, '): expected ', EXP_WAVE(I,J), ' but got ', WAVE(I,J)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 1

  ! Test Case 2: WAVE => GRID => WAVE with KWSKIP > 0 & JBEG = 0. Expect WAVE to remain unchanged.
  KWSKIP = MX
  JBEG = 0
  IDIR = 1
  WAVE = 0.0
  WAVE(1,1) = 12.0 * sqrt(2.0)
  GRIDMN = 0.0
  GRIDXN = 0.0
  GRIDXS = 0.0
  GRIDYN = 0.0
  GRIDYS = 0.0

  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  IDIR = -1
  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  res = 0
  do I = 1, MX
    do J = 1, KMAX
      if (abs(WAVE(I,J) - EXP_WAVE(I,J)) > TOL) then
        print *, 'Mismatch at WAVE(', I, ',', J, '): expected ', EXP_WAVE(I,J), ' but got ', WAVE(I,J)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 2

  ! Test Case 3: WAVE => GRID => WAVE with KWSKIP = 0 & JBEG > 0.
  KWSKIP = 0
  JBEG = 1
  IDIR = 1
  WAVE = 0.0
  WAVE(1,1) = 12.0 * sqrt(2.0)
  GRIDMN = 0.0
  GRIDXN = 0.0
  GRIDXS = 0.0
  GRIDYN = 0.0
  GRIDYS = 0.0

  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  IDIR = -1
  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  res = 0
  do I = 1, MX
    do J = 1, KMAX
      if (abs(WAVE(I,J) - EXP_WAVE(I,J)) > TOL) then
        print *, 'Mismatch at WAVE(', I, ',', J, '): expected ', EXP_WAVE(I,J), ' but got ', WAVE(I,J)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 3

  ! Test Case 4: WAVE => GRID => WAVE with KWSKIP > 0 & JBEG > 0.
  KWSKIP = MX
  JBEG = 1
  IDIR = 1
  WAVE = 0.0
  WAVE(1,1) = 12.0 * sqrt(2.0)
  GRIDMN = 0.0
  GRIDXN = 0.0
  GRIDXS = 0.0
  GRIDYN = 0.0
  GRIDYS = 0.0

  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  IDIR = -1
  call SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX, &
              IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,&
              JBEG,JEND,JCPU,&
              WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

  res = 0
  do I = 1, MX
    do J = 1, KMAX
      if (abs(WAVE(I,J) - EXP_WAVE(I,J)) > TOL) then
        print *, 'Mismatch at WAVE(', I, ',', J, '): expected ', EXP_WAVE(I,J), ' but got ', WAVE(I,J)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 4

  print *, "Success!"
end program test_sptrand