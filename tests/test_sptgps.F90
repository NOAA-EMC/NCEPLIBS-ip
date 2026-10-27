! This is a test for NCEPLIBS-ip library.
!
! This tests the sptgps() subroutine.
! 
! Alyson Stahl 9/2026
program test_sptgps
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-4
  integer, parameter :: IROMB = 0, MAXWV = 1
  integer, parameter :: NPS = 1
  integer, parameter :: MX = (MAXWV+1)*((IROMB+1)*MAXWV+2)/2
  !
  integer :: i, res=0
  integer :: KWSKIP=0, KGSKIP=0, NISKIP=0, NJSKIP=0
  integer :: KMAX
  real :: TRUE, XMESH, ORIENT
  real :: WAVE(4*MX)
  real :: GN(2), GS(2), EXP_GN(2), EXP_GS(2)

  TRUE = 60.0
  XMESH = 381000.0
  ORIENT = 80.0
  WAVE = (/ 396.0, 0.0, 5.0, 0.0, 2.0, 1.0, &
            424.0, 0.0, 4.0, 0.0, 3.0, 1.0 /)

  print *, "Testing with KMAX = 1"
  KMAX = 1

  GN = 0.0
  GS = 0.0
  EXP_GN = (/ 286.13800970683081, 0.0 /)
  EXP_GS = (/ 273.89056099291491, 0.0 /)

  call SPTGPS(IROMB,MAXWV,KMAX,NPS,         &
              KWSKIP,KGSKIP,NISKIP,NJSKIP,  &
              TRUE,XMESH,ORIENT,WAVE,GN,GS)

  do i = 1, 2
    if (abs(GN(i) - EXP_GN(i)) > TOL) then
      print *, "Error: GN(", i, ") = ", GN(i), " Expected = ", EXP_GN(i)
      res = 1
    end if
    if (abs(GS(i) - EXP_GS(i)) > TOL) then
      print *, "Error: GS(", i, ") = ", GS(i), " Expected = ", EXP_GS(i)
      res = 1
    end if
  end do

  if (res .ne. 0) stop 1

  print *, "Testing with KMAX = 2"
  KMAX = 2

  GN = 0.0
  GS = 0.0
  EXP_GN = (/ 286.13800970683081, 304.71225470866256 /)
  EXP_GS = (/ 273.89056099291491, 294.91429573752981 /)

  call SPTGPS(IROMB,MAXWV,KMAX,NPS,         &
              KWSKIP,KGSKIP,NISKIP,NJSKIP,  &
              TRUE,XMESH,ORIENT,WAVE,GN,GS)

  do i = 1, 2
    if (abs(GN(i) - EXP_GN(i)) > TOL) then
      print *, "Error: GN(", i, ") = ", GN(i), " Expected = ", EXP_GN(i)
      res = 1
    end if
    if (abs(GS(i) - EXP_GS(i)) > TOL) then
      print *, "Error: GS(", i, ") = ", GS(i), " Expected = ", EXP_GS(i)
      res = 1
    end if
  end do

  if (res .ne. 0) stop 2

  print *, "Success!"
end program test_sptgps