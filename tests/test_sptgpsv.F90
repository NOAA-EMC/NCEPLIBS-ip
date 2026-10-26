! This is a test for NCEPLIBS-ip library.
!
! This tests the sptgpsv() subroutine.
! 
! Alyson Stahl 9/2026
program test_sptgpsv
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
  real :: WAVED(4*MX), WAVEZ(4*MX)
  real :: UN(2), VN(2), US(2), VS(2)
  real :: EXP_UN(2), EXP_VN(2), EXP_US(2), EXP_VS(2)

  TRUE = 60.0
  XMESH = 381000.0
  ORIENT = 80.0
  WAVED = (/ 0.0, 0.0, 1.0e-6, 0.0, 2.0e-6, 1.0e-6, &
             0.0, 0.0, 2.0e-6, 0.0, 1.0e-6, 2.0e-6 /)
  WAVEZ = (/ 0.0, 0.0, 1.0e-5, 0.0, 3.0e-5, 2.0e-5, &
             0.0, 0.0, 2.0e-5, 0.0, 2.0e-5, 4.0e-5 /)

  print *, "Testing with KMAX = 1"
  KMAX = 1

  UN = 0.0
  VN = 0.0
  US = 0.0
  VS = 0.0
  EXP_UN = (/ 91.757891412722188, 0.0 /)
  EXP_VN = (/ 178.65883023094909, 0.0 /)
  EXP_US = (/ -68.106457767942189, 0.0 /)
  EXP_VS = (/ 185.69392284914306, 0.0 /)

  call SPTGPSV(IROMB,MAXWV,KMAX,NPS,        &
              KWSKIP,KGSKIP,NISKIP,NJSKIP,  &
              TRUE,XMESH,ORIENT,WAVED,WAVEZ,UN,VN,US,VS)

  do i = 1, 2
    if (abs(UN(i) - EXP_UN(i)) > TOL) then
      print *, "Error: UN(", i, ") = ", UN(i), " Expected = ", EXP_UN(i)
      res = 1
    end if
    if (abs(VN(i) - EXP_VN(i)) > TOL) then
      print *, "Error: VN(", i, ") = ", VN(i), " Expected = ", EXP_VN(i)
      res = 1
    end if
    if (abs(US(i) - EXP_US(i)) > TOL) then
      print *, "Error: US(", i, ") = ", US(i), " Expected = ", EXP_US(i)
      res = 1
    end if
    if (abs(VS(i) - EXP_VS(i)) > TOL) then
      print *, "Error: VS(", i, ") = ", VS(i), " Expected = ", EXP_VS(i)
      res = 1
    end if
  end do

  if (res .ne. 0) stop 1

  print *, "Testing with KMAX = 2"
  KMAX = 2

  UN = 0.0
  VN = 0.0
  US = 0.0
  VS = 0.0
  EXP_UN = (/ 91.757891412722188, 205.53938848668605 /)
  EXP_VN = (/ 178.65883023094909, 137.09144630836437 /)
  EXP_US = (/ -68.106457767942189, -190.83929714177845 /)
  EXP_VS = (/ 185.69392284914306, 156.91038058978759 /)

  call SPTGPSV(IROMB,MAXWV,KMAX,NPS,        &
              KWSKIP,KGSKIP,NISKIP,NJSKIP,  &
              TRUE,XMESH,ORIENT,WAVED,WAVEZ,UN,VN,US,VS)

  do i = 1, 2
    if (abs(UN(i) - EXP_UN(i)) > TOL) then
      print *, "Error: UN(", i, ") = ", UN(i), " Expected = ", EXP_UN(i)
      res = 1
    end if
    if (abs(VN(i) - EXP_VN(i)) > TOL) then
      print *, "Error: VN(", i, ") = ", VN(i), " Expected = ", EXP_VN(i)
      res = 1
    end if
    if (abs(US(i) - EXP_US(i)) > TOL) then
      print *, "Error: US(", i, ") = ", US(i), " Expected = ", EXP_US(i)
      res = 1
    end if
    if (abs(VS(i) - EXP_VS(i)) > TOL) then
      print *, "Error: VS(", i, ") = ", VS(i), " Expected = ", EXP_VS(i)
      res = 1
    end if
  end do

  if (res .ne. 0) stop 2

  print *, "Success!"
end program test_sptgpsv