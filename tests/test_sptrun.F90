! This is a test for NCEPLIBS-ip library.
!
! This tests the sptrun() subroutine.
! 
! Alyson Stahl 9/2026
program test_sptrun
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-5
  integer, parameter :: IROMB=0, MAXWV=8
  integer, parameter :: IMAXI=26, JMAXI=13, KMAX=1
  integer, parameter :: IMAXO=26, JMAXO=27
  ! 
  integer :: IDRTI=4, IDRTO=0
  integer :: ISKIPI=0, JSKIPI=0, KSKIPI=IMAXI*JMAXI
  integer :: ISKIPO=0, JSKIPO=0, KSKIPO=IMAXO*JMAXO
  integer :: IPRIME=0, JCPU=0
  integer :: I, J, res
  real :: GRIDI(IMAXI, JMAXI, KMAX)
  real :: GRIDO(IMAXO, JMAXO, KMAX)
  real :: EXP_GRIDO(IMAXO, JMAXO, KMAX)

  GRIDO = 0.0
  GRIDI = 1.0
  EXP_GRIDO = 1.0
  
  call SPTRUN(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMAXO,JMAXO,  &
              KMAX,IPRIME,ISKIPI,JSKIPI,KSKIPI,                 &
              ISKIPO,JSKIPO,KSKIPO,JCPU,GRIDI,GRIDO)

  res = 0
  do J = 1, JMAXO
    do I = 1, IMAXO
      if (abs(GRIDO(I, J, 1) - EXP_GRIDO(I, J, 1)) > TOL) then
        print *, "GRIDO(", I, ",", J, ",1) = ", GRIDO(I, J, 1), " EXP_GRIDO = ", EXP_GRIDO(I, J, 1)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 1

  print *, "Success!"
end program test_sptrun