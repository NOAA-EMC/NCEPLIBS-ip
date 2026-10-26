! This is a test for NCEPLIBS-ip library.
!
! This tests the sptrunv() subroutine.
! 
! Alyson Stahl 9/2026
program test_sptrunv
  use sp_mod
  implicit none

  real, parameter :: TOL = 1e-6
  integer, parameter :: IROMB=0, MAXWV=1
  integer, parameter :: IMAXI=4, JMAXI=7, KMAX=1
  integer, parameter :: IMAXO=4, JMAXO=3
  ! 
  integer :: IDRTI=0, IDRTO=4
  integer :: ISKIPI=0, JSKIPI=0, KSKIPI=IMAXI*JMAXI
  integer :: ISKIPO=0, JSKIPO=0, KSKIPO=IMAXO*JMAXO
  integer :: IPRIME=0, JCPU=0
  logical :: LUV, LDZ, LPS
  integer :: I, J, res=0
  real :: GRIDUI(IMAXI, JMAXI, KMAX), GRIDVI(IMAXI, JMAXI, KMAX)
  ! Output grids
  real :: GRIDUO(IMAXO, JMAXO, KMAX), GRIDVO(IMAXO, JMAXO, KMAX)
  real :: GRIDDO(IMAXO, JMAXO, KMAX), GRIDZO(IMAXO, JMAXO, KMAX)
  real :: GRIDPO(IMAXO, JMAXO, KMAX), GRIDSO(IMAXO, JMAXO, KMAX)
  !
  real :: EXP_GRIDUO(IMAXO, JMAXO, KMAX), EXP_GRIDVO(IMAXO, JMAXO, KMAX)
  real :: EXP_GRIDDO(IMAXO, JMAXO, KMAX), EXP_GRIDZO(IMAXO, JMAXO, KMAX)
  real :: EXP_GRIDPO(IMAXO, JMAXO, KMAX), EXP_GRIDSO(IMAXO, JMAXO, KMAX)

  GRIDUI = 0.0
  GRIDVI = 0.0
  
  GRIDUI(:, 1, 1) = 0.0
  GRIDUI(:, 2, 1) = 30.0
  GRIDUI(:, 3, 1) = 51.0
  GRIDUI(:, 4, 1) = 60.0
  GRIDUI(:, 5, 1) = 51.0
  GRIDUI(:, 6, 1) = 30.0
  GRIDUI(:, 7, 1) = 0.0

  GRIDVI(:, 1, 1) = 0.0
  GRIDVI(:, 2, 1) = 30.0
  GRIDVI(:, 3, 1) = 51.0
  GRIDVI(:, 4, 1) = 60.0
  GRIDVI(:, 5, 1) = 51.0
  GRIDVI(:, 6, 1) = 30.0
  GRIDVI(:, 7, 1) = 0.0

  ! Test Case 1: LUV, LDX, and LPS = FALSE
  GRIDUO = 0.0
  GRIDVO = 0.0
  GRIDDO = 0.0
  GRIDZO = 0.0
  GRIDPO = 0.0
  GRIDSO = 0.0
  
  LUV = .FALSE.
  LDZ = .FALSE.
  LPS = .FALSE.

  print *, "Testing SPTRUNV with LUV, LDZ, and LPS = FALSE"
  call SPTRUNV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,   &
              IDRTO,IMAXO,JMAXO,KMAX, &
              IPRIME,ISKIPI,JSKIPI,KSKIPI,  &
              ISKIPO,JSKIPO,KSKIPO,JCPU,GRIDUI,GRIDVI,  &
              LUV,GRIDUO,GRIDVO,LDZ,GRIDDO,GRIDZO,  &
              LPS,GRIDPO,GRIDSO)

  do I = 1, IMAXO
    do J = 1, JMAXO
      if (abs(GRIDUO(I, J, 1)) > TOL) then
        print *, "GRIDUO(", I, ",", J, ",1) = ", GRIDUO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDVO(I, J, 1)) > TOL) then
        print *, "GRIDVO(", I, ",", J, ",1) = ", GRIDVO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDDO(I, J, 1)) > TOL) then
        print *, "GRIDDO(", I, ",", J, ",1) = ", GRIDDO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDZO(I, J, 1)) > TOL) then
        print *, "GRIDZO(", I, ",", J, ",1) = ", GRIDZO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDPO(I, J, 1)) > TOL) then
        print *, "GRIDPO(", I, ",", J, ",1) = ", GRIDPO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDSO(I, J, 1)) > TOL) then
        print *, "GRIDSO(", I, ",", J, ",1) = ", GRIDSO(I, J, 1)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 1

  ! Test Case 2: LUV, LDX, and LPS = TRUE
  GRIDUO = 0.0
  GRIDVO = 0.0
  GRIDDO = 0.0
  GRIDZO = 0.0
  GRIDPO = 0.0
  GRIDSO = 0.0

  LUV = .TRUE.
  LDZ = .TRUE.
  LPS = .TRUE.

#if (LSIZE==4)
  EXP_GRIDUO(:,1,1) = 37.6313362
  EXP_GRIDUO(:,2,1) = 59.5003738
  EXP_GRIDUO(:,3,1) = 37.6313362 

  EXP_GRIDVO(:,1,1) = 37.6313362
  EXP_GRIDVO(:,2,1) = 59.5003738
  EXP_GRIDVO(:,3,1) = 37.6313362 

  EXP_GRIDDO(:,1,1) = -1.44678534E-05
  EXP_GRIDDO(:,2,1) = 0.0
  EXP_GRIDDO(:,3,1) = 1.44678534E-05

  EXP_GRIDZO(:,1,1) = 1.44678534E-05
  EXP_GRIDZO(:,2,1) = 0.0
  EXP_GRIDZO(:,3,1) = -1.44678534E-05
  
  EXP_GRIDPO(:,1,1) = 293640928. 
  EXP_GRIDPO(:,2,1) = 0.0
  EXP_GRIDPO(:,3,1) = -293640928. 

  EXP_GRIDSO(:,1,1) = -293640928. 
  EXP_GRIDSO(:,2,1) = 0.0
  EXP_GRIDSO(:,3,1) = 293640928. 
# else
  EXP_GRIDUO(:,1,1) = 37.631342816769383
  EXP_GRIDUO(:,2,1) = 59.500377355803828
  EXP_GRIDUO(:,3,1) = 37.631342816769383

  EXP_GRIDVO(:,1,1) = 37.631342816769383
  EXP_GRIDVO(:,2,1) = 59.500377355803828
  EXP_GRIDVO(:,3,1) = 37.631342816769383

  EXP_GRIDDO(:,1,1) = -1.4467853502767779E-005
  EXP_GRIDDO(:,2,1) = 0.0
  EXP_GRIDDO(:,3,1) = 1.4467853502767779E-005

  EXP_GRIDZO(:,1,1) = 1.4467853502767779E-005
  EXP_GRIDZO(:,2,1) = 0.0
  EXP_GRIDZO(:,3,1) = -1.4467853502767779E-005
  
  EXP_GRIDPO(:,1,1) = 293640925.08725864 
  EXP_GRIDPO(:,2,1) = 0.0
  EXP_GRIDPO(:,3,1) = -293640925.08725864 

  EXP_GRIDSO(:,1,1) = -293640925.08725864 
  EXP_GRIDSO(:,2,1) = 0.0
  EXP_GRIDSO(:,3,1) = 293640925.08725864 
# endif

  print *, "Testing SPTRUNV with LUV, LDZ, and LPS = TRUE"
  call SPTRUNV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,   &
              IDRTO,IMAXO,JMAXO,KMAX, &
              IPRIME,ISKIPI,JSKIPI,KSKIPI,  &
              ISKIPO,JSKIPO,KSKIPO,JCPU,GRIDUI,GRIDVI,  &
              LUV,GRIDUO,GRIDVO,LDZ,GRIDDO,GRIDZO,  &
              LPS,GRIDPO,GRIDSO)

  do I = 1, IMAXO
    do J = 1, JMAXO
      if (abs(GRIDUO(I, J, 1)-EXP_GRIDUO(I, J, 1)) > TOL) then
        print *, "GRIDUO(", I, ",", J, ",1) = ", GRIDUO(I, J, 1),   &
                  " EXPECTED: ", EXP_GRIDUO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDVO(I, J, 1)-EXP_GRIDVO(I, J, 1)) > TOL) then
        print *, "GRIDVO(", I, ",", J, ",1) = ", GRIDVO(I, J, 1),   &
                  " EXPECTED: ", EXP_GRIDVO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDDO(I, J, 1)-EXP_GRIDDO(I, J, 1)) > TOL) then
        print *, "GRIDDO(", I, ",", J, ",1) = ", GRIDDO(I, J, 1),   &
                  " EXPECTED: ", EXP_GRIDDO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDZO(I, J, 1)-EXP_GRIDZO(I, J, 1)) > TOL) then
        print *, "GRIDZO(", I, ",", J, ",1) = ", GRIDZO(I, J, 1),   &
                  " EXPECTED: ", EXP_GRIDZO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDPO(I, J, 1)-EXP_GRIDPO(I, J, 1)) > TOL) then
        print *, "GRIDPO(", I, ",", J, ",1) = ", GRIDPO(I, J, 1),   &
                  " EXPECTED: ", EXP_GRIDPO(I, J, 1)
        res = 1
      end if
      if (abs(GRIDSO(I, J, 1)-EXP_GRIDSO(I, J, 1)) > TOL) then
        print *, "GRIDSO(", I, ",", J, ",1) = ", GRIDSO(I, J, 1),   &
                  " EXPECTED: ", EXP_GRIDSO(I, J, 1)
        res = 1
      end if
    end do
  end do

  if (res .ne. 0) stop 2
  
  print *, "Success!"
end program test_sptrunv