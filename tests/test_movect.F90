! This is a test for the NCEPLIBS-ip library.
!
! Tests the MOVECT subroutine which computes rotation parameters
! to move a vector along a great circle.
!
! Hang Lei, Aug 2026

program test_movect
  implicit none

  real :: crot, srot
  real, parameter :: tol = 1.0e-6

  print *, 'Testing MOVECT...'

  ! Test 1: Same point - nearly coincident, expect cos(dlon)~1, sin~0
  call movect(45.0, 90.0, 45.0, 90.0, crot, srot)
  if (abs(crot - 1.0) > tol) then
    print *, 'FAILED test 1: crot =', crot, ' expected ~1.0'
    stop 1
  end if
  if (abs(srot) > tol) then
    print *, 'FAILED test 1: srot =', srot, ' expected ~0.0'
    stop 2
  end if

  ! Test 2: Points on equator, same latitude - crot^2 + srot^2 = 1
  call movect(0.0, 0.0, 0.0, 90.0, crot, srot)
  if (abs(crot**2 + srot**2 - 1.0) > tol) then
    print *, 'FAILED test 2: crot^2 + srot^2 =', crot**2 + srot**2
    stop 3
  end if

  ! Test 3: North pole to equator
  call movect(90.0, 0.0, 0.0, 0.0, crot, srot)
  if (abs(crot**2 + srot**2 - 1.0) > tol) then
    print *, 'FAILED test 3: crot^2 + srot^2 =', crot**2 + srot**2
    stop 4
  end if

  ! Test 4: Equator to North pole - by symmetry same magnitude as test 3
  call movect(0.0, 0.0, 90.0, 0.0, crot, srot)
  if (abs(crot**2 + srot**2 - 1.0) > tol) then
    print *, 'FAILED test 4: crot^2 + srot^2 =', crot**2 + srot**2
    stop 5
  end if

  ! Test 5: Known values - flat=35N, flon=90W, tlat=35N, tlon=90W (same)
  ! Should be the coincident case: crot ~ cos(0)=1, srot ~ sin(0)*sin(lat)
  call movect(35.0, -90.0, 35.0, -90.0, crot, srot)
  if (abs(crot - 1.0) > tol) then
    print *, 'FAILED test 5: crot =', crot
    stop 6
  end if

  ! Test 6: Moving along a meridian (same longitude, different latitudes)
  ! No rotation around the great circle when travelling along a meridian.
  call movect(0.0, 45.0, 45.0, 45.0, crot, srot)
  if (abs(crot**2 + srot**2 - 1.0) > tol) then
    print *, 'FAILED test 6: crot^2 + srot^2 =', crot**2 + srot**2
    stop 7
  end if

  print *, 'SUCCESS!'

end program test_movect
