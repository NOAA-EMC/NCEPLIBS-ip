! This is a test for the NCEPLIBS-ip library.
!
! Tests ip_constants_mod for expected constant values.
!
! Hang Lei, Aug 2026

program test_ip_constants
  use ip_constants_mod
  implicit none

  real, parameter :: tol = 1.0e-6

  print *, 'Testing ip_constants_mod...'

  ! pi should be close to the mathematical value
  if (abs(pi - 3.14159265358979) > tol) then
    print *, 'FAILED: pi =', pi
    stop 1
  end if

  ! dpr = 180/pi
  if (abs(dpr - 180.0/pi) > tol) then
    print *, 'FAILED: dpr =', dpr
    stop 2
  end if

  ! pi2 = pi/2
  if (abs(pi2 - pi/2.0) > tol) then
    print *, 'FAILED: pi2 =', pi2
    stop 3
  end if

  ! pi4 = pi/4
  if (abs(pi4 - pi/4.0) > tol) then
    print *, 'FAILED: pi4 =', pi4
    stop 4
  end if

  ! WGS-84 Earth radius in meters
  if (abs(RERTH_WGS84 - 6.378137e6) > 1.0) then
    print *, 'FAILED: RERTH_WGS84 =', RERTH_WGS84
    stop 5
  end if

  ! WGS-84 eccentricity squared
  if (abs(E2_WGS84 - 0.00669437999013) > tol) then
    print *, 'FAILED: E2_WGS84 =', E2_WGS84
    stop 6
  end if

  ! DLON_OFFSET
  if (abs(DLON_OFFSET - 1.0e-2) > tol) then
    print *, 'FAILED: DLON_OFFSET =', DLON_OFFSET
    stop 7
  end if

  ! Sanity: 2*pi2 == pi
  if (abs(2.0*pi2 - pi) > tol) then
    print *, 'FAILED: 2*pi2 /= pi'
    stop 8
  end if

  ! Sanity: 4*pi4 == pi
  if (abs(4.0*pi4 - pi) > tol) then
    print *, 'FAILED: 4*pi4 /= pi'
    stop 9
  end if

  print *, 'SUCCESS!'

end program test_ip_constants
