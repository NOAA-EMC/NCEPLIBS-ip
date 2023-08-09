#if (LSIZE==D)
#define REALSIZE 8
#elif (LSIZE==4)
#define REALSIZE 4
#endif

program test_earth_radius
  use earth_radius_mod
  implicit none

  integer :: igdtlen
  parameter(igdtlen = 7)
  integer :: igdtmpl(igdtlen)
  real(KIND=REALSIZE) :: eccen_squared
  real(KIND=REALSIZE) :: radius
  
  print *,'Testing earth_radius_mod...'

  ! Case 0.
  igdtmpl(1) = 0
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 6367470.0) > epsilon(radius)) stop 2
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 3

  ! Case 1.
  igdtmpl(1) = 1
  igdtmpl(2) = 1
  igdtmpl(3) = 1
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 0.1) > epsilon(radius)) stop 4
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 5

  ! Another case 1.
  igdtmpl(1) = 1
  igdtmpl(2) = 2
  igdtmpl(3) = 63600000
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 636000.0) > epsilon(radius)) stop 4
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 5

  ! Case 2.
  igdtmpl(1) = 2
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 6378160) > epsilon(radius)) stop 6
  if (abs(eccen_squared - 6.7226700223333219E-003) > epsilon(eccen_squared)) stop 7

  ! Case 3.
  igdtmpl(1) = 3
  igdtmpl(4) = 1
  igdtmpl(5) = 1
  igdtmpl(6) = 1
  igdtmpl(7) = 1
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 100) > epsilon(radius)) stop 20
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 20

  ! Case 4.
  igdtmpl(1) = 4
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 6378137) > epsilon(radius)) stop 22
  if (abs(eccen_squared - 6.6943805181245335E-003) > epsilon(eccen_squared)) stop 23

  ! Case 5.
  igdtmpl(1) = 5
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 6378137) > epsilon(radius)) stop 24
  if (abs(eccen_squared - 6.6943799901300000E-003) > epsilon(eccen_squared)) stop 25

  ! Case 6.
  igdtmpl(1) = 6
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 6371229.0) > epsilon(radius)) stop 26
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 27

  ! Case 7.
  igdtmpl(1) = 7
  igdtmpl(4) = 1
  igdtmpl(5) = 1
  igdtmpl(6) = 1
  igdtmpl(7) = 1
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 0.1) > epsilon(radius)) stop 28
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 29

  ! Case 8.
  igdtmpl(1) = 8
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - 6371200.0) > epsilon(radius)) stop 28
  if (abs(eccen_squared - 0.0) > epsilon(eccen_squared)) stop 29

  ! Case 9.
  igdtmpl(1) = 9
  call earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
  if (abs(radius - (-9999)) > epsilon(radius)) stop 40
  if (abs(eccen_squared - (-9999)) > epsilon(eccen_squared)) stop 41
  
  print *,'SUCCESS!'

end program test_earth_radius
