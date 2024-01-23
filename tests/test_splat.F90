! This is a test from the NCEPLIBS-sp project.
!
! This test tests the splat() subrroutine.
!
! Kyle Gerheiser
program test_splat
  use iso_fortran_env, only: real64
  implicit none

  integer :: j, jj, jmax, ref_j(5)
  real :: slat(584), wlat(584), ref_slat(5), ref_wlat(5)
  real :: tini=1e-5
  
  jmax = 584  ! t382 grid

  call splat(0, jmax, slat, wlat)

  if (slat(1) /= 1d0) then
     error stop "slat(1) should equal 1.0"
  endif

  if(slat(jmax) /= -1d0) then
     error stop "slat(jmax) should equal -1.0"
  endif

  if(wlat(1) /= 0d0) then
     error stop "wlat(1) should equal 0.0"
  endif

  if(wlat(jmax) /= 0d0) then
     error stop "wlat(jmax) should equal 0.0"
  endif

  do j = 2, jmax-1
     if (slat(j) < slat(j+1)) then
        error stop "slat should be monotonically decreasing"
     endif
  end do

  call splat(256, jmax, slat, wlat)
  ref_j = (/1, 20, 100, 292, 584/)
  ref_slat = (/0.999996364, 0.994503140, 0.860138953, 2.68967217E-03, -0.999996364/)
  ref_wlat = (/1.25922097E-05, 5.63323090E-04, 2.74388562E-03, 5.37929870E-03, 1.25922097E-05/)

  do jj = 1, 5
     if (abs(ref_slat(jj)-slat(ref_j(jj))) .gt. tini) error stop "slat mismatch for IDRT=256"
     if (abs(ref_wlat(jj)-wlat(ref_j(jj))) .gt. tini) error stop "wlat mismatch for IDRT=256"
  enddo

end program test_splat
