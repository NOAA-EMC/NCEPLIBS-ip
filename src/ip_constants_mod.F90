!> @file
!! @brief Constants for use in NCEPLIBS-ip.
!! @author Kyle Gerheiser, George Gayno, Alex Richert
!! @date July 2021

!> @brief Module containing common constants.
!!
!! @author Kyle Gerheiser, George Gayno, Alex Richert
module ip_constants_mod
  implicit none

  public
  
  real, parameter :: pi=3.14159265358979 !< PI
  real, parameter :: dpr=180.0/pi !< Radians to degrees
  real, parameter :: pi2=pi/2.0 !< PI / 2.0
  real, parameter :: pi4=pi/4.0 !< PI / 4.0
  real, parameter :: RERTH_WGS84=6.378137E6 !< Radius of the Earth defined by WGS-84
  real, parameter :: E2_WGS84 = 0.00669437999013 !< Eccentricity squared of Earth defined by WGS-84
  
end module ip_constants_mod

