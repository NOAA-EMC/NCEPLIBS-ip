!> @file
!! @brief Constants for use in NCEPLIBS-ip.
!! @author Kyle Gerheiser, George Gayno, Alex Richert
!! @date July 2021

!> @brief Module containing common constants.
!!
!! @author Kyle Gerheiser, George Gayno, Alex Richert
module constants_mod
  use iso_fortran_env, only: real64
  implicit none

  public
  
  real(real64), parameter :: pi=3.14159265358979_real64 !< PI
  real(real64), parameter :: dpr=180_real64/pi !< Radians to degrees
  real(real64), parameter :: pi2=pi/2.0 !< PI / 2.0
  real(real64), parameter :: pi4=pi/4.0 !< PI / 4.0
  real(real64), parameter :: RERTH_WGS84=6.378137E6_real64 !< Radius of the Earth defined by WGS-84
  real(real64), parameter :: E2_WGS84 = 0.00669437999013_real64 !< Eccentricity squared of Earth defined by WGS-84 
  
end module constants_mod

