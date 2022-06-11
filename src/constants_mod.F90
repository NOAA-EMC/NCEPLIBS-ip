!> @file
!! @brief Constants for use in NCEPLIBS-ip.
!! @author Kyle Gerheiser, George Gayno
!! @date July 2021

!> Module containing common constants.
!!
!! @author Kyle Gerheiser, George Gayno
module constants_mod
  use iso_fortran_env, only: real64
  implicit none

  public
  
  real(real64), parameter :: pi=3.14159265358979d0 !< PI
  real(real64), parameter :: dpr=180d0/pi !< Radians to degrees
  real(real64), parameter :: pi2=pi/2.0d0 !< PI / 2.0
  real(real64), parameter :: pi4=pi/4.0d0 !< PI / 4.0
  real(real64), parameter :: RERTH_WGS84=6.378137E6 !< Radius of the Earth defined by WGS-84
  real(real64), parameter :: E2_WGS84 = 0.00669437999013 !< Eccentricity squared of Earth defined by WGS-84 
  
end module constants_mod

