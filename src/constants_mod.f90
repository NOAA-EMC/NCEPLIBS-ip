module constants_mod
  use iso_fortran_env, only: real64
  implicit none

  public
  
  real(real64), parameter :: pi=3.14159265358979d0
  real(real64), parameter :: dpr=180d0/pi
  real(real64), parameter :: pi2=pi/2.0d0
  real(real64), parameter :: pi4=pi/4.0d0
  real(real64), parameter :: RERTH_WGS84=6.378137E6
  real(real64), parameter :: E2_WGS84 = 0.00669437999013
  
end module constants_mod

