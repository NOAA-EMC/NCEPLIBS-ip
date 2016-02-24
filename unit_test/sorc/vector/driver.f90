 program driver

!-----------------------------------------------------------------------
! Interpolate a global lat/lon grid of vector wind to several
! grids of various projections using all ipolatev 
! interpolation options.
!-----------------------------------------------------------------------

 use get_input_data

 implicit none

 call input_data

 call interp

!print*,"- NORMAL TERMINATION"

 stop
 end program driver
