 program driver

!--------------------------------------------------------------------
! regression test for iplib routines polateg0, polateg1 and polateg4
!--------------------------------------------------------------------

 use get_input_data

 implicit none

 call degrib_input_data

 call interp

 print*,"- NORMAL TERMINATION"

 stop 0
 end program driver
