 module get_input_data       

!------------------------------------------------------------------------
! Read the data that will be interpolated.  Data is a global one-degree
! grid of albedo with no bitmap.
!------------------------------------------------------------------------

 integer, public                :: input_kgds(200)

 integer, parameter, public     :: i_input = 360
 integer, parameter, public     :: j_input = 180

 logical*1, public :: input_bitmap(i_input,j_input)

 real, public      :: input_data(i_input,j_input)

 data input_kgds /0,  360,  180, -89500, -180000, 128,  &
                  89500, 179000, 1000, 1000,  64,  0,   &
                  6*-1, 0, 255, 180*-1/
 contains

 subroutine read_input_data

 implicit none

 character*100      :: input_file

 integer            :: iret
 integer, parameter :: iunit=9

 real(kind=4)       :: dummy(i_input,j_input)

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 open(iunit, file=input_file, access='direct', recl=i_input*j_input*4,  &
      iostat=iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 read(iunit, rec=1, iostat=iret) dummy
 input_data=dummy

 if (iret /= 0) then
   print*,"- BAD READ OF DATA. IRET IS ", iret
   stop 4
 end if

 close (iunit)

 input_bitmap=.true.

 return

 end subroutine read_input_data

 end module get_input_data
