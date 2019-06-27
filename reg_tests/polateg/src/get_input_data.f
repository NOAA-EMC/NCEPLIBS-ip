 module get_input_data       

!------------------------------------------------------------------------
! read the data that will be interpolated.  data is a global one-degree
! grid of land substrate temperatures.  water points are masked out
! with a bitmap.
!------------------------------------------------------------------------

 integer, public                :: input_kgds(200)
 integer, public                :: i_input, j_input

 logical*1, allocatable, public :: input_bitmap(:,:)

 real, allocatable, public      :: input_data(:,:)

 contains

 subroutine degrib_input_data

 implicit none

 character*100      :: input_file

 integer            :: iret, lskip, jpds(200), jgds(200), kpds(200)
 integer            :: lugi, numbytes, numpts, message_num
 integer, parameter :: iunit=9

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 call baopenr (iunit, input_file, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 lugi       = 0
 lskip      = -1
 jpds       = -1
 jgds       = -1
 jpds(5)    = 11    ! temperature
 kpds       = jpds
 input_kgds = jgds

 print*,"- GET GRIB HEADER"
 call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
             numpts, message_num, kpds, input_kgds, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF HEADER. IRET IS ", iret
   stop 3
 else
   print*,'- SUCCESSFULL DEGRIB OF HEADER'
 end if

 i_input = input_kgds(2)
 j_input = input_kgds(3)

 allocate(input_data(i_input,j_input))
 allocate(input_bitmap(i_input,j_input))

 print*,"- DEGRIB DATA"
 call getgb(iunit, lugi, numpts, lskip, jpds, jgds, &
            numpts, lskip, kpds, input_kgds, input_bitmap, input_data, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

 return

 end subroutine degrib_input_data

 end module get_input_data
