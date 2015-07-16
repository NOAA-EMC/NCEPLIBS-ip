 module get_input_data       

!------------------------------------------------------------------------
! Read the data that will be interpolated.  Data is a global one-degree
! grid of land substrate temperatures.  Water points are masked out
! with a bitmap.
!------------------------------------------------------------------------

 integer(kind=4), allocatable, public :: gdtmpl_input(:)
 integer, public                      :: gdtlen_input, gdtnum_input
 integer, public                      :: i_input, j_input

 logical*1, allocatable, public       :: input_bitmap(:)

 real, allocatable, public            :: input_data(:)

 contains

 subroutine degrib_input_data

 use grib_mod  ! ncep grib 2 library

 implicit none

 character*100      :: input_file

 integer            :: j, jdisc, jpdtn, jgdtn, k
 integer            :: jids(200), jgdt(200), jpdt(200)
 integer            :: iret
 integer            :: lugi, numpts
 integer, parameter :: iunit=9

 logical            :: unpack

 type(gribfield)    :: gfld

 nullify(gfld%idsect)
 nullify(gfld%local)
 nullify(gfld%list_opt)
 nullify(gfld%igdtmpl)
 nullify(gfld%ipdtmpl)
 nullify(gfld%coord_list)
 nullify(gfld%idrtmpl)
 nullify(gfld%bmap)
 nullify(gfld%fld)

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 call baopenr (iunit, input_file, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 j       = 0      ! search at beginning of file
 jdisc   = -1     ! search for any discipline
 jpdtn   = -1     ! search for any product definition template number
 jgdtn   =  0     ! search for any grid definition template number 0 - lat/lon grid
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 unpack  = .true. ! unpack data
 lugi    = 0      ! no index file

 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 gdtnum_input = gfld%igdtnum
 gdtlen_input = gfld%igdtlen
 allocate(gdtmpl_input(gdtlen_input))
 gdtmpl_input=gfld%igdtmpl

 i_input = gdtmpl_input(8)
 j_input = gdtmpl_input(9)
 numpts=i_input*j_input

 allocate(input_data(numpts))
 input_data=gfld%fld
 allocate(input_bitmap(numpts))
 input_bitmap=gfld%bmap

 call baclose (iunit, iret)

 if (associated(gfld%idsect)) deallocate(gfld%idsect)
 if (associated(gfld%local)) deallocate(gfld%local)
 if (associated(gfld%list_opt)) deallocate(gfld%list_opt)
 if (associated(gfld%igdtmpl)) deallocate(gfld%igdtmpl)
 if (associated(gfld%ipdtmpl)) deallocate(gfld%ipdtmpl)
 if (associated(gfld%coord_list)) deallocate(gfld%coord_list)
 if (associated(gfld%idrtmpl)) deallocate(gfld%idrtmpl)
 if (associated(gfld%bmap)) deallocate(gfld%bmap)
 if (associated(gfld%fld)) deallocate(gfld%fld)
 
 return

 end subroutine degrib_input_data

 end module get_input_data
