 module get_input_data       

!------------------------------------------------------------------------
! Read the data that will be interpolated.  Data is a global one-degree
! grid of 500 mb u and v-component winds.  There is no bitmap, so
! create one from the landmask record.
!------------------------------------------------------------------------

 integer,         allocatable, public   :: gdtmpl_input(:)
 integer, public                        :: gdtlen_input, gdtnum_input
 integer, public                        :: i_input, j_input

 logical*1, allocatable, public         :: input_bitmap(:)

 real, allocatable, public              :: input_u_data(:)
 real, allocatable, public              :: input_v_data(:)

 contains

 subroutine degrib_input_data

 use grib_mod

 implicit none

 character*100      :: input_file

 integer            :: iret, jdisc, jpdtn, jgdtn
 integer            :: jids(200), jgdt(200), jpdt(200)
 integer            :: j, k, lugi, numpts
 integer, parameter :: iunit=9

 logical            :: unpack

 type(gribfield)    :: gfld

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 call baopenr (iunit, input_file, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 nullify(gfld%idsect)
 nullify(gfld%local)
 nullify(gfld%list_opt)
 nullify(gfld%igdtmpl)
 nullify(gfld%ipdtmpl)
 nullify(gfld%coord_list)
 nullify(gfld%idrtmpl)
 nullify(gfld%bmap)
 nullify(gfld%fld)

 lugi    = 0
 j       = 0      ! search at beginning of file
 jdisc   = 0      ! search for discipline; 0 - meteorological products
 jpdtn   = 0      ! search for any product definition template 4.0 - anal at 1 lvl
 jgdtn   = -1     ! search for any grid definition template number
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 jpdt(1) = 2      ! search for parameter category - momentum
 jpdt(2) = 2      ! search for parameter number - u wind
 unpack  = .true. ! unpack data

 print*,"- DEGRIB U-WIND DATA"
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

 allocate(input_u_data(numpts))
 input_u_data=gfld%fld

 if (associated(gfld%idsect)) deallocate(gfld%idsect)
 if (associated(gfld%local)) deallocate(gfld%local)
 if (associated(gfld%list_opt)) deallocate(gfld%list_opt)
 if (associated(gfld%igdtmpl)) deallocate(gfld%igdtmpl)
 if (associated(gfld%ipdtmpl)) deallocate(gfld%ipdtmpl)
 if (associated(gfld%coord_list)) deallocate(gfld%coord_list)
 if (associated(gfld%idrtmpl)) deallocate(gfld%idrtmpl)
 if (associated(gfld%bmap)) deallocate(gfld%bmap)
 if (associated(gfld%fld)) deallocate(gfld%fld)

 lugi    = 0
 j       = 0      ! search at beginning of file
 jdisc   = 0      ! search for discipline; 0 - meteorological products
 jpdtn   = 0      ! search for any product definition template 4.0 - anal at 1 lvl
 jgdtn   = -1     ! search for any grid definition template number
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 jpdt(1) = 2      ! search for parameter category - momentum
 jpdt(2) = 3      ! search for parameter number - v wind
 unpack  = .true. ! unpack data

 print*,"- DEGRIB V-WIND DATA"
 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 5
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 allocate(input_v_data(numpts))
 input_v_data=gfld%fld

 if (associated(gfld%idsect)) deallocate(gfld%idsect)
 if (associated(gfld%local)) deallocate(gfld%local)
 if (associated(gfld%list_opt)) deallocate(gfld%list_opt)
 if (associated(gfld%igdtmpl)) deallocate(gfld%igdtmpl)
 if (associated(gfld%ipdtmpl)) deallocate(gfld%ipdtmpl)
 if (associated(gfld%coord_list)) deallocate(gfld%coord_list)
 if (associated(gfld%idrtmpl)) deallocate(gfld%idrtmpl)
 if (associated(gfld%bmap)) deallocate(gfld%bmap)
 if (associated(gfld%fld)) deallocate(gfld%fld)

 lugi    = 0
 j       = 0      ! search at beginning of file
 jdisc   = 2      ! search for discipline; 2 - land products
 jpdtn   = 0      ! search for any product definition template 4.0 - anal at 1 lvl
 jgdtn   = -1     ! search for any grid definition template number
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 jpdt(1) = 0      ! search for parameter category - veg/biomass
 jpdt(2) = 0      ! search for parameter number - land cover
 unpack  = .true. ! unpack data

 print*,"- DEGRIB MASK DATA"
 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

 allocate(input_bitmap(numpts))

 input_bitmap=.false.
 do j = 1, numpts
   if (gfld%fld(j) < 0.5) then  ! define wind over water
     input_bitmap(j)=.true.
   endif
 enddo

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
