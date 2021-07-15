 program driver

!----------------------------------------------------------------------------
! Test iplib routine ipxetas as follows:
!
!  1) Convert land/sea mask on a rotated lat/lon "e"-staggered "h" point
!     grid to a rotated lat/lon unstaggerd "full" grid.
!  2) Convert u-component wind on a rotated lat/lon "e"-staggered "v" point
!     grid to a rotated lat/lon unstaggerd "full" grid.
!  3) Convert land/sea mask on a rotated lat/lon unstaggered "full" grid
!     to a rotated lat/lon "e"-staggered "h" point grid.
!  4) Convert land/sea mask on a rotated lat/lon unstaggered "full" grid
!     to a rotated lat/lon "e"-staggered "v" point grid.
!
! All input files are in grib 2 format.  The converted output data is
! also grib 2.
!----------------------------------------------------------------------------

 use grib_mod

 implicit none

 character(len=2)                     :: option
 character(len=100)                   :: input_file, output_file

 integer                              :: lugi, iunit, iret
 integer                              :: j, jdisc, jpdtn, jgdtn, k
 integer                              :: jids(200), jgdt(200), jpdt(200)
 integer                              :: idir, npts_input, npts_output
 integer                              :: igdtnum_input, igdtnum_output, istat
 integer                              :: igdtlen
 integer(kind=4)                      :: i1
 integer        , allocatable, target :: igdtmpl_output(:)
 integer        , allocatable         :: igdtmpl_input(:)

 logical                              :: unpack
 logical*1, allocatable, target       :: bitmap_input(:), bitmap_output(:)
 
 real, allocatable, target            :: data_input(:), data_output(:)

 type(gribfield)                      :: gfld_input, gfld_output

 option="xx"
 i1=1
 call getarg(i1, option)

 select case (trim(option))
   case ('0')
     idir=0       ! "H" or "V" grid to full grid
   case ('-1')
     idir=-1      ! full grid to "H"
   case ('-2')
     idir=-2      ! full grid to "V"
   case default
     print*,"ERROR: ENTER VALID OPTION."
     stop 24
 end select

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 iunit=9
 call baopenr (iunit, input_file, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 j       = 0      ! search at beginning of file
 jdisc   = -1     ! search for any discipline
 jpdtn   = -1     ! search for any product definition template number
 jgdtn   =  1     ! search for grid definition template number 1 - rot lat/lon grid
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 unpack  = .true. ! unpack data
 lugi    = 0      ! no index file

 nullify(gfld_input%idsect)
 nullify(gfld_input%local)
 nullify(gfld_input%list_opt)
 nullify(gfld_input%igdtmpl)
 nullify(gfld_input%ipdtmpl)
 nullify(gfld_input%coord_list)
 nullify(gfld_input%idrtmpl)
 nullify(gfld_input%bmap)
 nullify(gfld_input%fld)

 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld_input, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

 print*,'- SCAN MODE IS: ', gfld_input%igdtmpl(19)
 print*,'- DIMENSION OF GRID: ', gfld_input%igdtmpl(8), gfld_input%igdtmpl(9)

 if (idir==0)then
   npts_output = gfld_input%igdtmpl(9) * (gfld_input%igdtmpl(8)*2-1)
 else
   npts_output = gfld_input%igdtmpl(9) * (gfld_input%igdtmpl(8)+1)/2
 endif

 npts_input = gfld_input%ngrdpts

 allocate(igdtmpl_input(gfld_input%igdtlen))
 igdtmpl_input=gfld_input%igdtmpl

 allocate(igdtmpl_output(gfld_input%igdtlen))
 igdtmpl_output=igdtmpl_input

 allocate(data_input(npts_input))
 data_input=gfld_input%fld

 allocate(bitmap_input(npts_input))
 bitmap_input=.true.

 allocate(data_output(npts_output))
 data_output=-999.

 allocate(bitmap_output(npts_output))
 bitmap_input=.true.

 igdtnum_input=gfld_input%igdtnum
 igdtlen=gfld_input%igdtlen

 call ipxetas(idir, igdtnum_input, igdtlen, igdtmpl_input, & 
              npts_input, bitmap_input, data_input, igdtnum_output,  &
              igdtmpl_output, npts_output, bitmap_output, data_output, istat)

 if (istat /= 0) then
   print*,"- PROBLEM IN IPXETAS. ISTAT IS: ", istat
   stop 5
 endif

 nullify(gfld_output%idsect)
 nullify(gfld_output%local)
 nullify(gfld_output%list_opt)
 nullify(gfld_output%igdtmpl)
 nullify(gfld_output%ipdtmpl)
 nullify(gfld_output%coord_list)
 nullify(gfld_output%idrtmpl)
 nullify(gfld_output%bmap)
 nullify(gfld_output%fld)

 gfld_output=gfld_input
 gfld_output%ngrdpts=npts_output
 gfld_output%igdtmpl=>igdtmpl_output
 gfld_output%idrtmpl=>gfld_input%idrtmpl
 gfld_output%idrtmpl=0
 gfld_output%idrtmpl(3)=gfld_input%idrtmpl(3)
 gfld_output%fld=>data_output

 if(any(.not.bitmap_output)) then
   gfld_output%ibmap=0
   gfld_output%bmap=>bitmap_output
 else
   gfld_output%ibmap=255
 endif

 output_file="./output.grb2"
 print*,"- OPEN AND WRITE FILE ", trim(output_file)
 iunit=10
 call baopenw(iunit, output_file, iret)
 if (iret /= 0) then
   print*,"- PROBLEM OPENING FILE. IRET: ", iret
   stop 21
 endif

 call putgb2(iunit, gfld_output, iret)
 if (iret /= 0) then
   print*,"- PROBLEM WRITING FILE. IRET: ", iret
   stop 22
 endif

 call baclose(iunit, iret)

 end program driver
