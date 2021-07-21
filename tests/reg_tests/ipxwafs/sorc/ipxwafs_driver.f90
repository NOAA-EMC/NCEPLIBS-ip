 program driver

!----------------------------------------------------------------------------
! Test iplib routine ipxwafs, ipxwafs2 and ipxwafs3 by transforming
! data on a thinned/full WAFS grid to its full/thinned counterpart.
!
! All input files are in grib 2 format.  The converted output data is
! also grib 2.
!----------------------------------------------------------------------------

 use grib_mod

 implicit none

 character(len=1)               :: option
 character(len=1), allocatable  :: cgrib(:)
 character(len=150)             :: input_file, output_file

 integer(kind=4)                :: i1
 integer                        :: listsec0(2), lengrib
 integer                        :: igds(5), lcgrib
 integer                        :: num_opt_pts, istat
 integer                        :: idir, km, m1, m2, igdtlen
 integer                        :: ib1, ib2
 integer                        :: j, jdisc, jpdtn, jgdtn, k
 integer                        :: jids(200), jgdt(200), jpdt(200)
 integer                        :: lugi, iret, iunit
 integer, allocatable           :: igdtmpl1(:)
 integer, allocatable, target   :: opt_pts(:), igdtmpl2(:)
 
 logical                        :: unpack
 logical*1, allocatable, target :: bmap1(:), bmap2(:)

 real, allocatable, target      :: data1(:), data2(:)

 type(gribfield)                :: gfld1, gfld2

 option="xx"
 i1=1
 call getarg(i1, option)

 select case (option)
   case ('1')
     print*,"- WILL TEST ROUTINE IPXWAFS"
   case ('2')
     print*,"- WILL TEST ROUTINE IPXWAFS2"
   case ('3')
     print*,"- WILL TEST ROUTINE IPXWAFS3"
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
 jgdtn   =  0     ! search for any grid definition template number 0 - lat/lon grid
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 unpack  = .true. ! unpack data
 lugi    = 0      ! no index file

 nullify(gfld1%idsect)
 nullify(gfld1%local)
 nullify(gfld1%list_opt)
 nullify(gfld1%igdtmpl)
 nullify(gfld1%ipdtmpl)
 nullify(gfld1%coord_list)
 nullify(gfld1%idrtmpl)
 nullify(gfld1%bmap)
 nullify(gfld1%fld)

 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld1, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

 if(associated(gfld1%list_opt)) then
   print*,"- THIS IS A THINNED GRID."
   print*,"- NUMBER OF OPTIONAL ENTRIES: ", gfld1%num_opt
   print*,"- NUMBER OF OCTETS REQUIRED OPTIONAL ENTRIES: ", gfld1%numoct_opt
   print*,"- INTERPRETATION OF OPTIONAL LIST: ", gfld1%interp_opt
   print*,"- NUMBER OF POINTS PER ROW: ", gfld1%list_opt
   print*,"- TOTAL NUMBER OF GRID POINTS: ", gfld1%ngrdpts
   idir=1
   m1=gfld1%ngrdpts
   m2=73*73
   num_opt_pts=gfld1%num_opt
   allocate(opt_pts(num_opt_pts))
   opt_pts=gfld1%list_opt
 else
   print*,"- THIS IS A FULL GRID."
   print*,"- TOTAL NUMBER OF GRID POINTS: ", gfld1%ngrdpts
   print*,"- NUMBER OF OPTIONAL ENTRIES: ", gfld1%num_opt
   idir=-1
   m1=gfld1%ngrdpts
   m2=3447
   num_opt_pts=73
   allocate(opt_pts(num_opt_pts))
   opt_pts=-999
 endif

 if (gfld1%ibmap==0) then
   print*,"- FILE HAS A BITMAP"
   ib1=1  ! has bitmap
   allocate(bmap1(m1))
   bmap1=gfld1%bmap
   allocate(bmap2(m2))
   bmap2=.false.
 else
   ib1=0
 endif
 ib2=0  ! output from ipxwafs2/3, 1 if bitmap

 igdtlen=gfld1%igdtlen
 allocate(igdtmpl1(igdtlen))
 igdtmpl1 = gfld1%igdtmpl

 allocate(igdtmpl2(igdtlen))
 igdtmpl2 = -999

 allocate(data1(m1))
 data1 = gfld1%fld

 allocate(data2(m2))
 data2=-9999.

 km=1
 if (idir == 1) then         ! create full grid
   if (gfld1%ibmap==0) then  ! input data has bitmap
     if (option=='2') then
       print*,"- CALL IPXWAFS2 TO CREATE FULL GRID"
       call ipxwafs2(idir, m1, m2, km, num_opt_pts, opt_pts, igdtlen, igdtmpl1, &
                     data1, ib1, bmap1, igdtmpl2, data2, ib2, bmap2, istat)
     elseif (option=='3') then
       print*,"- CALL IPXWAFS3 TO CREATE FULL GRID"
       call ipxwafs3(idir, m1, m2, km, num_opt_pts, opt_pts, igdtlen, igdtmpl1, &
                     data1, ib1, bmap1, igdtmpl2, data2, ib2, bmap2, istat)
     else  ! must use ipxwafs2 or ipxwafs3 for bitmap data
       print*,"- BAD OPTION"
       stop 67
     endif
   else  ! input data has no bitmap
     if (option=='1') then
       print*,"- CALL IPXWAFS TO CREATE FULL GRID"
       call ipxwafs(idir, m1, m2, km, num_opt_pts, opt_pts, igdtlen, igdtmpl1, &
                    data1, igdtmpl2, data2, istat)
     else  ! must use ipxwafs for non-bitmapped data
       print*,"- BAD OPTION"
       stop 68
     endif
   endif
 else                         ! create thin grid
   if (gfld1%ibmap==0) then   ! input data has bitmap
     if (option=='2') then
       print*,"- CALL IPXWAFS2 TO CREATE THIN GRID"
       call ipxwafs2(idir, m2, m1, km, num_opt_pts, opt_pts, igdtlen, igdtmpl2, &
                     data2, ib2, bmap2, igdtmpl1, data1, ib1, bmap1, istat)
     elseif (option=='3') then
       print*,"- CALL IPXWAFS3 TO CREATE THIN GRID"
       call ipxwafs3(idir, m2, m1, km, num_opt_pts, opt_pts, igdtlen, igdtmpl2, &
                     data2, ib2, bmap2, igdtmpl1, data1, ib1, bmap1, istat)
     else  ! must use ipxwafs2 or ipxwafs3 for bitmap data
       print*,"- BAD OPTION"
       stop 64
     endif
   else                     ! input data has no bitmap
     if (option=='1') then
       print*,"- CALL IPXWAFS TO CREATE THIN GRID"
       call ipxwafs(idir, m2, m1, km, num_opt_pts, opt_pts, igdtlen, igdtmpl2, &
                    data2, igdtmpl1, data1, istat)
     else  ! must use ipxwafs for non-bitmapped data
       print*,"- BAD OPTION"
       stop 63
     endif
   endif
 endif

 if (istat /= 0) then
   print*,'- ERROR IN IPXWAFS.  ISTAT: ', istat
   stop 17
 endif

 nullify(gfld2%idsect)
 nullify(gfld2%local)
 nullify(gfld2%list_opt)
 nullify(gfld2%igdtmpl)
 nullify(gfld2%ipdtmpl)
 nullify(gfld2%coord_list)
 nullify(gfld2%idrtmpl)
 nullify(gfld2%bmap)
 nullify(gfld2%fld)

! Section 0

 gfld2%version=gfld1%version
 gfld2%discipline=gfld1%discipline

! Section 1

 gfld2%idsect=>gfld1%idsect

! Section 2

 gfld2%locallen=0   ! local section not used

! Section 3

 gfld2%igdtlen=gfld1%igdtlen
 gfld2%griddef=gfld1%griddef
 gfld2%ngrdpts=m2
 if (idir == 1) then
   gfld2%num_opt=0
   gfld2%list_opt=>opt_pts
 else
   gfld2%num_opt=num_opt_pts
   gfld2%numoct_opt=1
   gfld2%list_opt=>opt_pts
   gfld2%interp_opt=1
 endif
 gfld2%igdtnum=gfld1%igdtnum
 gfld2%igdtmpl=>igdtmpl2

! Section 4

 gfld2%num_coord=gfld1%num_coord
 allocate(gfld2%coord_list(gfld2%num_coord))
 gfld2%ipdtnum=gfld1%ipdtnum
 gfld2%ipdtlen=gfld1%ipdtlen
 gfld2%ipdtmpl=>gfld1%ipdtmpl

! Section 5

 gfld2%idrtnum=gfld1%idrtnum
 gfld2%idrtlen=gfld1%idrtlen
 gfld2%idrtmpl=>gfld1%idrtmpl

 if (ib2==1)then
   gfld2%ibmap=0
   gfld2%bmap=>bmap2
 else
   gfld2%ibmap=255
   allocate(gfld2%bmap(gfld2%ngrdpts))
 endif

 gfld2%fld=>data2

 lcgrib=gfld2%ngrdpts*4
 allocate(cgrib(lcgrib))

 listsec0(1)=gfld2%discipline
 listsec0(2)=gfld2%version
 
 call gribcreate(cgrib,lcgrib,listsec0,gfld2%idsect,iret)
 if (iret /= 0) then
   print*,"- ERROR IN ROUTINE GRIBCREATE"
   stop 8
 endif

 igds(1)=gfld2%griddef
 igds(2)=gfld2%ngrdpts
 igds(3)=gfld2%numoct_opt
 igds(4)=gfld2%interp_opt
 igds(5)=gfld2%igdtnum

 call addgrid(cgrib,lcgrib,igds,gfld2%igdtmpl,gfld2%igdtlen,  &
              gfld2%list_opt,gfld2%num_opt,iret)
 if (iret /= 0) then
   print*,"- ERROR IN ROUTINE ADDGRID"
   stop 9
 endif

 call addfield(cgrib,lcgrib,gfld2%ipdtnum,gfld2%ipdtmpl,   &
                 gfld2%ipdtlen,gfld2%coord_list,gfld2%num_coord,  &
                 gfld2%idrtnum,gfld2%idrtmpl,gfld2%idrtlen,  &
                 gfld2%fld,gfld2%ngrdpts,gfld2%ibmap,gfld2%bmap, &
                 iret)
 if (iret /= 0) then
   print*,"- ERROR IN ROUTINE ADDFIELD"
   stop 10
 endif

 call gribend(cgrib,lcgrib,lengrib,iret)
 if (iret /= 0) then
   print*,"- ERROR IN ROUTINE GRIBEND"
   stop 10
 endif

 output_file="./output.grb2"
 print*,"- OPEN AND WRITE FILE ", trim(output_file)
 iunit=10
 call baopenw(iunit, output_file, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 call wryte(iunit, lengrib, cgrib)

 call baclose(iunit, iret)

 end program driver
