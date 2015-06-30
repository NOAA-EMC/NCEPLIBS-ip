 program driver

 use grib_mod

 implicit none

 character(len=1), allocatable :: cgrib(:)

 character*150   :: input_file, output_file

 integer         :: listsec0(2), lengrib
 integer(kind=4) :: igds(5), lcgrib
 integer         :: num_opt_pts
 integer         :: idir, km, m1, m2
 integer         :: j, jdisc, jpdtn, jgdtn, k
 integer         :: jids(200), jgdt(200), jpdt(200)
 integer         :: lugi, iret, iunit
 integer(kind=4), allocatable, target :: opt_pts(:), igdtmpl1(:), igdtmpl2(:)
 
 logical         :: unpack

 real, allocatable, target :: data1(:), data2(:)

 type(gribfield) :: gfld1, gfld2

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

 allocate(igdtmpl1(gfld1%igdtlen))
 igdtmpl1 = gfld1%igdtmpl

 allocate(igdtmpl2(gfld1%igdtlen))
 igdtmpl2 = -999

 allocate(data1(m1))
 data1 = gfld1%fld

 allocate(data2(m2))
 data2=-9999.

 km=1
 if (idir == 1) then
   call ipxwafs(idir, m1, m2, km, num_opt_pts, opt_pts, gfld1%igdtlen, igdtmpl1, &
                data1, igdtmpl2, data2, iret)
 else
   call ipxwafs(idir, m2, m1, km, num_opt_pts, opt_pts, gfld1%igdtlen, igdtmpl2, &
                data2, igdtmpl1, data1, iret)
 endif

 if (iret /= 0) then
   print*,'error in ipxwafs. iret: ', iret
   stop
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
 gfld2%coord_list=gfld1%coord_list
 gfld2%ipdtnum=gfld1%ipdtnum
 gfld2%ipdtlen=gfld1%ipdtlen
 gfld2%ipdtmpl=>gfld1%ipdtmpl

! Section 5

 gfld2%idrtnum=gfld1%idrtnum
 gfld2%idrtlen=gfld1%idrtlen
 gfld2%idrtmpl=>gfld1%idrtmpl

 gfld2%ibmap=255
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
                   opt_pts,gfld2%num_opt,iret)
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
