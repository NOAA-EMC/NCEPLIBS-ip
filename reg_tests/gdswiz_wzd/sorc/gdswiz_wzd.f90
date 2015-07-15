 program gdswzd_driver
 
!------------------------------------------------------------------
! Test gdswzd routines for all map projections.  
!
! Routine is called twice to test both the (1) i/j to lat/lon 
! and the (2) lat/lon to i/j conversions.  This should be
! reversable.  If not, a warning message is printed to
! standard output.  (1) is invoked by setting the IOPT argument
! to '0' and (2) is invoked by setting it to '-1'.
!
! This program takes one argument:  The grid number.  The valid
! grids are defined by data statements below.  They are:
!
! grid #       description
! ======       ===========
! 003          one-degree global lat/lon (ncep grid 3)
! 008          mercator (ncep grid 8)
! 127          t254 gaussian (ncep grid 127)
! 203h         rotated lat/lon e-staggered (number meaningless)
!              this is the old 12km eta grid - 'h' points
! 203v         rotated lat/lon e-staggered (number meaningless)
!              this is the old 12km eta grid - 'v' points
! 205h         rotated lat/lon b-staggered (number meaningless)
!              this is the 12km nam grid - 'h' points
! 205v         rotated lat/lon b-staggered (number meaningless)
!              this is the 12km nam grid - 'v' points
! 212          nh polar stereographic, spherical earth (number meaningless)
! 213          sh polar stereographic, spherical earth (number meaningless)
! 218          lambert conformal (ncep grid 218)
! 222          nh polar stereographic, elliptical earth (number meaningless)
!
! All computed fields are output to a direct access file so
! they may be compared for bit-idenicalness with other test runs,
! or so they may be visualized.  Separate files are written
! to output the data from each call to gdswzd.  The
! file naming convention is:
!
! grid${gridnum}.iopt${0 or m1}.bin
!------------------------------------------------------------------

 implicit none

 character*4   :: grid
 character*3   :: routine
 character*100 :: outfile

 integer*4 :: i1
 integer   :: nret, lrot, lmap, iopt, npts, imdl, jmdl
 integer   :: i, j, n, iret, kgds(200), nscan, kscan, is1, nm
 integer   :: ii, jj, iii, jjj, badpts, i_offset_odd

 real :: diff, fill, maxdiffx, maxdiffy
 real, allocatable :: xpts(:,:), ypts(:,:)
 real, allocatable :: rlat(:,:), rlon(:,:)

 real, allocatable :: crot(:,:), srot(:,:)
 real, allocatable :: xlon(:,:), xlat(:,:)
 real, allocatable :: ylon(:,:), ylat(:,:), area(:,:)

!-----------------------------------------------------------------------
! the grids that will be tested.  they are defined via the grib 2 grid
! definition template, which corresponds to the gfld%igdtmpl component
! of the ncep g2 library.
!-----------------------------------------------------------------------

 integer(kind=4), parameter   :: missing=b'11111111111111111111111111111111'

 integer(kind=4), allocatable :: igdtmpl(:)
 integer                      :: igdtlen
 integer                      :: igdtnum

 integer, parameter :: igdtlen3 = 19 ! ncep grid3; one-degree lat/lon
 integer(kind=4)    :: igdtmpl3(igdtlen3)
 data igdtmpl3 / 6, 255, missing, 255, missing, 255, missing, 360, 181, 0, missing, &
                 90000000, 0, 56, -90000000, 359000000, 1000000, 1000000, 0 /

 integer, parameter :: igdtlen8 = 19 ! ncep grid8; mercator 
 integer(kind=4)    :: igdtmpl8(igdtlen8)  
 data igdtmpl8 / 6, 255, missing, 255, missing, 255, missing, 116, 44, &
                 -48670000, 3104000, 56, 22500000, 61050000, 0, 64, 0, &
                  318830000, 318830000/

 integer, parameter :: igdtlen127=19  ! t254 gaussain
 integer(kind=4)    :: igdtmpl127(igdtlen127)
 data igdtmpl127 /6, 255, missing, 255, missing, 255, missing, 768, 384, &
                  0, missing, 89642000, 0, 48, -89642000, 359531000,  &
                  469000, 192, 0/

 integer, parameter :: igdtlen203h=22 ! 12km eta, h pts
 integer(kind=4)    :: igdtmpl203h(igdtlen203h)
 data igdtmpl203h/6, 255, missing, 255, missing, 255, missing, 669, 1165, &
                  0, missing, -7450000, 215860000, 56, 44560100, 14744800, &
                  179641, 77320, 68, -36000000, 254000000, 0 /

 integer, parameter :: igdtlen203v=22 ! 12km eta, v pts
 integer(kind=4)    :: igdtmpl203v(igdtlen203v)
 data igdtmpl203v/6, 255, missing, 255, missing, 255, missing, 669, 1165, &
                  0, missing, -7450000, 215860000, 56, 44560100, 14744800, &
                  179641, 77320, 72, -36000000, 254000000, 0 /

 integer, parameter :: igdtlen212=18 ! nh polar, spherical earth
 integer(kind=4)    :: igdtmpl212(igdtlen212)
 data igdtmpl212 /6, 255, missing, 255, missing, 255, missing, 512, 512, &
                  -20826000, 145000000, 56, 60000000, 280000000, 47625000, 47625000, &
                  0, 0/ 

 integer, parameter :: igdtlen222=18 ! nh polar, elliptical earth
 integer(kind=4)    :: igdtmpl222(igdtlen222)
 data igdtmpl222 /5, 255, missing, 255, missing, 255, missing, 512, 512, &
                  -20826000, 145000000, 56, 60000000, 280000000, 47625000, 47625000, &
                  0, 0/ 

 integer, parameter :: igdtlen213=18 ! sh polar, spherical earth
 integer(kind=4)    :: igdtmpl213(igdtlen213)
 data igdtmpl213 /6, 255, missing, 255, missing, 255, missing, 512, 512, &
                  20826000, 235000000, 56, -60000000, 100000000, 47625000, 47625000, &
                  128, 0/ 

 integer, parameter :: igdtlen205h=22 ! 12km nam, h pts
 integer(kind=4)    :: igdtmpl205h(igdtlen205h)
 data igdtmpl205h/6, 255, missing, 255, missing, 255, missing, 954, 835, &
                  0, missing, -7491200, 215866300, 56, 44539600, 14801500, &
                  126000, 108000, 64, -36000000, 254000000, 0 /

 integer, parameter :: igdtlen205v=22 ! 12km nam, v pts
 integer(kind=4)    :: igdtmpl205v(igdtlen205v)
 data igdtmpl205v/6, 255, missing, 255, missing, 255, missing, 954, 835, &
                  0, missing, -7491200, 215866300, 56, 44539600, 14801500, &
                  126000, 108000, 78, -36000000, 254000000, 0 /

 integer, parameter :: igdtlen218 = 22 ! ncep grid 218; lambert conf
 integer(kind=4)    :: igdtmpl218(igdtlen218)  
 data igdtmpl218 / 6, 255, missing, 255, missing, 255, missing, 614, 428, &
                   12190000, 226541000, 56, 25000000, 265000000, &
                   12191000, 12191000, 0, 64, 25000000, 25000000, -90000000, 0/

 i1=1
 call getarg(i1,grid)

 kgds=0
 select case (trim(grid))
   case ('3')
     igdtnum=0
     igdtlen=igdtlen3
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl3
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('8')
     igdtnum=10
     igdtlen=igdtlen8
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl8
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('203h')
     igdtnum=1
     igdtlen=igdtlen203h
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl203h
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('203v')
     igdtnum=1
     igdtlen=igdtlen203v
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl203v
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('127')
     igdtnum=40
     igdtlen=igdtlen127
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl127
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('205h')
     igdtnum=1
     igdtlen=igdtlen205h
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl205h
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('205v')
     igdtnum=1
     igdtlen=igdtlen205v
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl205v
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('212')
     igdtnum=20
     igdtlen=igdtlen212
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl212
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('222')
     igdtnum=20
     igdtlen=igdtlen222
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl222
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('213')
     igdtnum=20
     igdtlen=igdtlen213
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl213
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case ('218')
     igdtnum=30
     igdtlen=igdtlen218
     allocate(igdtmpl(igdtlen))
     igdtmpl=igdtmpl218
     imdl=igdtmpl(8)
     jmdl=igdtmpl(9)
   case default
     print*,"ENTER GRID NUMBER TO TEST"
     stop
 end select

 print*,"PROCESS GRID ", grid

!---------------------------------------------------------------------------
! first, call gdswzd to calculate lat/lon for each grid point.
!---------------------------------------------------------------------------

 fill = -9999.

 allocate (xpts(imdl,jmdl),ypts(imdl,jmdl))
 allocate (rlat(imdl,jmdl),rlon(imdl,jmdl))
 allocate (crot(imdl,jmdl),srot(imdl,jmdl))
 allocate (xlon(imdl,jmdl),xlat(imdl,jmdl))
 allocate (ylon(imdl,jmdl),ylat(imdl,jmdl),area(imdl,jmdl))

 xpts = fill
 ypts = fill
 rlat = fill
 rlon = fill
 crot = fill
 srot = fill
 xlon = fill
 xlat = fill
 ylon = fill
 ylat = fill
 area = fill

 lrot = 1
 lmap = 1
 iopt = 0
 npts = imdl * jmdl

 call gdswzd(igdtnum, igdtmpl, igdtlen, iopt, npts, fill, xpts, ypts, rlon, rlat, &
             nret, lrot, crot, srot, lmap, xlon, xlat, ylon, ylat, area)

 if (nret /= npts) then
   print*,'ERROR. WRONG NUMBER OF POINTS RETURNED ',nret,npts
   stop 33
 endif

 print*,'LAT/LON POINT(1,1):   ',rlat(1,1),rlon(1,1)
 print*,'LAT/LON POINT(1,JM):  ',rlat(1,jmdl),rlon(1,jmdl)
 print*,'LAT/LON POINT(IM,1):  ',rlat(imdl,1),rlon(imdl,1)
 print*,'LAT/LON POINT(IM,JM): ',rlat(imdl,jmdl),rlon(imdl,jmdl)

 outfile = "./grid" // trim(grid) // ".iopt0.bin"
 open (9, file=trim(outfile), access='direct', err=55, recl=imdl*jmdl*4)
 write(9, rec=1, err=55) real(rlat,4)
 write(9, rec=2, err=55) real(rlon,4)
 write(9, rec=3, err=55) real(xpts,4)
 write(9, rec=4, err=55) real(ypts,4)
 write(9, rec=5, err=55) real(crot,4)
 write(9, rec=6, err=55) real(srot,4)
 write(9, rec=7, err=55) real(xlon,4)
 write(9, rec=8, err=55) real(xlat,4)
 write(9, rec=9, err=55) real(ylon,4)
 write(9, rec=10, err=55) real(ylat,4)
 write(9, rec=11, err=55) real(area,4)
 close (9)

! the first call to gdswzd computed the lat/lon at each point.  now,
! given that lat/lon, compute the i/j.  it should be reversable.

 iopt = -1
 xpts=fill
 ypts=fill

 call gdswzd(igdtnum, igdtmpl, igdtlen, iopt, npts, fill, xpts, ypts, rlon, rlat, &
             nret, lrot, crot, srot, lmap, xlon, xlat, ylon, ylat, area)

 if (nret /= npts) then
   print*,'ERROR. WRONG NUMBER OF POINTS RETURNED ',nret,npts
   stop 34
 endif

 outfile = "./grid" // trim(grid) // ".ioptm1.bin"
 open (49, file=trim(outfile), access='direct', err=55, recl=imdl*jmdl*4)
 write(49, rec=1, err=55) real(rlat,4)
 write(49, rec=2, err=55) real(rlon,4)
 write(49, rec=3, err=55) real(xpts,4)
 write(49, rec=4, err=55) real(ypts,4)
 write(49, rec=5, err=55) real(crot,4)
 write(49, rec=6, err=55) real(srot,4)
 write(49, rec=7, err=55) real(xlon,4)
 write(49, rec=8, err=55) real(xlat,4)
 write(49, rec=9, err=55) real(ylon,4)
 write(49, rec=10, err=55) real(ylat,4)
 write(49, rec=11, err=55) real(area,4)
 close (49)

!------------------------------------------------------------------------------
! did the second call to gdswzd work?
!
! note: the gdswzdcb routine works on a grid that is tilted 45 degrees, so
! the internal i/j's do not match the normal convention.  account for this.
!------------------------------------------------------------------------------

 maxdiffx = -99999.
 maxdiffy = -99999.

 if (grid == "203h" .or. grid == "203v") then
   I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
   kscan=I_OFFSET_ODD
   if(kscan.eq.0) THEN
     is1=(jmdl+1)/2
   else
     is1=jmdl/2
   endif
   nm=imdl*jmdl
   nscan=mod(IGDTMPL(19)/32,2)
   badpts = 0
   do iii = 1, imdl    ! here iii/jjj are the conventional i/j
   do jjj = 1, jmdl
     n = (jjj-1)*imdl + iii
     if(nscan.eq.0) then
       j=(n-1)/imdl+1
       i=(n-imdl*(j-1))*2-mod(j+kscan,2)
     else
       i=(n-1)/jmdl+1
       j=(n-jmdl*(i-1))*2-mod(i+kscan,2)
     endif
     ii = (is1+(i-(j-kscan))/2)  ! here, ii/jj are from the tilted grid
     jj = ((i+(j-kscan))/2)      ! the tilted value is stored in xpts/ypts.
     diff = abs(float(ii)-xpts(iii,jjj))
     maxdiffx = max(maxdiffx, diff)
     if ( diff > .01) then
       print*,'BAD X POINT: ',iii,jjj,ii,jj,xpts(iii,jjj),ypts(iii,jjj)
       badpts=badpts+1
     endif 
     diff = abs(float(jj)-ypts(iii,jjj))
     maxdiffy = max(maxdiffy, diff)
     if ( diff > .01) then
       print*,'BAD Y POINT: ',iii,jjj,ii,jj,xpts(iii,jjj),ypts(iii,jjj)
       badpts=badpts+1
     endif 
   enddo
   enddo
 else
   badpts=0
   do j = 1, jmdl
   do i = 1, imdl
     diff = abs(float(i)-xpts(i,j))
     maxdiffx = max(maxdiffx, diff)
     if ( diff > .01) then
       print*,'BAD X POINT: ',i,j,xpts(i,j),ypts(i,j)
       badpts=badpts+1
     endif 
     diff = abs(float(j)-ypts(i,j))
     maxdiffy = max(maxdiffy, diff)
     if ( diff  > .01) then
       print*,'BAD Y POINT: ',i,j,xpts(i,j),ypts(i,j)
       badpts=badpts+1
     endif 
   enddo
   enddo
 endif

 if (badpts > 0) print*,"NUMBER OF BAD POINTS: ", badpts
 print*,'MAX DIFFERENCES IN X/Y CALCULATIONS: ', maxdiffx, maxdiffy

 deallocate (igdtmpl)
 deallocate (xpts,ypts,rlat,rlon,crot,srot,xlon,xlat,ylon,ylat,area)

 98 continue

 print*,'NORMAL TERMINATION'

 stop

 55 continue
 print*,'ERROR WRITING OUTPUT FILE.'
 stop 44

 end program gdswzd_driver
