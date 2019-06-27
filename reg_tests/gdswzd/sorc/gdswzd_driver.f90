 program gdswzd_driver
 
!------------------------------------------------------------------
! Test gdswzd routines for all map projections.  
!
! Routines are called twice to test both the (1) i/j to lat/lon 
! and the (2) lat/lon to i/j conversions.  This should be
! reversable.  If not, a warning message is printed to
! standard output.  (1) is invoked by setting the IOPT argument
! to '0' and (2) is invoked by setting it to '-1'.
!
! This program takes one argument: the grid number.  The valid 
! grids are defined by data statements below.  They are:
!
! grid #       description
! ======       ===========
! 003          one-degree global lat/lon (ncep grid 3)
! 008          mercator (ncep grid 8)
! 127          t254 gaussian (ncep grid 127)
! 203          rotated lat/lon e-staggered (number refers to gds octet 6)
!              tests routine gdswzdcb
! 205          rotated lat/lon b-staggered (number refers to gds octet 6)
!              tests routine gdswzdcd
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

 use gdswzd_mod

 implicit none

 character*3   :: grid
 character*100 :: outfile

 integer*4 :: i1
 integer   :: nret, lrot, lmap, iopt, npts, imdl, jmdl
 integer   :: i, j, n, iret, kgds(200), nscan, kscan, is1, nm
 integer   :: ii, jj, iii, jjj, badpts

 real :: diff, fill, maxdiffx, maxdiffy
 real, allocatable :: xpts(:,:), ypts(:,:)
 real, allocatable :: rlat(:,:), rlon(:,:)

 real, allocatable :: crot(:,:), srot(:,:)
 real, allocatable :: xlon(:,:), xlat(:,:)
 real, allocatable :: ylon(:,:), ylat(:,:), area(:,:)

! the grids that will be tested.

 integer :: grd3(200)    ! ncep grid3; one-degree lat/lon, for gdswzd00 routine
 data grd3 / 0, 360, 181, 90000, 0, 128, -90000,  &
            -1000, 1000, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd8(200)    ! ncep grid8; mercator, for gdswzd01 routine
 data grd8 / 1, 116, 44, -48670, 3104, 128, 61050,  &
             0, 22500, 0, 64, 318830, 318830, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd127(200)  ! ncep grid 127; gaussian (t254), for gdswzd04 routine
 data grd127 /4, 768, 384, 89642, 0, 128, -89642,  &
             -469, 469, 192, 0, 0, 255, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd203(200)  ! nam e-grid, for gdswzdcb routine
 data grd203 /203, 669, 1165, -7450, -144140, 136, 54000,  &
              -106000, 90, 77, 64, 0, 0, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd212(200)  ! afwa nh polar, spherical earth, for gdswzd05 routine
 data grd212 /5,2*512,-20826,145000,8,-80000,2*47625,0,  &
              9*0,255,180*0/

 integer :: grd222(200)  ! afwa nh polar, oblate spheroid earth for gdswzd05
                         ! routine.  
 data grd222 /5,2*512,-20826,145000,72,-80000,2*47625,0,  &
              9*0,255,180*0/

 integer :: grd213(200)  ! afwa sh polar, spherical earth for gdswzd05 routine.
 data grd213/5,2*512,20826,-125000,8,-80000,2*47625,128, &
             9*0,255,180*0/

 integer :: grd205(200)  ! nam 12km b-grid, for gdswzdcd routine
 data grd205 /205, 954, 835, -7491, -144134, 136, 54000,  &
             -106000, 126, 108, 64, 44540, 14800, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd218(200)  ! lambert conformal (ncep grid 218) for gdswzd03 routine
 data grd218 /3, 614, 428, 12190, -133459, 8, -95000,  &
              12191, 12191, 0, 64, 25000, 25000, 0, 0, 0, 0, 0, 0, 255, 180*0/

 i1=1
 call getarg(i1,grid)

 kgds=0
 select case (trim(grid))
   case ('3')
     kgds=grd3
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('8')
     kgds=grd8
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('203')
     kgds=grd203
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('127')
     kgds=grd127
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('205')
     kgds=grd205
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('212')
     kgds=grd212
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('222')
     kgds=grd222
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('213')
     kgds=grd213
     imdl=kgds(2)
     jmdl=kgds(3)
   case ('218')
     kgds=grd218
     imdl=kgds(2)
     jmdl=kgds(3)
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

 call gdswzd(kgds, iopt, npts, fill, xpts, ypts, rlon, rlat, &
             nret, crot, srot, xlon, xlat, ylon, ylat, area)

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

 call gdswzd(kgds, iopt, npts, fill, xpts, ypts, rlon, rlat,&
             nret, crot, srot, xlon, xlat, ylon, ylat, area)

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

 if (kgds(1) == 203) then
   kscan=mod(kgds(11)/256,2)
   if(kscan.eq.0) THEN
     is1=(jmdl+1)/2
   else
     is1=jmdl/2
   endif
   nm=imdl*jmdl
   nscan=mod(kgds(11)/32,2)
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

 deallocate (xpts,ypts,rlat,rlon,crot,srot,xlon,xlat,ylon,ylat,area)

 98 continue

 print*,'NORMAL TERMINATION'

 stop

 55 continue
 print*,'ERROR WRITING OUTPUT FILE.'
 stop 44

 end program gdswzd_driver
