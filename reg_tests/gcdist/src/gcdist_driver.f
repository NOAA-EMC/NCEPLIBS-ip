 program gcdist_driver

!---------------------------------------------------------------
! test iplib routines gcdist and movect.  both perform
! calculations related to great circle routes (distance
! and rotational parameters.  output is piped to 
! standard output.
!---------------------------------------------------------------

 implicit none

 real, external :: gcdist
 real           :: rlat1, rlat2, rlon1, rlon2, crot, srot, dist

 rlat1 = 65.0
 rlon1 = 65.0
 rlat2 = 0.0
 rlon2 = 170.0

 dist=gcdist(rlat1,rlon1,rlat2,rlon2)
 call movect(rlat1,rlon1,rlat2,rlon2,crot,srot)

 print*,'TEST 1'
 print*,'BEGIN LAT/LON: ',rlat1, rlon1
 print*,'END LAT/LON:   ',rlat2, rlon2
 print*,'GREAT CIRCLE DISTANCE: ',dist
 print*,'ROTATION COS/SIN: ', crot, srot

 rlat1 = -30.0
 rlon1 = -10.0
 rlat2 = 30.0
 rlon2 = 15.0

 dist=gcdist(rlat1,rlon1,rlat2,rlon2)
 call movect(rlat1,rlon1,rlat2,rlon2,crot,srot)

 print*,''
 print*,'TEST 2'
 print*,'BEGIN LAT/LON: ',rlat1, rlon1
 print*,'END LAT/LON:   ',rlat2, rlon2
 print*,'GREAT CIRCLE DISTANCE: ',dist
 print*,'ROTATION COS/SIN: ', crot, srot

! test special case where arc is zero length

 rlat1 = 0.0
 rlon1 = 0.0
 rlat2 = 0.0
 rlon2 = 0.0

 dist=gcdist(rlat1,rlon1,rlat2,rlon2)
 call movect(rlat1,rlon1,rlat2,rlon2,crot,srot)

 print*,''
 print*,'TEST 3'
 print*,'BEGIN LAT/LON: ',rlat1, rlon1
 print*,'END LAT/LON:   ',rlat2, rlon2
 print*,'GREAT CIRCLE DISTANCE: ',dist
 print*,'ROTATION COS/SIN: ', crot, srot

 print*,''
 print*,'NORMAL TERMINATION'

 end program gcdist_driver
