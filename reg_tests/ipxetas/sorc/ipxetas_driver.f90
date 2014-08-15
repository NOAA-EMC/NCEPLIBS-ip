 program ipxetas_driver     

!------------------------------------------------------------------------
! Test routine ipxetas, which transforms between staggered eta (rotated
! lat/lon) grids.
!
! Reads an input file of vegetation greenness on a 'filled'
! 12km eta grid, then calls routine ipxetas to do the following:
!
! 1) create a staggered mass grid from the full grid.
! 2) create a staggered velocity grid from the full grid.
! 3) create a full grid from the staggered mass grid created by step (1)
! 4) create a full grid from the staggered vel grid created by step (2)
!
! The output from the first two steps is written to binary file,
! "staggered.bin".  The output from the last two steps is written to
! binary file, "full.bin".
!
! Note: I don't think this routine works for full grids with an even
! number of points in the x-direction.  So, the input test 
! data (full grid) has an odd number of points.
!------------------------------------------------------------------------

 implicit none

 character*100                  :: file_full

 integer                        :: idir, km, m1_m, m1_v, m2, i, j, nn, iend
 integer                        :: kgds_full(200), kgds_stag(200)
 integer                        :: i_full, j_full, i_stag, j_stag
 integer                        :: iret, lskip, jpds(200), jgds(200), kpds(200)
 integer                        :: lugi, numbytes, numpts, message_num
 integer, parameter             :: iunit=9

 logical*1, allocatable         :: bitmap_full(:,:)

 real, allocatable              :: data_full(:,:), data_stag_m(:), data_stag_v(:)
 real, allocatable              :: data_full_m(:,:), data_full_v(:,:)
 real, allocatable              :: data_stag_m_2d(:,:), data_stag_v_2d(:,:)

!-------------------------------------------------------------------------------
! read in eta grid of vegetation greenness.  data is on the full 
! (not staggered) grid.
!-------------------------------------------------------------------------------

 file_full="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(file_full)
 call baopenr (iunit, file_full, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 lugi       = 0
 lskip      = -1
 jpds       = -1
 jgds       = -1
 jpds(5)    = 255   
 kpds       = jpds
 kgds_full  = jgds

 print*,"- GET GRIB HEADER"
 call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
             numpts, message_num, kpds, kgds_full, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF HEADER. IRET IS ", iret
   stop 3
 else
   print*,'- SUCCESSFULL DEGRIB OF HEADER'
 end if

 i_full = kgds_full(7)
 j_full = kgds_full(8)

 allocate(data_full(i_full,j_full))
 allocate(bitmap_full(i_full,j_full))

 print*,"- DEGRIB DATA"
 call getgb(iunit, lugi, numpts, lskip, jpds, jgds, &
            numpts, lskip, kpds, kgds_full, bitmap_full, data_full, iret)

 deallocate (bitmap_full)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

!-------------------------------------------------------------------------------
! call ipxetas to created staggered mass grid from the full grid.
!-------------------------------------------------------------------------------

 km   = 1 ! interpolate one field 
 m2   = i_full * j_full

 m1_m = (i_full * j_full + 1) / 2   ! number of mass points

 idir = -1  ! calculate mass grid

 allocate (data_stag_m(m1_m))

 call ipxetas(idir, m1_m, m2, km, kgds_stag, data_stag_m, kgds_full, data_full, iret)

 if (iret /= 0) then
   print*,"- BAD CALL TO IPXETAS. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL CALL TO IPXETAS FOR FULL TO MASS GRID."
 end if

 print*,'- KGDS_FULL ',kgds_full(1:20)
 print*,'- KGDS_STAG ',kgds_stag(1:20)

!-------------------------------------------------------------------------------
! call ipxetas to create staggered velocity grid.
!-------------------------------------------------------------------------------

 idir = -2  ! calculate velocity grid from full grid.

 m1_v = (i_full * j_full + 1) / 2 - 1  ! number of velocity points.

 allocate (data_stag_v(m1_v))

 call ipxetas(idir, m1_v, m2, km, kgds_stag, data_stag_v, kgds_full, data_full, iret)

 deallocate(data_full)

 if (iret /= 0) then
   print*,"- BAD CALL TO IPXETAS. IRET IS ", iret
   stop 46
 else
   print*,"- SUCCESSFULL CALL TO IPXETAS FOR FULL TO VELOCITY GRID."
 end if

 print*,'- KGDS_FULL ',kgds_full(1:20)
 print*,'- KGDS_STAG ',kgds_stag(1:20)

!-------------------------------------------------------------------------------
! convert staggered arrays to 2d for display in grads.
!-------------------------------------------------------------------------------

 i_stag = kgds_stag(7)
 j_stag = kgds_stag(8)

 allocate (data_stag_m_2d(i_stag,j_stag))
 data_stag_m_2d = -9999.

 nn = 0
 do j = 1, j_stag
 do i = 1, i_stag
   if (mod(j,2) == 0 .and. i == i_stag) cycle
   nn = nn + 1
   data_stag_m_2d(i,j) = data_stag_m(nn)
 enddo
 enddo

 allocate (data_stag_v_2d(i_stag,j_stag))
 data_stag_v_2d = -9999.

 nn = 0
 do j = 1, j_stag
 do i = 1, i_stag
   if (mod(j,2) /= 0 .and. i == i_stag) cycle
   nn = nn + 1
   data_stag_v_2d(i,j) = data_stag_v(nn)
 enddo
 enddo

 open (33, file="./staggered.bin", access='direct', recl=i_stag*j_stag*4, err=77)
 write (33, rec=1, err=77) real(data_stag_m_2d,4)
 write (33, rec=2, err=77) real(data_stag_v_2d,4)
 close (33)

 deallocate (data_stag_m_2d, data_stag_v_2d)

!-------------------------------------------------------------------------------
! call ipxetas to create full field from staggered mass field.
!-------------------------------------------------------------------------------

 idir = 1
 allocate (data_full_m(i_full,j_full))
 data_full_m = -9999.

 call ipxetas(idir, m1_m, m2, km, kgds_stag, data_stag_m, kgds_full, data_full_m, iret)

 deallocate (data_stag_m)

 if (iret /= 0) then
   print*,"- BAD CALL TO IPXETAS. IRET IS ", iret
   stop 47
 else
   print*,"- SUCCESSFULL CALL TO IPXETAS FOR MASS TO FULL GRID."
 end if

 print*,'- KGDS_FULL ',kgds_full(1:20)
 print*,'- KGDS_STAG ',kgds_stag(1:20)

!-------------------------------------------------------------------------------
! call ipxetas to create full field from staggered velocity field.
!-------------------------------------------------------------------------------

 idir = 2
 allocate (data_full_v(i_full,j_full))
 data_full_v = -9999.

 call ipxetas(idir, m1_v, m2, km, kgds_stag, data_stag_v, kgds_full, data_full_v, iret)

 deallocate (data_stag_v)

 if (iret /= 0) then
   print*,"- BAD CALL TO IPXETAS. IRET IS ", iret
   stop 48
 else
   print*,"- SUCCESSFULL CALL TO IPXETAS FOR VELOCITY TO FULL GRID."
 end if

 print*,'- KGDS_FULL ',kgds_full(1:20)
 print*,'- KGDS_STAG ',kgds_stag(1:20)

 open (38, file="./full.bin", access='direct', recl=i_full*j_full*4, err=77)
 write (38, rec=1, err=77) real(data_full_m,4)
 write (38, rec=2, err=77) real(data_full_v,4)
 close (38)

 deallocate (data_full_m, data_full_v)

 print*,'- NORMAL TERMINATION'

 stop 

 77 continue
 print*,'- ERROR WRITING BINARY FILE'
 stop 88

 end program ipxetas_driver
