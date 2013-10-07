 program ipsector_driver     

!----------------------------------------------------------------------------
! Read in a global dataset of substrate temperature,
! then call ipsector to create a subsector of
! the original grid.  Then call ipspaste to 'paste'
! the subsectored data (created by call to ipsector) back
! to the original grid.  The data returned from ipsector
! should match the original data.
!
! Three sets of calls to ipsector/ipspaste are made:
! - for a North America subsector
! - for a non-overlapping subsector
! - for an overlapping subsector
!
! The output from the three calls to ipsector are stored in
! files "ipsector.namer.bin", "ipsector.no.sect.bin" and
! "ipsector.ovlp.sect.bin". 
!
! The output from the three calls to ipspaste are stored in
! files "ipspaste.namer.bin", "ipspaste.no.sect.bin" and
! "ipspaste.ovlp.sect.bin". 
!----------------------------------------------------------------------------

 implicit none

 character*100                  :: orig_file

 integer                        :: i, j
 integer                        :: i1, i2, j1, j2, m, ms, nf
 integer                        :: i1a, i2a, j1a, j2a
 integer                        :: orig_kgds(200), paste_kgds(200), sector_kgds(200)
 integer                        :: i_orig, j_orig, i_paste, j_paste
 integer                        :: i_sector, j_sector
 integer                        :: iret, lskip, jpds(200), jgds(200), kpds(200)
 integer                        :: lugi, numbytes, numpts, message_num
 integer, parameter             :: iunit=9

 logical*1, allocatable         :: orig_bitmap(:,:)
 logical*1, allocatable         :: paste_bitmap(:,:)
 logical*1, allocatable         :: sector_bitmap(:,:)

 real, allocatable              :: orig_data(:,:)
 real, allocatable              :: paste_data(:,:)
 real, allocatable              :: sector_data(:,:)

!-------------------------------------------------------------------------------
! read in global grid of substrate temperatures.
!-------------------------------------------------------------------------------

 orig_file="./fort.9"
 print*,"- OPEN AND READ INPUT DATA: ", trim(orig_file)
 call baopenr (iunit, orig_file, iret)

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
 orig_kgds  = jgds

 print*,"- GET GRIB HEADER"
 call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
             numpts, message_num, kpds, orig_kgds, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF HEADER. IRET IS ", iret
   stop 3
 else
   print*,'- SUCCESSFULL DEGRIB OF HEADER'
 end if

 i_orig = orig_kgds(2)
 j_orig = orig_kgds(3)

 allocate(orig_data(i_orig,j_orig))
 allocate(orig_bitmap(i_orig,j_orig))

 print*,"- DEGRIB DATA"
 call getgb(iunit, lugi, numpts, lskip, jpds, jgds, &
            numpts, lskip, kpds, orig_kgds, orig_bitmap, orig_data, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

 open (24, file="./orig.bin", access='direct', err=77, recl=i_orig*j_orig*4)
 write (24, rec=1, err=77) real(orig_data,4)
 close (24)

!--------------------------------------------------------------------------------
! call ipsector to chop out north america.
!--------------------------------------------------------------------------------

 print*,'- CREATE NORTH AMERICA SECTOR.'

 i1 = 1          ! first i point of sector
 i2 = i_orig/2   ! last i point of sector
 j1 = 1          ! first j point of sector
 j2 = j_orig/2   ! last j point of sector

 i_sector = i2 - i1 + 1
 j_sector = j2 - j1 + 1
 ms       = i_sector * j_sector

 nf = 1  ! # of fields to cut

 m = i_orig * j_orig

 allocate (sector_bitmap(i_sector,j_sector))
 allocate (sector_data(i_sector,j_sector))

 call ipsector(i1, i2, j1, j2, nf, m, orig_kgds, orig_bitmap, &
               orig_data, ms, sector_kgds, sector_bitmap, sector_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSECTOR, IRET IS ", iret
   stop 44
 else
   print*,"- SUCCESSFULL CALL TO IPSECTOR."
 end if

 print*,"- ORIGINAL GRID KGDS: ",orig_kgds(1:20)
 print*,"- IPSECTOR GRID KGDS: ",sector_kgds(1:20)

 open (11, file="./ipsector.namer.bin", access='direct', err=77, recl=ms*4)
 write (11, rec=1, err=77) real(sector_data,4)
 close (11)

!--------------------------------------------------------------------------------
! now paste back north america.  the result should match the original data.
!--------------------------------------------------------------------------------

 nf = 1
 ms = i_sector * j_sector

 i_paste = i_orig
 j_paste = j_orig

 m = i_paste * j_paste

 paste_kgds = orig_kgds

 allocate(paste_bitmap(i_paste,j_paste))
 allocate(paste_data(i_paste,j_paste))

 paste_data = orig_data
 paste_bitmap = orig_bitmap

 call ipspaste(i1, i2, j1, j2, nf, ms, sector_bitmap, sector_data, &
               m, paste_kgds, paste_bitmap, paste_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSPASTE, IRET IS ", iret
   stop 47
 else
   print*,"- SUCCESSFULL CALL TO IPSPASTE."
 end if

 open (21, file="./ipspaste.namer.bin", access='direct', recl=m*4, err=77)
 write (21, rec=1, err=77) real(paste_data,4)   ! data returned from ipspaste
 close (21)

!--------------------------------------------------------------------------------
! test non-overlapping sectors logic.  first, call ipsector to create
! a subsector of data.
!--------------------------------------------------------------------------------

 print*,'- CREATE NON-OVERLAPPING SECTOR.'

 deallocate (sector_bitmap, sector_data)

 i1 = 4    ! number of sectors
 i2 = 3    ! sector number; goes from 0 to (#sectors - 1).
           ! i.e., uses "C" langauge indexing.

 I1A=MIN(I2*((orig_kgds(2)-1)/I1+1)+1,orig_kgds(2)+1)
 I2A=MIN((I2+1)*((orig_kgds(2)-1)/I1+1),orig_kgds(2))

 print*,"- NUMBER OF SUBSECTORS IN I DIRECTION: ",i1
 print*,"- SECTOR NUMBER IN I DIRECTION:        ",i2
 print*,"- START/END I INDICES                  ",i1a,i2a

 j1 = 5    ! number of sectors
 j2 = 0    ! sector number; goes from 0 to (#sectors - 1)

 J1A=MIN(J2*((orig_kgds(3)-1)/J1+1)+1,orig_kgds(3)+1)
 J2A=MIN((J2+1)*((orig_kgds(3)-1)/J1+1),orig_kgds(3))

 print*,"- NUMBER OF SUBSECTORS IN J DIRECTION: ",j1
 print*,"- SECTOR NUMBER IN J DIRECTION:        ",j2
 print*,"- START/END J INDICES                  ",j1a,j2a

 i_sector = i2a - i1a + 1
 j_sector = j2a - j1a + 1
 ms       = i_sector * j_sector

 nf = 1  ! # of fields to cut

 m = i_orig * j_orig

 allocate (sector_bitmap(i_sector,j_sector))
 allocate (sector_data(i_sector,j_sector))

 call ipsector(i1, i2, j1, j2, nf, m, orig_kgds, orig_bitmap, &
               orig_data, ms, sector_kgds, sector_bitmap, sector_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSECTOR, IRET IS ", iret
   stop 47
 else
   print*,"- SUCCESSFULL CALL TO IPSECTOR."
 end if

 print*,"- ORIGINAL GRID KGDS: ",orig_kgds(1:20)
 print*,"- IPSECTOR GRID KGDS: ",sector_kgds(1:20)

 open (35, file="./ipsector.no.sect.bin", access='direct', recl=ms*4, err=77)
 write (35, rec=1, err=77) real(sector_data,4)
 close (35)

!--------------------------------------------------------------------------------
! now paste back to original grid.
!--------------------------------------------------------------------------------

 i_paste = i_orig
 j_paste = j_orig

 m = i_paste * j_paste

 paste_kgds = orig_kgds
 paste_data = orig_data
 paste_bitmap = orig_bitmap

 call ipspaste(i1, i2, j1, j2, nf, ms, sector_bitmap, sector_data, &
               m, paste_kgds, paste_bitmap, paste_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSPASTE, IRET IS ", iret
   stop 47
 else
   print*,"- SUCCESSFULL CALL TO IPSPASTE."
 end if

 open (21, file="./ipspaste.no.sect.bin", access='direct', recl=m*4, err=77)
 write (21, rec=1, err=77) real(paste_data,4)   ! data returned from ipspaste
 close (21)

!--------------------------------------------------------------------------------
! test overlapping sectors logic.  first, call ipsector to create
! a subsector of of data.
!--------------------------------------------------------------------------------

 print*,'- CREATE OVERLAPPING SECTOR.'

 deallocate (sector_bitmap, sector_data)

 i1 = -4   ! number of sectors
 i2 = 3    ! sector number; goes from 0 to (#sectors - 1).
           ! i.e., uses "C" langauge indexing.

 I1A=MIN(I2*((orig_kgds(2)-2)/(-I1)+1)+1,orig_kgds(2)+1)
 I2A=MIN((I2+1)*((orig_kgds(2)-2)/(-I1)+1)+1,orig_kgds(2))

 print*,"- NUMBER OF SUBSECTORS IN I DIRECTION: ",abs(i1)
 print*,"- SECTOR NUMBER IN I DIRECTION:        ",i2
 print*,"- START/END I INDICES                  ",i1a,i2a

 j1 = -5   ! number of sectors
 j2 = 0    ! sector number; goes from 0 to (#sectors - 1)

 J1A=MIN(J2*((orig_kgds(3)-2)/(-J1)+1)+1,orig_kgds(3)+1)
 J2A=MIN((J2+1)*((orig_kgds(3)-2)/(-J1)+1)+1,orig_kgds(3))

 print*,"- NUMBER OF SUBSECTORS IN J DIRECTION: ",abs(j1)
 print*,"- SECTOR NUMBER IN J DIRECTION:        ",j2
 print*,"- START/END J INDICES                  ",j1a,j2a

 i_sector = i2a - i1a + 1
 j_sector = j2a - j1a + 1
 ms       = i_sector * j_sector

 nf = 1  ! # of fields to cut

 m = i_orig * j_orig

 allocate (sector_bitmap(i_sector,j_sector))
 allocate (sector_data(i_sector,j_sector))

 call ipsector(i1, i2, j1, j2, nf, m, orig_kgds, orig_bitmap, &
               orig_data, ms, sector_kgds, sector_bitmap, sector_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSECTOR, IRET IS ", iret
   stop 49
 else
   print*,"- SUCCESSFULL CALL TO IPSECTOR."
 end if

 print*,"- ORIGINAL GRID KGDS: ",orig_kgds(1:20)
 print*,"- IPSECTOR GRID KGDS: ",sector_kgds(1:20)

 open (35, file="./ipsector.ovlp.sect.bin", access='direct', recl=ms*4, err=77)
 write (35, rec=1, err=77) real(sector_data,4)
 close (35)

!--------------------------------------------------------------------------------
! now paste back to original grid.
!--------------------------------------------------------------------------------

 i_paste = i_orig
 j_paste = j_orig

 m = i_paste * j_paste

 paste_kgds = orig_kgds
 paste_data = orig_data
 paste_bitmap = orig_bitmap

 call ipspaste(i1, i2, j1, j2, nf, ms, sector_bitmap, sector_data, &
               m, paste_kgds, paste_bitmap, paste_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSPASTE, IRET IS ", iret
   stop 66
 else
   print*,"- SUCCESSFULL CALL TO IPSPASTE."
 end if

 open (21, file="./ipspaste.ovlp.sect.bin", access='direct', recl=m*4, err=77)
 write (21, rec=1, err=77) real(paste_data,4)   ! data returned from ipspaste
 close (21)

 deallocate (paste_bitmap, paste_data)
 deallocate (sector_bitmap, sector_data)
 deallocate (orig_bitmap, orig_data)

 print*,'- NORMAL TERMINATION'

 stop

 77 continue
 print*,'- ERROR WRITING A BINARY FILE.'
 stop 1

 end program ipsector_driver
