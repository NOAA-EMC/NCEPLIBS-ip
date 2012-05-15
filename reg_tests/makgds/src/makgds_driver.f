 program makgds_driver

!--------------------------------------------------------------
! test iplib routine makgds.
!--------------------------------------------------------------

 implicit none

 character*1    :: gds(400)

 integer        :: i, iopt, iret, lengds, kgds(200)

 print*,"MAKE GDS AND KGDS FOR NCEP GRID 3."

 iopt = 3  ! ncep grid 3
 kgds = -999
 gds  = repeat('X',400)
 lengds = -999
 iret   = -999

 call makgds(iopt, kgds, gds, lengds, iret)

 if (iret /= 0) then
   print*,"** ERROR IN MAKGDS **"
 endif

 print*,'LENGDS: ',lengds
 print*,'KGDS:   ',kgds(1:25)

! the gds array can have some weird looking characters that won't print
! to the screen. so convert to integers.

 do i = 1, lengds
   print*,'GDS:    ',i, ichar(gds(i:i))
 enddo

 print*,''
 print*,"MAKE KGDS FROM GDS."

 iopt   = -1
 kgds   = -999
 iret   = -999
 lengds = -999

 call makgds(iopt, kgds, gds, lengds, iret)

 if (iret /= 0) then
   print*,"** ERROR IN MAKGDS **"
 endif

 print*,'KGDS:   ',kgds(1:25)

 print*,''
 print*,"MAKE GDS FROM KGDS"

 iopt = 255
 gds = repeat('X',400)
 iret = -999
 lengds = -999

 call makgds(iopt, kgds, gds, lengds, iret)

 if (iret /= 0) then
   print*,"** ERROR IN MAKGDS **"
 endif

 do i = 1, lengds
   print*,'GDS:    ',i, ichar(gds(i:i))
 enddo

 print*,''
 print*,"NORMAL TERMINATION"

 end program makgds_driver
