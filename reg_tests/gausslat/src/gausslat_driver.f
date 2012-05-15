 program gausslat_driver

!--------------------------------------------------------------
! test iplib routine gausslat. 
!
! is this routine even used anymore?
!--------------------------------------------------------------

 implicit none

 integer :: j, jmax

 real, allocatable :: slat(:), wlat(:)

 jmax = 384  ! t382 grid

 allocate (slat(jmax))
 allocate (wlat(jmax))

 print*,'CALL ROUTINE GAUSSLAT'

 call gausslat(jmax,slat,wlat)

 do j = 1, jmax
   print*,'J/SLAT/WLAT ',j, slat(j), wlat(j)
 enddo

 deallocate (slat, wlat)

 print*,'NORMAL TERMINATION'

 end program gausslat_driver
