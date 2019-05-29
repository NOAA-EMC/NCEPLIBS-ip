 subroutine interp

!-----------------------------------------------------------------------
! call polateg series of routines to determine the vector gradients
! of a field of scalars interpolated to a high-res lat/lon grid.
!
! polateg0 - uses bilinear interpolation
! polateg1 - uses bicubic interpolation
! polateg4 - uses spectral interpolation
!
! the computed data from each routine call is placed in its own
! direct access file with the following names:
! "polateg0.bin", "polateg1.bin", "polateg4.bin"
!-----------------------------------------------------------------------

 use get_input_data, only : input_kgds, input_data, input_bitmap

 implicit none

 integer     :: i, j, iret, km, ibi, mi, mo, no, ibo
 integer     :: i_input, j_input
 integer     :: i_output, j_output
 integer     :: ipopt(20)  ! not used by polateg0
 integer     :: output_kgds(200)

 logical*1, allocatable :: output_bitmap(:,:)

 real, allocatable :: output_rlat(:,:), output_rlon(:,:)
 real, allocatable :: output_crot(:,:), output_srot(:,:)
 real, allocatable :: output_xo(:,:), output_yo(:,:)

! note: the spectal routine - polateg4 - requires a global grid with
! the x-axis starting at greenwich.

 integer :: grd0(200)    ! a global lat/lon grid.  grid to interpolate to.
 data grd0 / 0, 1440, 720, 89875, 0, 128, -89875, 359750,  &
             250, 250, 0, 0, -1, -1, -1, -1, -1, -1, 0, 255, 180*0/

 i_input = input_kgds(2)
 j_input = input_kgds(3)

 km  = 1   ! # of fields to interpolate
 ibi = 1   ! input data bitmap flag
 mi  = i_input * j_input

 output_kgds = grd0

 i_output = output_kgds(2)
 j_output = output_kgds(3)

 mo = i_output*j_output

 allocate (output_rlat(i_output,j_output))
 allocate (output_rlon(i_output,j_output))
 allocate (output_srot(i_output,j_output))
 allocate (output_crot(i_output,j_output))
 allocate (output_bitmap(i_output,j_output))
 allocate (output_xo(i_output,j_output))
 allocate (output_yo(i_output,j_output))

 print*,'- CALL POLATEG0'

 call polateg0(ipopt, input_kgds, output_kgds, mi, mo, km, ibi, &
               input_bitmap, input_data, no, &
               output_rlat, output_rlon, output_crot, output_srot, &
               ibo, output_bitmap, output_xo, output_yo, iret)

 if (iret /= 0) then
   print*,'- BAD CALL TO POLATEG0, IRET: ',iret
   stop 9
 endif

 do j = 1, j_output
 do i = 1, i_output
   if (.not.output_bitmap(i,j)) then
     output_xo(i,j) = -99999.
     output_yo(i,j) = -99999.
   endif
 enddo
 enddo

 open (19, file="./polateg0.bin", access='direct', recl=i_output*j_output*4)
 write(19, rec=1) real(output_rlat,4)
 write(19, rec=2) real(output_rlon,4)
 write(19, rec=3) real(output_srot,4)
 write(19, rec=4) real(output_crot,4)
 write(19, rec=5) real(output_xo,4)
 write(19, rec=6) real(output_yo,4)
 close (19)

 ibi  =0      ! spectral version does not work with bitmaps.
 ipopt=0
 ipopt(1)=0   ! triangular
 ipopt(2)=-1  ! code picks truncation number

 output_rlat=-9999.
 output_rlon=-9999.
 output_crot=-9999.
 output_srot=-9999.
 output_xo=-9999.
 output_yo=-9999.
 output_bitmap=.false.

 print*,'- CALL POLATEG4'

 call polateg4(ipopt, input_kgds, output_kgds, mi, mo, km, ibi, &
               input_bitmap, input_data, no, &
               output_rlat, output_rlon, output_crot, output_srot, &
               ibo, output_bitmap, output_xo, output_yo, iret)

 if (iret /= 0) then
   print*,'- BAD CALL TO POLATEG4, IRET: ',iret
   stop 10
 endif

 open (19, file="./polateg4.bin", access='direct', recl=i_output*j_output*4)
 write(19, rec=1) real(output_rlat,4)
 write(19, rec=2) real(output_rlon,4)
 write(19, rec=3) real(output_srot,4)
 write(19, rec=4) real(output_crot,4)
 write(19, rec=5) real(output_xo,4)
 write(19, rec=6) real(output_yo,4)
 close (19)

 ibi  =1      
 ipopt=0
 ipopt(1)=0   ! straight bicubic

 output_rlat=-9999.
 output_rlon=-9999.
 output_crot=-9999.
 output_srot=-9999.
 output_xo=-9999.
 output_yo=-9999.
 output_bitmap=.false.

 print*,'- CALL POLATEG1'

 call polateg1(ipopt, input_kgds, output_kgds, mi, mo, km, ibi, &
               input_bitmap, input_data, no, &
               output_rlat, output_rlon, output_crot, output_srot, &
               ibo, output_bitmap, output_xo, output_yo, iret)

 if (iret /= 0) then
   print*,'- BAD CALL TO POLATEG1, IRET: ',iret
   stop 12
 endif

 deallocate (input_bitmap, input_data)

 do j = 1, j_output
 do i = 1, i_output
   if (.not.output_bitmap(i,j)) then
     output_xo(i,j) = -99999.
     output_yo(i,j) = -99999.
   endif
 enddo
 enddo

 open (19, file="./polateg1.bin", access='direct', recl=i_output*j_output*4)
 write(19, rec=1) real(output_rlat,4)
 write(19, rec=2) real(output_rlon,4)
 write(19, rec=3) real(output_srot,4)
 write(19, rec=4) real(output_crot,4)
 write(19, rec=5) real(output_xo,4)
 write(19, rec=6) real(output_yo,4)
 close(19)

 deallocate (output_rlat)
 deallocate (output_rlon)
 deallocate (output_srot)
 deallocate (output_crot)
 deallocate (output_bitmap)
 deallocate (output_xo)
 deallocate (output_yo)

 return
 end subroutine interp
