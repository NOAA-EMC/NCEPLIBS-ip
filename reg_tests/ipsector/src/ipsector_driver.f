 program ipsector_driver     

!------------------------------------------------------------------------
! test routine ipsector, which creates a subset of a larger 
! two-dimensional field, and routine ipspaste, which does the
! opposite.
!------------------------------------------------------------------------

 implicit none

 character*100                  :: input_file

 integer                        :: i, j
 integer                        :: i1, i2, j1, j2, m, ms, nf
 integer                        :: input_kgds(200)
 integer                        :: output_kgds(200)
 integer                        :: i_input, j_input
 integer                        :: i_output, j_output
 integer                        :: iret, lskip, jpds(200), jgds(200), kpds(200)
 integer                        :: lugi, numbytes, numpts, message_num
 integer, parameter             :: iunit=9

 logical*1, allocatable         :: input_bitmap(:,:)
 logical*1, allocatable         :: output_bitmap(:,:)

 real, allocatable              :: input_data(:,:)
 real, allocatable              :: input_data_sav(:,:)
 real, allocatable              :: output_data(:,:)

!-------------------------------------------------------------------------------
! read in global grid of substrate temperatures.
!-------------------------------------------------------------------------------

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 call baopenr (iunit, input_file, iret)

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
 input_kgds = jgds

 print*,"- GET GRIB HEADER"
 call getgbh(iunit, lugi, lskip, jpds, jgds, numbytes,  &
             numpts, message_num, kpds, input_kgds, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF HEADER. IRET IS ", iret
   stop 3
 else
   print*,'- SUCCESSFULL DEGRIB OF HEADER'
 end if

 i_input = input_kgds(2)
 j_input = input_kgds(3)

 allocate(input_data(i_input,j_input))
 allocate(input_bitmap(i_input,j_input))

 print*,"- DEGRIB DATA"
 call getgb(iunit, lugi, numpts, lskip, jpds, jgds, &
            numpts, lskip, kpds, input_kgds, input_bitmap, input_data, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 allocate(input_data_sav(i_input,j_input))
 input_data_sav=input_data

 call baclose (iunit, iret)

!--------------------------------------------------------------------------------
! call ipsector to chop out north america.
!--------------------------------------------------------------------------------

 i1 = 1           ! first i point of sector
 i2 = i_input/2   ! last i point of sector
 j1 = 1           ! first j point of sector
 j2 = j_input/2   ! last j point of sector

 i_output = i2 - i1 + 1
 j_output = j2 - j1 + 1
 ms       = i_output * j_output

 nf = 1  ! # of fields to cut

 m = i_input * j_input

 allocate (output_bitmap(i_output,j_output))
 allocate (output_data(i_output,j_output))

 call ipsector(i1, i2, j1, j2, nf, m, input_kgds, input_bitmap, &
               input_data, ms, output_kgds, output_bitmap, output_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSECTOR, IRET IS ", iret
   stop 44
 else
   print*,"- SUCCESSFULL CALL TO IPSECTOR."
 end if

 print*,"- INPUT GRID KGDS:  ",input_kgds(1:20)
 print*,"- OUTPUT GRID KGDS: ",output_kgds(1:20)

 do j = 1, j_output
 do i = 1, i_output
   if (.not.output_bitmap(i,j)) then
!     output_data(i,j) = -9999.
   endif
 enddo
 enddo

 open (11, file="./ipsector.bin", access='direct', recl=ms*4)
 write (11, rec=1) real(output_data,4)
 close (11)

!--------------------------------------------------------------------------------
! now paste back north america.  the result should match the original data.
!--------------------------------------------------------------------------------

 call ipspaste(i1, i2, j1, j2, nf, ms, output_bitmap, output_data, &
               m, input_kgds, input_bitmap, input_data, iret)

 if (iret /= 0) then
   print*,"- BAD STATUS FROM IPSPASTE, IRET IS ", iret
   stop 47
 else
   print*,"- SUCCESSFULL CALL TO IPSPASTE."
 end if

 open (21, file="./ipspaste.bin", access='direct', recl=m*4)
 write (21, rec=1) real(input_data_sav,4)
 write (21, rec=2) real(input_data,4)
 close (21)

 deallocate (output_bitmap, output_data)
 deallocate (input_bitmap, input_data, input_data_sav)

 print*,'- NORMAL TERMINATION'

 stop

 end program ipsector_driver
