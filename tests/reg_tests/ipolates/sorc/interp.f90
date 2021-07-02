 subroutine interp

!-------------------------------------------------------------------------
! Call the scalar polates routines to interpolate the input data
! using all available interpolation methods (neighbor, bilinear, etc.).
! Several output grids of various map projections are tested.
!
! The routine reads in two arguments from stnd input.  The first is
! the grid to which you want to interpolate the data.
! The valid choices are:
!
!    3 -  one-degree global lat/lon (ncep grid 3)
!    8 -  mercator (ncep grid 8)
!  127 -  t254 gaussian (ncep grid 127)
!  203 -  rotated lat/lon e-staggered (number meaningless).
!         this is the old 12km eta grid - 'h' pts
!  205 -  rotated lat/lon b-staggered (number meaningless).
!         this is the 12km nam grid - 'h' pts
!  212 -  nh polar stereographic, spherical earth (number meaningless)
!  218 -  lambert conformal (ncep grid 218)
!
! The second argument is the interpolation option.  The valid choices:
!
!  0 - bilinear
!  1 - bicubic
!  2 - neighbor
!  3 - budget
!  4 - spectral
!  6 - budget-neighbor
!
! The interpolated data is output in a binary file with the following
! naming convention:  "grid${grid_num}.opt${interp_opt}.bin"
!-------------------------------------------------------------------------

 use get_input_data, only : input_data, &
                            gdtmpl_input, &
                            gdtlen_input, &
                            gdtnum_input, & 
                            input_bitmap, &
                            i_input, j_input

 implicit none

 character*1                  :: interp_opt
 character*3                  :: grid
 character*100                :: output_file

 integer(kind=4)              :: i1
 integer                      :: ip, ipopt(20)
 integer                      :: km, ibi, mi, iret, i, j
 integer                      :: i_output, j_output, mo, no, ibo
 integer        , allocatable :: gdtmpl_output(:)
 integer                      :: gdtlen_output, gdtnum_output

 logical*1, allocatable       :: output_bitmap(:,:)

 real, allocatable            :: output_rlat(:,:), output_rlon(:,:)
 real, allocatable            :: output_data(:,:)

 integer        , parameter   :: missing=b'11111111111111111111111111111111'

 integer, parameter :: gdtlen3 = 19 ! one-degree lat/lon
 integer            :: gdtmpl3(gdtlen3)
 data gdtmpl3 / 6, 255, missing, 255, missing, 255, missing, 360, 181, 0, missing, &
                90000000, 0, 56, -90000000, 359000000, 1000000, 1000000, 0 /

 integer, parameter :: gdtlen8 = 19 ! ncep grid8; mercator
 integer            :: gdtmpl8(gdtlen8)
 data gdtmpl8 / 6, 255, missing, 255, missing, 255, missing, 116, 44, &
                -48670000, 3104000, 56, 22500000, 61050000, 0, 64, 0, &
                 318830000, 318830000/

 integer, parameter :: gdtlen127=19  ! t254 gaussain
 integer            :: gdtmpl127(gdtlen127)
 data gdtmpl127 /6, 255, missing, 255, missing, 255, missing, 768, 384, &
                 0, missing, 89642000, 0, 48, -89642000, 359531000,  &
                 469000, 192, 0/

 integer, parameter :: gdtlen203=22  ! 12km eta, 'h' pts
 integer            :: gdtmpl203(gdtlen203)
 data gdtmpl203/6, 255, missing, 255, missing, 255, missing, 669, 1165, &
                0, missing, -45000000, 300000000, 56, 45000000, 60000000, &
                179641, 77320, 68, -36000000, 254000000, 0 /

 integer, parameter :: gdtlen205=22  ! 12km nam, 'h' pts
 integer            :: gdtmpl205(gdtlen205)
 data gdtmpl205/6, 255, missing, 255, missing, 255, missing, 954, 835, &
                0, missing, -45036000, 299961000, 56, 45036000, 60039000, &
                126000, 108000, 64, -36000000, 254000000, 0 /

 integer, parameter :: gdtlen212=18  ! nh polar
 integer            :: gdtmpl212(gdtlen212)
 data gdtmpl212 /6, 255, missing, 255, missing, 255, missing, 513, 513, &
                 -20826000, 235000000, 56, 60000000, 280000000, 47625000, 47625000, &
                 0, 64/

 integer, parameter :: gdtlen218 = 22 ! ncep grid 218; lambert conf
 integer            :: gdtmpl218(gdtlen218)
 data gdtmpl218 / 6, 255, missing, 255, missing, 255, missing, 614, 428, &
                  12190000, 226541000, 56, 25000000, 265000000, &
                  12191000, 12191000, 0, 64, 25000000, 25000000, -90000000, 0/

 i1=1
 call getarg(i1, grid)
 i1=2
 call getarg(i1, interp_opt)

 select case (trim(grid))
   case ('3')
     gdtnum_output=0
     gdtlen_output=gdtlen3
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl3
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case ('8')
     gdtnum_output=10
     gdtlen_output=gdtlen8
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl8
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case ('127')
     gdtnum_output=40
     gdtlen_output=gdtlen127
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl127
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case ('203')
     gdtnum_output=1
     gdtlen_output=gdtlen203
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl203
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case ('205')
     gdtnum_output=1
     gdtlen_output=gdtlen205
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl205
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case ('212')
     gdtnum_output=20
     gdtlen_output=gdtlen212
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl212
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case ('218')
     gdtnum_output=30
     gdtlen_output=gdtlen218
     allocate(gdtmpl_output(gdtlen_output))
     gdtmpl_output=gdtmpl218
     i_output = gdtmpl_output(8)
     j_output = gdtmpl_output(9)
   case default
     print*,"ERROR: ENTER VALID GRID NUMBER."
     stop 55
 end select
 print*,"- CALL IPOLATES FOR GRID: ", grid

 select case (interp_opt)
! bi-linear
   case ('0')           
     ip = 0
     ipopt = 0  
     ibi   = 1  ! use bitmap
! bi-cubic
   case ('1')          
     ip = 1
     ipopt = 0  
     ibi   = 1  ! use bitmap
! neighbor
   case ('2')       
     ip = 2
     ipopt = 0  
     ibi   = 1  ! use bitmap
! budget
   case ('3')       
     ip = 3
     ipopt = -1  
     ibi   = 1  ! use bitmap
! spectral 
   case ('4')         
     ip = 4
     ipopt(1) = 0  ! triangular
     ipopt(2) = -1 ! code chooses wave number
     ipopt(3:20)=0
     ibi   = 0     ! can't use bitmap with spectral
! neighbor-budget
   case ('6')            
     ip = 6
     ipopt = -1  
     ibi   = 1  ! use bitmap
   case default
     print*,"ERROR: ENTER VALID INTERP OPTION."
     stop 56
 end select
 print*,"- USE INTERP OPTION: ", interp_opt

 km = 1  ! number of fields to interpolate
 mi = i_input * j_input ! dimension of input grids

 mo = i_output * j_output 

 allocate (output_rlat(i_output,j_output))
 allocate (output_rlon(i_output,j_output))
 allocate (output_data(i_output,j_output))
 allocate (output_bitmap(i_output,j_output))

 print*,'- CALL IPOLATES'
 call ipolates(ip, ipopt, gdtnum_input, gdtmpl_input, gdtlen_input, & 
               gdtnum_output, gdtmpl_output, gdtlen_output, &
               mi, mo, km, ibi, input_bitmap, input_data, &
               no, output_rlat, output_rlon, ibo, output_bitmap, output_data, iret)

 deallocate(input_bitmap, input_data)

 if (iret /= 0) then
   print*,'- BAD STATUS FROM IPOLATES: ', iret
   stop 6
 end if

 if (no /= mo) then
   print*,'- ERROR: WRONG # OF POINTS RETURNED FROM IPOLATES'
   stop 7
 end if

 print*,'- SUCCESSFULL CALL TO IPOLATES'

 do j = 1, j_output
 do i = 1, i_output
    if (.not. output_bitmap(i,j)) then
      output_data(i,j) = -9.
    endif
 enddo
 enddo

 output_file = "./grid" // trim(grid) // ".opt" // interp_opt // ".bin" 
 open (12, file=output_file, access="direct", err=38, recl=mo*4)
 write (12, err=38, rec=1) real(output_data,4)
 close (12)

 deallocate (output_rlat, output_rlon, output_data, output_bitmap)

 return

 38 continue

 print*,'-ERROR WRITING BINARY FILE.'
 stop 77

 end subroutine interp
