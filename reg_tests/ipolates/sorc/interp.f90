 subroutine interp

!-------------------------------------------------------------------------
! Call the scalar polates routines to interpolate the input data
! using all available interpolation methods (neighbor, bilinear, etc.)
! several output grids of various map projections are tested.
!
! The routine reads in two arguments from stnd input.  The first is
! the grid to which you want to interpolate the data.
! The valid choices are:
!
!    3 -  one-degree global lat/lon (ncep grid 3)
!    8 -  mercator (ncep grid 8)
!  127 -  t254 gaussian (ncep grid 127)
!  203 -  rotated lat/lon e-staggered (number refers to gds octet 6)
!  205 -  rotated lat/lon b-staggered (number refers to gds octet 6)
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

 character*1             :: interp_opt
 character*3             :: grid
 character*100           :: output_file

 integer*4   :: i1
 integer     :: ip, ipopt(20), output_kgds(200)
 integer     :: km, ibi, mi, iret, i, j, n
 integer     :: i_output, j_output, mo, no, ibo
 integer(kind=4), allocatable :: gdtmpl_output(:)
 integer     :: gdtlen_output, gdtnum_output

 logical*1, allocatable :: output_bitmap(:,:)

 real, allocatable :: output_rlat(:,:), output_rlon(:,:)
 real, allocatable :: output_data(:,:)

 integer(kind=4), parameter   :: missing=b'11111111111111111111111111111111'

 integer, parameter :: gdtlen3 = 19 ! one-degree lat/lon
 integer(kind=4)    :: gdtmpl3(gdtlen3)
 data gdtmpl3 / 6, 255, missing, 255, missing, 255, missing, 360, 181, 0, missing, &
                90000000, 0, 56, -90000000, 359000000, 1000000, 1000000, 0 /

 integer :: grd8(200)    ! mercator
 data grd8 / 1, 116, 44, -48670, 3104, 128, 61050,  &
             0, 22500, 0, 64, 318830, 318830, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd127(200)  ! gaussian (t254)
 data grd127 /4, 768, 384, 89642, 0, 128, -89642,  &
             -469, 469, 192, 0, 0, 255, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd203(200)  ! nam 12km e-grid
 data grd203 /203, 669, 1165, -7450, -144140, 136, 54000,  &
              -106000, 90, 77, 64, 0, 0, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd205(200)  ! nam 12km b-grid
 data grd205 /205, 954, 835, -7491, -144134, 136, 54000,  &
             -106000, 126, 90, 64, 44540, 14800, 0, 0, 0, 0, 0, 0, 255, 180*0/

 integer :: grd212(200)  ! afwa nh polar, spherical earth
 data grd212 /5,2*512,-20826,145000,8,-80000,2*47625,0,  &
              9*0,255,180*0/

 integer :: grd218(200)  ! lambert conformal (ncep grid 218) 
 data grd218 /3, 614, 428, 12190, -133459, 8, -95000,  &
              12191, 12191, 0, 64, 25000, 25000, 0, 0, 0, 0, 0, 0, 255, 180*0/
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
     output_kgds = grd8
     i_output = output_kgds(2)
     j_output = output_kgds(3)
   case ('127')
     output_kgds = grd127
     i_output = output_kgds(2)
     j_output = output_kgds(3)
   case ('203')
     output_kgds = grd203
     i_output = output_kgds(2)
     j_output = output_kgds(3)
   case ('205')
     output_kgds = grd205
     i_output = output_kgds(2)
     j_output = output_kgds(3)
   case ('212')
     output_kgds = grd212
     i_output = output_kgds(2)
     j_output = output_kgds(3)
   case ('218')
     output_kgds = grd218
     i_output = output_kgds(2)
     j_output = output_kgds(3)
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

 do n=1,2

 print*,'- CALL IPOLATES, ITERATION: ', n
 call ipolates(ip, ipopt, gdtnum_input, gdtmpl_input, gdtlen_input, & 
               gdtnum_output, gdtmpl_output, gdtlen_output, &
               mi, mo, km, ibi, input_bitmap, input_data, &
               no, output_rlat, output_rlon, ibo, output_bitmap, output_data, iret)

!deallocate(input_bitmap, input_data)

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

 if(n==1) output_file = "./grid" // trim(grid) // ".opt" // interp_opt // ".bin1" 
 if(n==2) output_file = "./grid" // trim(grid) // ".opt" // interp_opt // ".bin2" 
 open (12, file=output_file, access="direct", err=38, recl=mo*4)
 write (12, err=38, rec=1) real(output_data,4)
 close (12)

 output_bitmap=.false.
 output_data=-999.

 enddo

 deallocate (output_rlat, output_rlon, output_data, output_bitmap)

 return

 38 continue

 print*,'-ERROR WRITING BINARY FILE.'
 stop 77


 end subroutine interp
