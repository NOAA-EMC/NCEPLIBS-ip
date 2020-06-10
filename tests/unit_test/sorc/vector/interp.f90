 subroutine interp

!-------------------------------------------------------------------------
! Call the vector polates routines to interpolate the input data
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
! The interpolated data is compared to a baseline set of data.
! Differences are printed to standard output.
!-------------------------------------------------------------------------

 use get_input_data, only : input_u_data, input_v_data, &
                            input_kgds, &
                            input_bitmap, &
                            i_input, j_input

 implicit none

 character*1               :: interp_opt
 character*3               :: grid
 character*100             :: output_file

 integer(kind=4)           :: i1
 integer                   :: ip, ipopt(20), output_kgds(200)
 integer                   :: km, ibi, mi, iret, i, j
 integer                   :: i_output, j_output, mo, no, ibo
 integer                   :: num_upts_diff, num_vpts_diff

 logical*1, allocatable    :: output_bitmap(:,:)

 real, allocatable         :: output_rlat(:,:), output_rlon(:,:)
 real, allocatable         :: output_crot(:,:), output_srot(:,:)
 real, allocatable         :: output_u_data(:,:), output_v_data(:,:)
 real                      :: avg_u_diff, avg_v_diff
 real                      :: max_u_diff, max_v_diff
 real(kind=4)              :: output_data4
 real(kind=4), allocatable :: baseline_u_data(:,:)
 real(kind=4), allocatable :: baseline_v_data(:,:)

 integer :: grd3(200)    ! global one-degree lat/lon
 data grd3 / 0, 360, 181, 90000, 0, 128, -90000,  &
            -1000, 1000, 1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 255, 180*0/

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
     output_kgds = grd3
     i_output = output_kgds(2)
     j_output = output_kgds(3)
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
     ibi   = 0  ! no bitmap
! bi-cubic
   case ('1')          
     ip = 1
     ipopt = 0  
     ibi   = 0  ! no bitmap
! neighbor
   case ('2')       
     ip = 2
     ipopt = 0  
     ibi   = 0  ! no bitmap
! budget
   case ('3')       
     ip = 3
     ipopt = -1  
     ibi   = 0  ! no bitmap
! spectral 
   case ('4')         
     ip = 4
     ipopt(1) = 0  ! triangular
     ipopt(2) = -1 ! code chooses wave number
     ipopt(3:20)=0
     ibi   = 0     ! no bitmap
! neighbor-budget
   case ('6')            
     ip = 6
     ipopt = -1  
     ibi   = 0  ! no bitmap
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
 allocate (output_u_data(i_output,j_output))
 allocate (output_v_data(i_output,j_output))
 allocate (output_bitmap(i_output,j_output))
 allocate (output_srot(i_output,j_output))
 allocate (output_crot(i_output,j_output))

 call ipolatev(ip, ipopt, input_kgds, output_kgds, mi, mo,   &
               km, ibi, input_bitmap, input_u_data, input_v_data,  &
               no, output_rlat, output_rlon, output_crot, output_srot, &
               ibo, output_bitmap, output_u_data, output_v_data, iret)

 deallocate(input_bitmap, input_u_data, input_v_data)

 if (iret /= 0) then
   print*,'- BAD STATUS FROM IPOLATES: ', iret
   stop 6
 end if

 if (no /= mo) then
   print*,'- ERROR: WRONG # OF POINTS RETURNED FROM IPOLATES'
   stop 7
 end if

!print*,'- SUCCESSFULL CALL TO IPOLATES'

!-------------------------------------------------------------------------
! polatev4 does not always initialize the ibo or output_bitmap variables,
! so they can't be used.  this should be fixed.
! according to the comments, polatev4 only outputs a bitmap in areas
! that extend outside the domain of the input grid. 
!-------------------------------------------------------------------------
 
 if (ip /= 4) then
   do j = 1, j_output
   do i = 1, i_output
      if (.not. output_bitmap(i,j)) then
        output_u_data(i,j) = -9999.
        output_v_data(i,j) = -9999.
      endif
   enddo
   enddo
 endif

!-------------------------------------------------------------------------
! Compared data from ipolatev to its corresponding baseline
! data.  
!-------------------------------------------------------------------------

 if (kind(output_u_data) == 8) then
   output_file = "../baseline_data/vector/8_byte_bin/grid" // trim(grid) // ".opt" // interp_opt // ".bin_8"
 else
   output_file = "../baseline_data/vector/4_byte_bin/grid" // trim(grid) // ".opt" // interp_opt // ".bin_4"
 endif

 allocate (baseline_u_data(i_output,j_output))
 allocate (baseline_v_data(i_output,j_output))

 open (12, file=output_file, access="direct", err=38, recl=mo*4)
 read (12, err=38, rec=1) baseline_u_data
 read (12, err=38, rec=2) baseline_v_data
 close (12)

 avg_u_diff=0.0
 max_u_diff=0.0
 num_upts_diff=0
 avg_v_diff=0.0
 max_v_diff=0.0
 num_vpts_diff=0

 do j = 1, j_output
 do i = 1, i_output
   output_data4 = real(output_u_data(i,j),4)
   if (abs(output_data4 - baseline_u_data(i,j)) > 0.0001) then
     avg_u_diff = avg_u_diff + abs(output_data4-baseline_u_data(i,j))
     num_upts_diff = num_upts_diff + 1
     if (abs(output_data4-baseline_u_data(i,j)) > abs(max_u_diff))then
       max_u_diff = output_data4-baseline_u_data(i,j)
     endif
   endif
   output_data4 = real(output_v_data(i,j),4)
   if (abs(output_data4 - baseline_v_data(i,j)) > 0.0001) then
     avg_v_diff = avg_v_diff + abs(output_data4-baseline_v_data(i,j))
     num_vpts_diff = num_vpts_diff + 1
     if (abs(output_data4-baseline_v_data(i,j)) > abs(max_v_diff))then
       max_v_diff = output_data4-baseline_v_data(i,j)
     endif
   endif
 enddo
 enddo

 print*,'- MAX/MIN OF U-WIND DATA: ', maxval(output_u_data),minval(output_u_data)
 print*,'- NUMBER OF PTS DIFFERENT: ',num_upts_diff
 print*,'- PERCENT OF TOTAL: ',(float(num_upts_diff)/float(i_output*j_output))*100.
 print*,'- MAX DIFFERENCE: ', max_u_diff
 if (num_upts_diff > 0) then
   avg_u_diff = avg_u_diff / float(num_upts_diff)
 endif
 print*,'- AVG DIFFERENCE: ', avg_u_diff

 print*,'- MAX/MIN OF V-WIND DATA: ', maxval(output_v_data),minval(output_v_data)
 print*,'- NUMBER OF PTS DIFFERENT: ',num_vpts_diff
 print*,'- PERCENT OF TOTAL: ',(float(num_vpts_diff)/float(i_output*j_output))*100.
 print*,'- MAX DIFFERENCE: ', max_v_diff
 if (num_vpts_diff > 0) then
   avg_v_diff = avg_v_diff / float(num_vpts_diff)
 endif
 print*,'- AVG DIFFERENCE: ', avg_v_diff

 deallocate (baseline_u_data, baseline_v_data)
 deallocate (output_rlat, output_rlon, output_u_data, output_bitmap)
 deallocate (output_v_data, output_crot, output_srot)

 return

 38 continue
 print*,'- ERROR READING BASELINE DATA FILE'
 stop 7

 end subroutine interp
