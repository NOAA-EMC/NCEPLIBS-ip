 subroutine interp

!-------------------------------------------------------------------------
! call the vector polates routines to interpolate the input data
! using all available interpolation methods (neighbor, bilinear, etc.)
! several output grids of various map projections are tested.
!-------------------------------------------------------------------------

 use get_input_data, only : input_u_data, input_v_data, &
                            input_kgds, &
                            input_bitmap, &
                            i_input, j_input

 implicit none

 character*1             :: interp_opt
 character*3             :: grid
 character*100           :: output_file

 integer*4   :: i1
 integer     :: ip, ipopt(20), output_kgds(200)
 integer     :: km, ibi, mi, iret, i, j
 integer     :: i_output, j_output, mo, no, ibo

 logical*1, allocatable :: output_bitmap(:,:)

 real, allocatable :: output_rlat(:,:), output_rlon(:,:)
 real, allocatable :: output_crot(:,:), output_srot(:,:)
 real, allocatable :: output_u_data(:,:), output_v_data(:,:)

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
     ibi   = 1  ! use input bitmap
! bi-cubic
   case ('1')          
     ip = 1
     ipopt = 0  
     ibi   = 1  ! use input bitmap
! neighbor
   case ('2')       
     ip = 2
     ipopt = 0  
     ibi   = 1  ! use input bitmap
! budget
   case ('3')       
     ip = 3
     ipopt = -1  
     ibi   = 1  ! use input bitmap
! spectral 
   case ('4')         
     ip = 4
     ipopt(1) = 0  ! triangular
     ipopt(2) = -1 ! code chooses wave number
     ipopt(3:20)=0
     ibi   = 0     ! can't use input bitmap with spectral
! neighbor-budget
   case ('6')            
     ip = 6
     ipopt = -1  
     ibi   = 1  ! use input bitmap
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

 print*,'- SUCCESSFULL CALL TO IPOLATES'

! polatev4 does not always initialize the ibo or output_bitmap variables,
! so they can't be used.  this should be fixed.
! according to the comments, polatev4 only outputs a bitmap in areas
! that extend outside the domain of the input grid. 
 
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

 output_file = "./grid" // trim(grid) // ".opt" // interp_opt // ".bin" 
 open (12, file=output_file, access="direct", recl=mo*4)
 write (12, rec=1) real(output_u_data,4)
 write (12, rec=2) real(output_v_data,4)
 close (12)

 deallocate (output_rlat, output_rlon, output_u_data, output_bitmap)
 deallocate (output_v_data, output_crot, output_srot)

 return
 end subroutine interp
