  use ip_mod
  implicit none

contains

  subroutine interp(grid, interp_opt)
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
    !  203 -  rotated lat/lon e-staggered (number refers to the old
    !                                      grib 1 gds octet 6)
    !  205 -  rotated lat/lon b-staggered (number refers to the old
    !                                      gds octet 6)
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
    ! The interpolated data is compared against a baseline binary
    ! file.  Any differences are written to standard output.
    !-------------------------------------------------------------------------
#if(LSIZE==4)
    use input_data_mod_grib1_4, &
#elif(LSIZE==D)
    use input_data_mod_grib1_d, &
#elif(LSIZE==8)
    use input_data_mod_grib1_8, &
#endif
      only: input_data, input_kgds, input_bitmap, i_input, j_input

    implicit none

    character(*), intent(in) :: grid, interp_opt
    character*100 :: baseline_file


    integer     :: ip, ipopt(20), output_kgds(200)
    integer     :: km, ibi(1), mi, iret, i, j, ibi_scalar=0
    integer     :: i_output, j_output, mo, no, ibo(1), ibo_scalar
    integer                   :: num_pts_diff, which_func, ntol

    logical*1, allocatable    :: output_bitmap(:,:)

    real, allocatable         :: output_rlat(:), output_rlon(:)
    real, allocatable         :: output_data(:,:)
    real(kind=4), allocatable :: baseline_data(:,:)
    real                      :: avgdiff, maxdiff
    real(kind=4)              :: output_data4
    real                      :: abstol

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

#if (LSIZE==4)
    abstol=0.05
    ntol = 10
#else
    abstol=0.0001
    ntol = 0
#endif

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
       ibi   = 0     ! can't use bitmap with spectral
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

    allocate(output_rlat(i_output * j_output))
    allocate(output_rlon(i_output * j_output))

    allocate (output_data(i_output,j_output))
    allocate (output_bitmap(i_output,j_output))

    allocate (baseline_data(i_output,j_output))

    do which_func=1,2
        if (which_func .eq. 1) then
            call ipolates(ip, ipopt, input_kgds, output_kgds, mi, mo,&
                 km, ibi, input_bitmap, input_data, &
                 no, output_rlat, output_rlon, ibo, output_bitmap, output_data, iret)
        elseif (which_func .eq. 2) then
            call ipolates_grib1_single_field(ip, ipopt, input_kgds, output_kgds, mi, mo,&
                 km, ibi_scalar, input_bitmap, input_data, &
                 no, output_rlat, output_rlon, ibo_scalar, output_bitmap, output_data, iret)
        endif

! Uncomment to generate new baseline file:
!        open (13, file="grid"//trim(grid)//".opt"//trim(interp_opt)//".bin_"//"LSIZE", access="direct", recl=mo*4)
!        write (13, rec=1) real(output_data, kind=4)
!        close (13)

        if (iret /= 0) then
           print*,'- BAD STATUS FROM IPOLATES: ', iret
           stop 6
        end if

        if (no /= mo) then
           print*,'- ERROR: WRONG # OF POINTS RETURNED FROM IPOLATES'
           stop 7
        end if

        !print*,'- SUCCESSFULL CALL TO IPOLATES'

        do j = 1, j_output
           do i = 1, i_output
              if (.not. output_bitmap(i,j)) then
                 output_data(i,j) = -9.
              endif
           enddo
        enddo

        if (kind(output_data) == 8) then
           baseline_file = "./data/baseline_data/grib1/scalar/8_byte_bin/grid" // &
               trim(grid) // ".opt" // trim(interp_opt) // ".bin_8"
        else
           baseline_file = "./data/baseline_data/grib1/scalar/4_byte_bin/grid" // &
               trim(grid) // ".opt" // trim(interp_opt) // ".bin_4"
        endif

        open (12, file=baseline_file, access="direct", err=38, recl=mo*4)
        read (12, err=38, rec=1) baseline_data
        close (12)

        avgdiff=0.0
        maxdiff=0.0
        num_pts_diff=0

        do j = 1, j_output
           do i = 1, i_output
              output_data4 = real(output_data(i,j),4)
              if ( abs(output_data4 - baseline_data(i,j)) > abstol) then
                 avgdiff = avgdiff + abs(output_data4-baseline_data(i,j))
                 num_pts_diff = num_pts_diff + 1
                 if (abs(output_data4-baseline_data(i,j)) > abs(maxdiff))then
                    maxdiff = output_data4-baseline_data(i,j)
                 endif
              endif
           enddo
        enddo

        print*,'- MAX/MIN OF DATA: ', maxval(output_data),minval(output_data)
        print*,'- NUMBER OF PTS DIFFERENT: ',num_pts_diff
        print*,'- PERCENT OF TOTAL: ',(float(num_pts_diff)/float(i_output*j_output))*100.
        print*,'- MAX DIFFERENCE: ', maxdiff
        if (num_pts_diff > 0) then
           avgdiff = avgdiff / float(num_pts_diff)
        endif
        print*,'- AVG DIFFERENCE: ', avgdiff

        if (num_pts_diff > ntol) then
           print *, "# DIFFERING POINTS: ", num_pts_diff
           error stop
        endif
    enddo ! which_func

    deallocate (output_rlat, output_rlon, output_data, output_bitmap, baseline_data)

    return

38  continue

    print*,'-ERROR READING BASELINE DATA FILE FILE.'
    stop 77

  end subroutine interp



  subroutine interp_vector(grid, interp_opt)

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
    !  203 -  rotated lat/lon e-staggered (number refers to the old
    !                                      gds octet 6)
    !  205 -  rotated lat/lon b-staggered (number refers to the old
    !                                      gds octet 6)
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

#if(LSIZE==4)
    use input_data_mod_grib1_4, &
#elif(LSIZE==D)
    use input_data_mod_grib1_d, &
#elif(LSIZE==8)
    use input_data_mod_grib1_8, &
#endif
      only: input_u_data, input_v_data, vector_input_kgds, input_bitmap, &
         i_input, j_input

    implicit none

    character(*), intent(in) :: grid, interp_opt

    character*100             :: baseline_file

    integer     :: ip, ipopt(20), output_kgds(200)
    integer     :: km, ibi(1), mi, iret, i, j, which_func
    integer     :: i_output, j_output, mo, no, ibo(1)
    integer     :: ibi_scalar=0, ibo_scalar
    integer                   :: num_upts_diff, num_vpts_diff, ntol

    logical*1, allocatable    :: output_bitmap(:,:)

    real, allocatable         :: output_rlat(:), output_rlon(:)
    real, allocatable         :: output_crot(:), output_srot(:)
    real, allocatable         :: output_u_data(:,:), output_v_data(:,:)
    real                      :: avg_u_diff, avg_v_diff
    real                      :: max_u_diff, max_v_diff
    real(kind=4)              :: output_data4
    real(kind=4), allocatable :: baseline_u_data(:,:)
    real(kind=4), allocatable :: baseline_v_data(:,:)
    real                      :: abstol

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

#if (LSIZE==4)
    abstol=0.1
    ntol = 15
#else
    abstol=0.0001
    ntol = 0
#endif

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

    allocate (output_rlat(i_output * j_output))
    allocate (output_rlon(i_output * j_output))
    allocate (output_u_data(i_output,j_output))
    allocate (output_v_data(i_output,j_output))
    allocate (output_bitmap(i_output,j_output))
    allocate (output_srot(i_output * j_output))
    allocate (output_crot(i_output * j_output))
    allocate (baseline_u_data(i_output,j_output))
    allocate (baseline_v_data(i_output,j_output))

    do which_func=1,2
        if (which_func .eq. 1) then
            call ipolatev(ip, ipopt, vector_input_kgds, output_kgds, mi, mo,   &
                 km, ibi, input_bitmap, input_u_data, input_v_data,  &
                 no, output_rlat, output_rlon, output_crot, output_srot, &
                 ibo, output_bitmap, output_u_data, output_v_data, iret)
        elseif (which_func .eq. 2) then
            call ipolatev_grib1_single_field(ip, ipopt, vector_input_kgds, output_kgds, mi, mo,   &
                 km, ibi_scalar, input_bitmap, input_u_data, input_v_data,  &
                 no, output_rlat, output_rlon, output_crot, output_srot, &
                 ibo_scalar, output_bitmap, output_u_data, output_v_data, iret)
        endif

! Uncomment to generate new baseline file:
!        open (13, file="grid"//trim(grid)//".opt"//trim(interp_opt)//".bin_"//"LSIZE", access="direct", recl=mo*4)
!        write (13, rec=1) real(output_u_data, kind=4)
!        write (13, rec=2) real(output_v_data, kind=4)
!        close (13)

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
           baseline_file = "./data/baseline_data/grib1/vector/8_byte_bin/grid" // &
               trim(grid) // ".opt" // trim(interp_opt) // ".bin_8"
        else
           baseline_file = "./data/baseline_data/grib1/vector/4_byte_bin/grid" // &
               trim(grid) // ".opt" // trim(interp_opt) // ".bin_4"
        endif

        print *,"baseline_file = ", baseline_file

        open (12, file=baseline_file, access="direct", err=38, recl=mo*4)
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
              if (abs(output_data4 - baseline_u_data(i,j)) > abstol) then
                 avg_u_diff = avg_u_diff + abs(output_data4-baseline_u_data(i,j))
                 num_upts_diff = num_upts_diff + 1
                 if (abs(output_data4-baseline_u_data(i,j)) > abs(max_u_diff))then
                    max_u_diff = output_data4-baseline_u_data(i,j)
                 endif
              endif
              output_data4 = real(output_v_data(i,j),4)
              if (abs(output_data4 - baseline_v_data(i,j)) > abstol) then
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

        if (num_vpts_diff + num_upts_diff > ntol) then
           print *, "# DIFFERING POINTS: ", num_upts_diff + num_vpts_diff
           error stop
        endif
    enddo ! which_func

    deallocate(input_bitmap, input_u_data, input_v_data)
    deallocate (baseline_u_data, baseline_v_data)
    deallocate (output_rlat, output_rlon, output_u_data, output_bitmap)
    deallocate (output_v_data, output_crot, output_srot)

    return

38  continue
    print*,'- ERROR READING BASELINE DATA FILE'
    stop 7

  end subroutine interp_vector
