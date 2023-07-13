! This is test code from the NCEPLIBS-ip project.
!
! This is a helper module for tests which load some GRIB2 input data
! to be interpolated.
!
! Kyle Gerheiser June, 2021

#if (LSIZE==D)
#define REALSIZE 8
#elif (LSIZE==4)
#define REALSIZE 4
#endif

module input_data_mod_grib2
  implicit none

  !------------------------------------------------------------------------
  ! Read the data that will be interpolated.  Data is a global one-degree
  ! grid of albedo with no bitmap.
  !------------------------------------------------------------------------

  integer,  public     :: i_input
  integer,  public     :: j_input

  integer, parameter, public     :: input_gdtnum=0
  integer, parameter, public     :: input_gdtlen=19
  integer, public                :: input_gdtmpl(input_gdtlen)
  integer, public                :: vector_input_gdtmpl(input_gdtlen)
  
  integer, parameter :: missing=4294967296

  real(KIND=REALSIZE), allocatable, public      :: input_data(:,:)
  real(KIND=REALSIZE), allocatable, public      :: input_u_data(:,:)
  real(KIND=REALSIZE), allocatable, public      :: input_v_data(:,:)

  logical*1, allocatable, public :: input_bitmap(:,:)
 
  data input_gdtmpl /6, 255, missing, 255, missing, 255, missing, &
       360, 180, 0, missing, -89500000, -180000000, &
       48, 89500000, 179000000, 1000000, 1000000, 64/

  
  data vector_input_gdtmpl /6, 255, missing, 255, missing, 255, missing, &
       360, 181, 0, missing, 90000000, 0, &
       48, -90000000, 359000000, 1000000, 1000000, 0/
  

contains

  subroutine read_input_data

    implicit none

    character*100      :: input_file

    integer            :: iret
    integer, parameter :: iunit=9

    real(kind=4), allocatable :: dummy(:,:)

    i_input = input_gdtmpl(8)
    j_input = input_gdtmpl(9)
    
    allocate(input_data(i_input, j_input), dummy(i_input, j_input), input_bitmap(i_input, j_input))
    
    input_file="data/input_data/scalar/global_snoalb.bin"
    open(iunit, file=input_file, access='direct', recl=i_input*j_input*4,  &
         iostat=iret)

    if (iret /= 0) then
       print*,'- BAD OPEN OF FILE, IRET IS ', iret
       stop 2
    end if

    read(iunit, rec=1, iostat=iret) dummy
    input_data=dummy

    if (iret /= 0) then
       print*,"- BAD READ OF DATA. IRET IS ", iret
       stop 4
    end if

    close (iunit)

    input_bitmap=.true.

    return

  end subroutine read_input_data

  subroutine read_vector_input_data()
    implicit none

    character*100             :: input_file

    integer                   :: iret
    integer, parameter        :: iunit=9

    real(kind=4), allocatable :: dummy(:,:)

    i_input = vector_input_gdtmpl(8)
    j_input = vector_input_gdtmpl(9)

    print *, i_input, j_input

    input_file="data/input_data/vector/global_uv_wind.bin"
    open(iunit, file=input_file, access='direct', recl=i_input*j_input*4, iostat=iret)

    if (iret /= 0) then
       print*,'- BAD OPEN OF INPUT DATA FILE, IRET IS ', iret
       stop 2
    end if

    allocate(dummy(i_input,j_input))

    read(iunit, rec=1, iostat=iret) dummy

    if (iret /= 0) then
       print*,"- BAD READ OF INPUT U-WIND DATA. IRET IS ", iret
       stop 4
    endif

    allocate(input_u_data(i_input,j_input))
    input_u_data=dummy

    read(iunit, rec=2, iostat=iret) dummy

    if (iret /= 0) then
       print*,"- BAD READ OF INPUT V-WIND DATA. IRET IS ", iret
       stop 4
    end if

    allocate(input_v_data(i_input,j_input))
    input_v_data=dummy

    close (iunit)

    deallocate(dummy)

    allocate(input_bitmap(i_input,j_input))
    input_bitmap=.true.

    ! impose pure north wind
    !input_u_data=0.
    !input_v_data=-1.

    return
  end subroutine read_vector_input_data

end module input_data_mod_grib2
