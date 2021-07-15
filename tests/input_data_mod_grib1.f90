module input_data_mod_grib1
  implicit none

  !------------------------------------------------------------------------
  ! Read the data that will be interpolated.  Data is a global one-degree
  ! grid of albedo with no bitmap.
  !------------------------------------------------------------------------

  integer,  public     :: i_input
  integer,  public     :: j_input

  integer, public                :: input_kgds(200)
  integer, public                :: vector_input_kgds(200)

  !  integer, parameter, public     :: input_gdtnum=0
  !  integer, parameter, public     :: input_gdtlen=19
  ! integer, public                :: input_gdtmpl(input_gdtlen)
  ! integer, public                :: vector_input_gdtmpl(input_gdtlen)
  
  ! integer, parameter :: missing=b'11111111111111111111111111111111'
  
  real, allocatable, public      :: input_data(:,:)
  logical*1, allocatable, public :: input_bitmap(:,:)
  real, allocatable, public      :: input_u_data(:,:)
  real, allocatable, public      :: input_v_data(:,:)


  data input_kgds /0,  360,  180, -89500, -180000, 128,  &
       89500, 179000, 1000, 1000,  64,  0,   &
       6*-1, 0, 255, 180*-1/

  data vector_input_kgds /0,  360,  181,  90000,   0,   128,   &
       -90000,  -1000,  1000,  1000,   0,   0,  &
       6*-1,   0,   255,  180*-1/

  
contains

  subroutine read_input_data

    implicit none

    character*100      :: input_file

    integer            :: iret
    integer, parameter :: iunit=9

    real(kind=4), allocatable :: dummy(:,:)

    i_input = input_kgds(2)
    j_input = input_kgds(3)

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

    i_input = vector_input_kgds(2)
    j_input = vector_input_kgds(3)

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

end module input_data_mod_grib1
