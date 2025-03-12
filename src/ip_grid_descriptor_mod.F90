!> @file
!! @brief Users derived type grid descriptor objects to abstract away
!! the raw GRIB1 and GRIB2 grid definitions.
!! @author Kyle Gerheiser

!> @brief Users derived type grid descriptor objects to abstract away
!! the raw GRIB1 and GRIB2 grid definitions.
!!
!! An abstract ip_grid_descriptor class is subclassed by the GRIB1 and
!! GRIB2 descriptor types which contain the raw integer grib arrays. A
!! comparison operator to test whether two grid descriptors are the
!! same is also included.
!!
!! @author Kyle Gerheiser @date July 2021
module ip_grid_descriptor_mod
  implicit none

  private

  public :: ip_grid_descriptor
  public :: grib1_descriptor, grib2_descriptor
  public :: init_descriptor, init_grib1_descriptor, init_grib2_descriptor

  public :: operator(==)

  !> Abstract descriptor object which represents a grib1 or grib2 descriptor.
  !! @date July 2021
  type, abstract :: ip_grid_descriptor
     integer :: grid_num !< Integer representing the grid type (see *_GRID_ID_GRIB1/2 in ip_grid_mod).
   contains
     !> Test whether two grid descriptors are the same.
     !> @return N/A @copydoc ip_grid_descriptor_mod::is_same_grid
     procedure :: is_same_grid
  end type ip_grid_descriptor

  !> Descriptor representing a grib1 grib descriptor section (GDS)
  !> with an integer array
  !! @date July 2021
  type, extends(ip_grid_descriptor) :: grib1_descriptor
     integer :: gds(200) !< Grib-1 grib descriptor section (GDS)
   contains
     !> Test whether two grid descriptors are the same.
     !> @return N/A @copydoc ip_grid_descriptor_mod::is_same_grid_grib1
     procedure :: is_same_grid_grib1
  end type grib1_descriptor

  !> Grib-2 descriptor containing a grib2 GDT represented by an integer array
  !! @date July 2021
  type, extends(ip_grid_descriptor) :: grib2_descriptor
     integer :: gdt_num !< Grid number which represents grid type.
     integer :: gdt_len !< Length of the template.
     integer, allocatable :: gdt_tmpl(:) !< Grib-2 grid definition template.
   contains
     !> Test whether two grid descriptors are the same.
     !> @return N/A @copydoc ip_grid_descriptor_mod::is_same_grid_grib2
     procedure :: is_same_grid_grib2
  end type grib2_descriptor

  interface operator (==)
     module procedure is_same_grid
  end interface operator (==)

  interface init_descriptor
     module procedure init_grib1_descriptor
     module procedure init_grib2_descriptor
  end interface init_descriptor

contains

  !> Initialize grib-1 descriptor from integer grid definition section (GDS).
  !! @param[in] gds Grib-1 grid definition section.
  !!
  !! @return Initialized Grib-1 descriptor.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  function init_grib1_descriptor(gds) result(desc)
    type(grib1_descriptor) :: desc
    integer, intent(in) :: gds(:)
    desc%gds = gds
    desc%grid_num = gds(1)

    !call desc%decode_template()

  end function init_grib1_descriptor

  !> Initialize grib-2 descriptor from integer grid definition template (GDT).
  !! @param[in] gdt_num Grib-2 grid number.
  !! @param[in] gdt_len Lenght of the grid definition template.
  !! @param[in] gdt_tmpl Grib-2 grid definition template.
  !!
  !! @return Initialized Grib-2 descriptor.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  function init_grib2_descriptor(gdt_num, gdt_len, gdt_tmpl) result(desc)
    type(grib2_descriptor) :: desc
    integer, intent(in) :: gdt_num, gdt_len, gdt_tmpl(:)
    desc%grid_num = gdt_num

    desc%gdt_num = gdt_num
    desc%gdt_len = gdt_len
    allocate(desc%gdt_tmpl(gdt_len))
    desc%gdt_tmpl = gdt_tmpl

    !call desc%decode_template()

  end function init_grib2_descriptor

  !> Test whether two grid descriptors are the same.
  !! @param[in] grid1 An ip_grid_descriptor.
  !! @param[in] grid2 Another ip_grid_descriptor.
  !!
  !! @return True if the grids are the same, false if they are not.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  logical function is_same_grid(grid1, grid2)
    class(ip_grid_descriptor), intent(in) :: grid1, grid2

    select type(grid1)
    type is(grib1_descriptor)
       select type(grid2)
       type is(grib1_descriptor)
          is_same_grid = grid1%is_same_grid_grib1(grid2)
          class default
          is_same_grid = .false.
       end select
    type is(grib2_descriptor)
       select type(grid2)
       type is(grib2_descriptor)
          is_same_grid = grid1%is_same_grid_grib2(grid2)
          class default
          is_same_grid = .false.
       end select
    end select

  end function is_same_grid

  !> Test whether two grib1_descriptors are the same.
  !! @param[in] self The grib1_descriptor which this routine was called on.
  !! @param[in] grid_desc A grib1_descriptor to compare.
  !!
  !! @return  True if the grids are the same, false if they are not.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  logical function is_same_grid_grib1(self, grid_desc) result(same_grid)
    class(grib1_descriptor), intent(in) :: self, grid_desc

    if (all(self%gds == grid_desc%gds)) then
       same_grid = .true.
    else
       same_grid = .false.
    end if

  end function is_same_grid_grib1

  !> Test whether two grib2_descriptors are the same.
  !! @param[in] self The grib2_descriptor which this routine was called on.
  !! @param[in] grid_desc grib2_descriptor to compare.
  !!
  !! @return  True if the grids are the same, false if they are not.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  logical function is_same_grid_grib2(self, grid_desc) result(same_grid)
    class(grib2_descriptor), intent(in) :: self, grid_desc

    same_grid = .false.
    if (self%grid_num == grid_desc%grid_num) then
       if (self%gdt_len == grid_desc%gdt_len) then
          if (all(self%gdt_tmpl == grid_desc%gdt_tmpl)) then
             same_grid = .true.
          end if
       end if
    end if

  end function is_same_grid_grib2



  ! subroutine decode_template_grib1(self)
  !   type(grib1_descriptor), intent(inout) :: self

  !   integer                            :: im, jm, iwrap, jg
  !   integer                            :: iscan, kscan, nscan, nscan_field_pos
  !   integer                            :: jwrap1, jwrap2

  !   real                               :: dlat, dlon
  !   real                               :: rlat1, rlat2
  !   real                               :: rlon1, rlon2

  !   ! set a default value for this to check if it has changed
  !   ! for rotated grids nscan is set = 3, but in other places it's not, so use this special value
  !   ! just for field_position routine
  !   nscan_field_pos = -1

  !   im = self%gds(2)
  !   jm = self%gds(3)
  !   iwrap = 0
  !   jwrap1 = 0
  !   jwrap2 = 0
  !   nscan = mod(self%gds(11) / 32, 2)
  !   kscan = 0

  !   select case(self%gds(1))
  !   case(0)
  !      rlon1=self%gds(5)*1.e-3
  !      rlon2=self%gds(8)*1.e-3
  !      iscan=mod(self%gds(11)/128,2)
  !      if(iscan.eq.0) then
  !         dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
  !      else
  !         dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
  !      endif
  !      iwrap=nint(360/abs(dlon))
  !      if(im.lt.iwrap) iwrap=0
  !      if(iwrap.gt.0.and.mod(iwrap,2).eq.0) then
  !         rlat1=self%gds(4)*1.e-3
  !         rlat2=self%gds(7)*1.e-3
  !         dlat=abs(rlat2-rlat1)/(jm-1)
  !         if(abs(rlat1).gt.90-0.25*dlat) then
  !            jwrap1=2
  !         elseif(abs(rlat1).gt.90-0.75*dlat) then
  !            jwrap1=1
  !         endif
  !         if(abs(rlat2).gt.90-0.25*dlat) then
  !            jwrap2=2*jm
  !         elseif(abs(rlat2).gt.90-0.75*dlat) then
  !            jwrap2=2*jm+1
  !         endif
  !      endif
  !   case(1)
  !      rlon1=self%gds(5)*1.e-3
  !      rlon2=self%gds(8)*1.e-3
  !      iscan=mod(self%gds(11)/128,2)
  !      if(iscan.eq.0) then
  !         dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
  !      else
  !         dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
  !      endif
  !      iwrap=nint(360/abs(dlon))
  !      if(im.lt.iwrap) iwrap=0
  !   case(4)
  !      rlon1=self%gds(5)*1.e-3
  !      rlon2=self%gds(8)*1.e-3
  !      iscan=mod(self%gds(11)/128,2)
  !      if(iscan.eq.0) then
  !         dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
  !      else
  !         dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
  !      endif
  !      iwrap=nint(360/abs(dlon))
  !      if(im.lt.iwrap) iwrap=0
  !      if(iwrap.gt.0.and.mod(iwrap,2).eq.0) then
  !         jg=self%gds(10)*2
  !         if(jm.eq.jg) then
  !            jwrap1=1
  !            jwrap2=2*jm+1
  !         endif
  !      endif
  !   case(203)
  !      ! why is nscan set to 3 for this grid?
  !      nscan_field_pos = 3
  !      kscan=mod(self%gds(11) / 256, 2)
  !   case default
  !      print *, "grib1 grid type ", self%gds(1), " not recognized"
  !      error stop
  !   end select

  !   self%im = im
  !   self%jm = jm
  !   self%nm = im * jm
  !   self%iwrap = iwrap
  !   self%jwrap1 = jwrap1
  !   self%jwrap2 = jwrap2
  !   self%nscan = nscan

  !   if (nscan_field_pos == -1) then
  !      ! just use regular value of nscan
  !      self%nscan_field_pos = nscan
  !   else
  !      ! nscan = 3 for rotated grids and is set to 3 for use in field_position
  !      self%nscan_field_pos = nscan_field_pos
  !   end if

  !   self%kscan = kscan

  ! end subroutine decode_template_grib1

  ! subroutine decode_template_grib2(self)
  !   type(grib2_descriptor), intent(inout) :: self

  !   integer                            :: im, jm, iwrap, jg
  !   integer                            :: i_offset_odd, i_offset_even
  !   integer                            :: iscan, kscan, nscan, nscan_field_pos
  !   integer                            :: jwrap1, jwrap2, iscale

  !   real                               :: dlat, dlon
  !   real                               :: rlat1, rlat2
  !   real                               :: rlon1, rlon2

  !   nscan_field_pos = -1

  !   select case(self%gdt_num)
  !   ! EQUIDISTANT CYLINDRICAL
  !   case(0)
  !      im=self%gdt_tmpl(8)
  !      jm=self%gdt_tmpl(9)
  !      nscan=mod(self%gdt_tmpl(19)/32,2)
  !      kscan=0
  !      iscale=self%gdt_tmpl(10)*self%gdt_tmpl(11)
  !      if(iscale==0) iscale=10**6
  !      rlon1=float(self%gdt_tmpl(13))/float(iscale)
  !      rlon2=float(self%gdt_tmpl(16))/float(iscale)
  !      iscan=mod(self%gdt_tmpl(19)/128,2)
  !      if(iscan.eq.0) then
  !         dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
  !      else
  !         dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
  !      endif
  !      iwrap=nint(360/abs(dlon))
  !      if(im.lt.iwrap) iwrap=0
  !      jwrap1=0
  !      jwrap2=0
  !      if(iwrap.gt.0.and.mod(iwrap,2).eq.0) then
  !         rlat1=float(self%gdt_tmpl(12))/float(iscale)
  !         rlat2=float(self%gdt_tmpl(15))/float(iscale)
  !         dlat=abs(rlat2-rlat1)/(jm-1)
  !         if(abs(rlat1).gt.90-0.25*dlat) then
  !            jwrap1=2
  !         elseif(abs(rlat1).gt.90-0.75*dlat) then
  !            jwrap1=1
  !         endif
  !         if(abs(rlat2).gt.90-0.25*dlat) then
  !            jwrap2=2*jm
  !         elseif(abs(rlat2).gt.90-0.75*dlat) then
  !            jwrap2=2*jm+1
  !         endif
  !      endif
  !   case(1)
  !      i_offset_odd=mod(self%gdt_tmpl(19)/8,2)
  !      i_offset_even=mod(self%gdt_tmpl(19)/4,2)
  !      im=self%gdt_tmpl(8)
  !      jm=self%gdt_tmpl(9)
  !      iwrap=0
  !      jwrap1=0
  !      jwrap2=0
  !      kscan=0
  !      nscan=mod(self%gdt_tmpl(19)/32,2)
  !      if(i_offset_odd/=i_offset_even)then
  !         kscan=i_offset_odd
  !         nscan_field_pos=3
  !      endif
  !   ! MERCATOR CYLINDRICAL
  !   case(10)
  !      im=self%gdt_tmpl(8)
  !      jm=self%gdt_tmpl(9)
  !      rlon1=float(self%gdt_tmpl(11))*1.0e-6
  !      rlon2=float(self%gdt_tmpl(15))*1.0e-6
  !      iscan=mod(self%gdt_tmpl(16)/128,2)
  !      if(iscan.eq.0) then
  !         dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
  !      else
  !         dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
  !      endif
  !      iwrap=nint(360/abs(dlon))
  !      if(im.lt.iwrap) iwrap=0
  !      jwrap1=0
  !      jwrap2=0
  !      kscan=0
  !      nscan=mod(self%gdt_tmpl(16)/32,2)
  !   ! POLAR STEREOGRAPHIC AZIMUTHAL
  !   case(20)
  !      im=self%gdt_tmpl(8)
  !      jm=self%gdt_tmpl(9)
  !      nscan=mod(self%gdt_tmpl(18)/32,2)
  !      iwrap=0
  !      jwrap1=0
  !      jwrap2=0
  !      kscan=0
  !   ! LAMBERT CONFORMAL CONICAL
  !   case(30)
  !      im=self%gdt_tmpl(8)
  !      jm=self%gdt_tmpl(9)
  !      nscan=mod(self%gdt_tmpl(18)/32,2)
  !      iwrap=0
  !      jwrap1=0
  !      jwrap2=0
  !      kscan=0
  !   ! GAUSSIAN CYLINDRICAL
  !   case(40)
  !      im=self%gdt_tmpl(8)
  !      jm=self%gdt_tmpl(9)
  !      iscale=self%gdt_tmpl(10)*self%gdt_tmpl(11)
  !      if(iscale==0) iscale=10**6
  !      rlon1=float(self%gdt_tmpl(13))/float(iscale)
  !      rlon2=float(self%gdt_tmpl(16))/float(iscale)
  !      iscan=mod(self%gdt_tmpl(19)/128,2)
  !      if(iscan.eq.0) then
  !         dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
  !      else
  !         dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
  !      endif
  !      iwrap=nint(360/abs(dlon))
  !      if(im.lt.iwrap) iwrap=0
  !      jwrap1=0
  !      jwrap2=0
  !      if(iwrap.gt.0.and.mod(iwrap,2).eq.0) then
  !         jg=self%gdt_tmpl(18)*2
  !         if(jm.eq.jg) then
  !            jwrap1=1
  !            jwrap2=2*jm+1
  !         endif
  !      endif
  !      nscan=mod(self%gdt_tmpl(19)/32,2)
  !      kscan=0
  !   case default
  !      print *, "gdt_num ", self%gdt_num, " not recognized"
  !      error stop
  !   end select

  !   self%im = im
  !   self%jm = jm
  !   self%nm = im*jm
  !   self%iwrap = iwrap
  !   self%jwrap1 = jwrap1
  !   self%jwrap2 = jwrap2
  !   self%nscan = nscan
  !   if (nscan_field_pos == -1) then
  !      self%nscan_field_pos = nscan
  !   else
  !      self%nscan_field_pos = nscan_field_pos
  !   end if
  !   self%kscan = kscan

  ! end subroutine decode_template_grib2

  ! subroutine earth_radius(self, radius, eccen_squared)
  !   class(ip_grid_descriptor), intent(in) :: self
  !   real, intent(out) :: radius, eccen_squared

  !   real :: flat, major_axis, minor_axis
  !   logical :: elliptical

  !   select type(self)
  !   type is(grib1_descriptor)
  !      associate(gds => self%gds)
  !        select case(gds(1))
  !        case(0)
  !           radius = 6.3712d6
  !           eccen_squared = 0d0
  !        case(1)
  !           radius = 6.3712d6
  !           eccen_squared = 0d0
  !        case(3)
  !           radius = 6.3712d6
  !           eccen_squared = 0d0
  !        case(4)
  !           radius = 6.3712d6
  !           eccen_squared = 0d0
  !        case(5)
  !           elliptical = mod(gds(6) / 64, 2) == 1
  !           if (.not. elliptical) then
  !              radius = 6.3712d6
  !              eccen_squared = 0d0
  !           else
  !              radius = 6.378137E6 ! WGS-84
  !              eccen_squared = 0.00669437999013d0
  !           end if
  !        case(203)
  !           radius = 6.3712d6
  !           eccen_squared = 0d0
  !        case(205)
  !           radius = 6.3712d6
  !           eccen_squared = 0d0
  !        case default
  !           print *, "grib1 grid not recognized. gds(1) = ", gds(1)
  !           error stop
  !        end select
  !      end associate
  !   type is(grib2_descriptor)
  !      associate(gdt_tmpl => self%gdt_tmpl)
  !        select case (gdt_tmpl(1))
  !        case (0)
  !           radius = 6367470.0
  !           eccen_squared = 0.0
  !        case (1)  ! user specified spherical
  !           radius = float(gdt_tmpl(3))/float(10**gdt_tmpl(2))
  !           eccen_squared = 0.0
  !        case (2)  ! iau 1965
  !           radius = 6378160.0      ! semi major axis
  !           flat = 1.0/297.0      ! flattening
  !           eccen_squared = (2.0*flat) - (flat**2)
  !        case (3)  ! user specified elliptical (km)
  !           major_axis = float(gdt_tmpl(5))/float(10**gdt_tmpl(4))
  !           major_axis = major_axis * 1000.0
  !           minor_axis = float(gdt_tmpl(7))/float(10**gdt_tmpl(6))
  !           minor_axis = minor_axis * 1000.0
  !           eccen_squared = 1.0 - (minor_axis**2 / major_axis**2)
  !           radius = major_axis
  !        case (4)  ! iag-grs80 model
  !           radius = 6378137.0      ! semi major axis
  !           flat = 1.0/298.2572   ! flattening
  !           eccen_squared = (2.0*flat) - (flat**2)
  !        case (5)  ! wgs84 datum
  !           radius = 6378137.0      ! semi major axis
  !           eccen_squared = 0.00669437999013
  !        case (6)
  !           radius = 6371229.0
  !           eccen_squared = 0.0
  !        case (7)  ! user specified elliptical (m)
  !           major_axis = float(gdt_tmpl(5))/float(10**gdt_tmpl(4))
  !           minor_axis = float(gdt_tmpl(7))/float(10**gdt_tmpl(6))
  !           eccen_squared = 1.0 - (minor_axis**2 / major_axis**2)
  !           radius = major_axis
  !        case (8)
  !           radius = 6371200.0
  !           eccen_squared = 0.0
  !        case default
  !           radius = -9999.
  !           eccen_squared = -9999.
  !           error stop
  !        end select
  !      end associate
  !   class default
  !      print *, "unknown descriptor type"
  !      error stop
  !   end select
  ! end subroutine earth_radius

end module ip_grid_descriptor_mod
