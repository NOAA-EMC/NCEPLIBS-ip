!> @file
!! @brief Abstract ip_grid type.
!!
!! @author Kyle Gerheiser @date July 2021

!> Abstract ip_grid type.
!!
!! @author Kyle Gerheiser
!! @date July 2021
module ip_grid_mod
  use ip_grid_descriptor_mod
  implicit none

  integer, public, parameter :: EQUID_CYLIND_GRID_ID_GRIB1 = 0 !< Integer grid number for equidistant cylindrical grid in grib1
  integer, public, parameter :: MERCATOR_GRID_ID_GRIB1 = 1 !< Integer grid number for Mercator grid in grib1
  integer, public, parameter :: LAMBERT_CONF_GRID_ID_GRIB1 = 3 !< Integer grid number for Lambert Conformal grid in grib1
  integer, public, parameter :: GAUSSIAN_GRID_ID_GRIB1 = 4 !< Integer grid number for Gaussian grid in grib1
  integer, public, parameter :: POLAR_STEREO_GRID_ID_GRIB1 = 5 !< Integer grid number for polar stereo grid in grib1
  integer, public, parameter :: ROT_EQUID_CYLIND_E_GRID_ID_GRIB1 = 203 !< Integer grid number for rotated equidistant cylindrical E-stagger grid
  integer, public, parameter :: ROT_EQUID_CYLIND_B_GRID_ID_GRIB1 = 205 !< Integer grid number for rotated equidistant cylindrical B-stagger grid

  integer, public, parameter :: EQUID_CYLIND_GRID_ID_GRIB2 = 0 !< Integer grid number for equidistant cylindrical grid in grib2
  integer, public, parameter :: ROT_EQUID_CYLIND_GRID_ID_GRIB2 = 1 !< Integer grid number for rotated equidistant cylindrical grid in grib2
  integer, public, parameter :: MERCATOR_GRID_ID_GRIB2 = 10 !< Integer grid number for Mercator grid in grib2
  integer, public, parameter :: POLAR_STEREO_GRID_ID_GRIB2 = 20 !< Integer grid number for polar stereo grid in grib2
  integer, public, parameter :: LAMBERT_CONF_GRID_ID_GRIB2 = 30 !< Integer grid number for Lambert conformal grid in grib2
  integer, public, parameter :: GAUSSIAN_GRID_ID_GRIB2 = 40 !< Integer grid number for Gaussian grid in grib2

  private
  public :: ip_grid
  !> Interface to gdswzd().  
  public :: gdswzd_interface
  public :: operator(==)

  !> Abstract grid that holds fields and methods common to all grids.
  !! ip_grid is meant to be subclassed when implementing a new grid.
  !!
  !! There are three methods that must be implemented:
  !! - init_grib1
  !! - init_grib2
  !! - gdswzd
  !!
  !! The init methods are responsible for setting up the grid
  !! using grib1/grib2 descriptors.
  !!
  !! Gdswzd performs transformations to/from Earth coordinates and grid coordinates.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  type, abstract :: ip_grid
     class(ip_grid_descriptor), allocatable :: descriptor !< Descriptor.
     
     integer :: im !< Number of x points
     integer :: jm !< Number of y points
     integer :: nm !< Total number of points

     !> @param Scanning mode.
     !! 0 if x first then y;
     !! 1 if y first then x;
     !! 3 if staggered diagonal like projection 203.
     integer :: nscan 
     integer :: kscan !< Mass/wind flag for staggered diagonal (0 if mass; 1 if wind)

     integer :: nscan_field_pos !< nscan for field_pos routine. Can be different than nscan due to differences in grib/grib2.
     
     integer :: iwrap !< x wraparound increment (0 if no wraparound).
     integer :: jwrap1 !< y wraparound lower pivot point (0 if no wraparound).
     integer :: jwrap2 !< y wraparound upper pivot point (0 if no wraparound).
     real :: rerth !< Radius of the Earth.
     real :: eccen_squared !< Eccentricity of the Earth squared (e^2)
   contains
     !> Initializer for grib1 input descriptor. @return N/A
     procedure(init_grib1_interface), deferred :: init_grib1
     !> Initializer for grib2 input descriptor. @return N/A
     procedure(init_grib2_interface), deferred :: init_grib2
     !> Coordinate transformations for the grid. @return N/A
     procedure(gdswzd_interface), deferred :: gdswzd
     !> Field position for a given grid point. @return Integer
     !> position in grib field to locate grid point.
     procedure :: field_pos
     !> Init subprogram. @return N/A
     generic :: init => init_grib1, init_grib2
  end type ip_grid

  abstract interface
     !> Interface to gdswzd().
     !>
     !> @param[in] self ???
     !> @param[in] iopt ???
     !> @param[in] npts ???
     !> @param[in] fill ???
     !> @param[in] xpts ???
     !> @param[in] ypts ???
     !> @param[in] rlon ???
     !> @param[in] rlat ???
     !> @param[in] nret ???
     !> @param[in] crot ???
     !> @param[in] srot ???
     !> @param[in] xlon ???
     !> @param[in] xlat ???
     !> @param[in] ylon ???
     !> @param[in] ylat ???
     !> @param[in] area ???
     !>
     !> @author Kyle Gerheiser @date July 2021
     subroutine gdswzd_interface(self, iopt, npts, fill, xpts, ypts, rlon, rlat, nret, crot, srot, &
          xlon, xlat, ylon, ylat, area)
       import
       class(ip_grid), intent(in) :: self
       INTEGER,          INTENT(IN   ) :: IOPT, NPTS
       INTEGER,          INTENT(  OUT) :: NRET
       !
       REAL,             INTENT(IN   ) :: FILL
       REAL,             INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
       REAL,             INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
       REAL, OPTIONAL,   INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
       REAL, OPTIONAL,   INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
       REAL, OPTIONAL,   INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
     end subroutine gdswzd_interface

     !> Init GRIB1 interface.
     !!
     !! @param[inout] self ???
     !! @param[in] g1_desc ???
     !!
     !! @author Kyle Gerheiser
     !! @date July 2021
     subroutine init_grib1_interface(self, g1_desc)
       import
       class(ip_grid), intent(inout) :: self
       type(grib1_descriptor), intent(in) :: g1_desc
     end subroutine init_grib1_interface

     !> Init GRIB2 interface.
     !!
     !! @param[inout] self ???
     !! @param[in] g2_desc ???
     !!
     !! @author Kyle Gerheiser
     !! @date July 2021
     subroutine init_grib2_interface(self, g2_desc)
       import
       class(ip_grid), intent(inout) :: self
       type(grib2_descriptor), intent(in) :: g2_desc
     end subroutine init_grib2_interface

  end interface

  !> Check equality.
  !> @author Kyle Gerheiser @date July 2021
  interface operator (==)
     module procedure is_same_grid
  end interface operator (==)
  

contains

  !> Compares two grids.
  !!
  !! @param[in] grid1 An ip_grid
  !! @param[in] grid2 Another ip_grid
  !!
  !! @return True if the grids are the same, false if not.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  logical function is_same_grid(grid1, grid2)
    class(ip_grid), intent(in) :: grid1, grid2
    is_same_grid = grid1%descriptor == grid2%descriptor
  end function is_same_grid

  !> Returns the field position for a given grid point.
  !!
  !! @param[in] self
  !! @param[in] i 
  !! @param[in] j
  !!
  !! @return Integer position in grib field to locate grid point.
  !!
  !! @author Mark Iredell, George Gayno, Kyle Gerheiser
  !! @date April 1996
  function field_pos(self, i, j)
    class(ip_grid), intent(in) :: self
    integer, intent(in) :: i, j
    integer :: field_pos

    integer :: ii, jj, im, jm
    integer :: iif, jjf, is1, iwrap
    integer :: jwrap1, jwrap2, kscan, nscan

    ! extract from navigation parameter array
    im=self%im
    jm=self%jm
    iwrap=self%iwrap
    jwrap1=self%jwrap1
    jwrap2=self%jwrap2
    nscan=self%nscan_field_pos
    kscan=self%kscan

    ! compute wraparounds in x and y if necessary and possible
    ii=i
    jj=j
    if(iwrap.gt.0) then
       ii=mod(i-1+iwrap,iwrap)+1
       if(j.lt.1.and.jwrap1.gt.0) then
          jj=jwrap1-j
          ii=mod(ii-1+iwrap/2,iwrap)+1
       elseif(j.gt.jm.and.jwrap2.gt.0) then
          jj=jwrap2-j
          ii=mod(ii-1+iwrap/2,iwrap)+1
       endif
    endif

    ! compute position for the appropriate scanning mode
    field_pos=0
    if(nscan.eq.0) then
       if(ii.ge.1.and.ii.le.im.and.jj.ge.1.and.jj.le.jm) field_pos=ii+(jj-1)*im
    elseif(nscan.eq.1) then
       if(ii.ge.1.and.ii.le.im.and.jj.ge.1.and.jj.le.jm) field_pos=jj+(ii-1)*jm
    elseif(nscan.eq.2) then
       is1=(jm+1-kscan)/2
       iif=jj+(ii-is1)
       jjf=jj-(ii-is1)+kscan
       if(iif.ge.1.and.iif.le.2*im-1.and.jjf.ge.1.and.jjf.le.jm) &
            field_pos=(iif+(jjf-1)*(2*im-1)+1-kscan)/2
    elseif(nscan.eq.3) then
       is1=(jm+1-kscan)/2
       iif=jj+(ii-is1)
       jjf=jj-(ii-is1)+kscan
       if(iif.ge.1.and.iif.le.2*im-1.and.jjf.ge.1.and.jjf.le.jm) field_pos=(iif+1)/2+(jjf-1)*im
    endif
  end function field_pos


end module ip_grid_mod

