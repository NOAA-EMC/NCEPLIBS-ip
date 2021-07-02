module ip_grid_mod
  use ip_grid_descriptor_mod
  implicit none

  integer, public, parameter :: EQUID_CYLIND_GRID_ID_GRIB1 = 0
  integer, public, parameter :: MERCATOR_GRID_ID_GRIB1 = 1
  integer, public, parameter :: LAMBERT_CONF_GRID_ID_GRIB1 = 3
  integer, public, parameter :: GAUSSIAN_GRID_ID_GRIB1 = 4
  integer, public, parameter :: POLAR_STEREO_GRID_ID_GRIB1 = 5
  integer, public, parameter :: ROT_EQUID_CYLIND_E_GRID_ID_GRIB1 = 203
  integer, public, parameter :: ROT_EQUID_CYLIND_B_GRID_ID_GRIB1 = 205

  integer, public, parameter :: EQUID_CYLIND_GRID_ID_GRIB2 = 0
  integer, public, parameter :: ROT_EQUID_CYLIND_GRID_ID_GRIB2 = 1
  integer, public, parameter :: MERCATOR_GRID_ID_GRIB2 = 10
  integer, public, parameter :: POLAR_STEREO_GRID_ID_GRIB2 = 20
  integer, public, parameter :: LAMBERT_CONF_GRID_ID_GRIB2 = 30
  integer, public, parameter :: GAUSSIAN_GRID_ID_GRIB2 = 40

  private
  public :: ip_grid, gdswzd_interface, operator(==)

  type, abstract :: ip_grid
     class(ip_grid_descriptor), allocatable :: descriptor

     integer :: im, jm, nm
     integer :: nscan, kscan, nscan_field_pos
     integer :: iwrap, jwrap1, jwrap2
     
     real :: rerth, eccen_squared
   contains
     procedure(init_grib1_interface), deferred :: init_grib1
     procedure(init_grib2_interface), deferred :: init_grib2
     procedure(gdswzd_interface), deferred :: gdswzd
     procedure :: field_pos
     generic :: init => init_grib1, init_grib2
  end type ip_grid

  abstract interface
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

     subroutine init_grib1_interface(self, g1_desc)
       import
       class(ip_grid), intent(inout) :: self
       type(grib1_descriptor), intent(in) :: g1_desc
     end subroutine init_grib1_interface

     subroutine init_grib2_interface(self, g2_desc)
       import
       class(ip_grid), intent(inout) :: self
       type(grib2_descriptor), intent(in) :: g2_desc
     end subroutine init_grib2_interface

  end interface

  interface operator (==)
     module procedure is_same_grid
  end interface operator (==)
  

contains

  logical function is_same_grid(grid1, grid2)
    class(ip_grid), intent(in) :: grid1, grid2
    is_same_grid = grid1%descriptor == grid2%descriptor
  end function is_same_grid
  
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

