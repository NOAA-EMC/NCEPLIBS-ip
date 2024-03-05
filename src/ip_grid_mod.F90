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

    integer,public,parameter :: equid_cylind_grid_id_grib1=0 !< Integer grid number for equidistant cylindrical grid in grib1
    integer,public,parameter :: mercator_grid_id_grib1=1 !< Integer grid number for Mercator grid in grib1
    integer,public,parameter :: lambert_conf_grid_id_grib1=3 !< Integer grid number for Lambert Conformal grid in grib1
    integer,public,parameter :: gaussian_grid_id_grib1=4 !< Integer grid number for Gaussian grid in grib1
    integer,public,parameter :: polar_stereo_grid_id_grib1=5 !< Integer grid number for polar stereo grid in grib1
    integer,public,parameter :: rot_equid_cylind_e_grid_id_grib1=203 !< Integer grid number for rotated equidistant cylindrical E-stagger grid
    integer,public,parameter :: rot_equid_cylind_b_grid_id_grib1=205 !< Integer grid number for rotated equidistant cylindrical B-stagger grid

    integer,public,parameter :: equid_cylind_grid_id_grib2=0 !< Integer grid number for equidistant cylindrical grid in grib2
    integer,public,parameter :: rot_equid_cylind_grid_id_grib2=1 !< Integer grid number for rotated equidistant cylindrical grid in grib2
    integer,public,parameter :: mercator_grid_id_grib2=10 !< Integer grid number for Mercator grid in grib2
    integer,public,parameter :: polar_stereo_grid_id_grib2=20 !< Integer grid number for polar stereo grid in grib2
    integer,public,parameter :: lambert_conf_grid_id_grib2=30 !< Integer grid number for Lambert conformal grid in grib2
    integer,public,parameter :: gaussian_grid_id_grib2=40 !< Integer grid number for Gaussian grid in grib2
    integer,public,parameter :: rot_equid_cylind_e_grid_id_grib2=32768 !< Integer grid number for rotated equidistant cylindrical E-stagger grid (grib2)
    integer,public,parameter :: rot_equid_cylind_b_grid_id_grib2=32769 !< Integer grid number for rotated equidistant cylindrical B-stagger grid (grib2)

    private
    public :: ip_grid
    public :: gdswzd_interface
    public :: operator(==)

    !> Abstract grid that holds fields and methods common to all grids.
  !! ip_grid is meant to be subclassed when implementing a new grid.
  !!
  !! There are three methods that must be implemented:
  !! - init_grib1()
  !! - init_grib2()
  !! - gdswzd()
  !!
  !! The init methods are responsible for setting up the grid
  !! using GRIB1/GRIB2 descriptors.
  !!
  !! gdswzd() performs transformations to and from Earth coordinates
  !! and grid coordinates.
  !!
  !! A good reference for all the map projection equations used by
  !! NCEPLIBS-ip can be found here: https://doi.org/10.3133/pp1395.
  !!
  !! @author Kyle Gerheiser @date July 2021
    type,abstract :: ip_grid
        class(ip_grid_descriptor),allocatable :: descriptor !< Descriptor.

        integer :: im !< Number of x points
        integer :: jm !< Number of y points
        integer :: nm !< Total number of points

        !> @param Scanning mode.
     !! - 0 if x first then y;
     !! - 1 if y first then x;
     !! - 3 if staggered diagonal like projection 203.
        integer :: nscan
        integer :: kscan !< Mass/wind flag for staggered diagonal (0 if mass; 1 if wind).

        integer :: nscan_field_pos !< nscan for field_pos routine. Can be different than nscan due to differences in grib/grib2.

        integer :: iwrap !< x wraparound increment (0 if no wraparound).
        integer :: jwrap1 !< y wraparound lower pivot point (0 if no wraparound).
        integer :: jwrap2 !< y wraparound upper pivot point (0 if no wraparound).
        real :: rerth !< Radius of the Earth.
        real :: eccen_squared !< Eccentricity of the Earth squared (e^2).
    contains
        !> Initializer for grib1 input descriptor. @return N/A
        procedure(init_grib1_interface),deferred :: init_grib1
        !> Initializer for grib2 input descriptor. @return N/A
        procedure(init_grib2_interface),deferred :: init_grib2
        !> Coordinate transformations for the grid. @return N/A
        procedure(gdswzd_interface),deferred :: gdswzd
        !> Field position for a given grid point. @return Integer
        !> position in grib field to locate grid point.
        procedure :: field_pos
        !> Init subprogram. @return N/A
        generic :: init=>init_grib1,init_grib2
    endtype ip_grid

    abstract interface

        !> @fn ip_grid_mod::gdswzd_interface::gdswzd_interface(self,
        !> iopt, npts, fill, xpts, ypts, rlon, rlat, nret, crot, srot,
        !> xlon, xlat, ylon, ylat, area)
        !> Interface to gdswzd().
        !>
        !> @param[in] self ip_grid_mod object.
        !> @param[in] iopt option flag
        !> - 1 to compute earth coords of selected grid coords
        !> - -1 to compute grid coords of selected earth coords
        !> @param[in] npts maximum number of coordinates
        !> @param[in] fill fill value to set invalid output data (must be
        !> impossible value; suggested value: -9999.)
        !> @param[inout] xpts (npts) grid x point coordinates if iopt>0
        !> @param[inout] ypts (npts) grid y point coordinates if iopt>0
        !> @param[inout] rlon (npts) earth longitudes in degrees e if iopt<0
        !> (acceptable range: -360. to 360.)
        !> @param[inout] rlat (npts) earth latitudes in degrees n if iopt<0
        !> (acceptable range: -90. to 90.)
        !> @param[out] nret number of valid points computed
        !> @param[out] crot optional (npts) clockwise vector rotation cosines
        !> @param[out] srot optional (npts) clockwise vector rotation sines
        !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
        !> @param[out] xlon optional (npts) dx/dlon in 1/degrees
        !> @param[out] xlat optional (npts) dx/dlat in 1/degrees
        !> @param[out] ylon optional (npts) dy/dlon in 1/degrees
        !> @param[out] ylat optional (npts) dy/dlat in 1/degrees
        !> @param[out] area optional (npts) area weights in m**2
        !> (proportional to the square of the map factor)
        !>
        !> @author Kyle Gerheiser @date July 2021
        subroutine gdswzd_interface(self,iopt,npts,fill,xpts,ypts,rlon,rlat,nret,crot,srot, &
                                    xlon,xlat,ylon,ylat,area)
            import
            class(ip_grid),intent(in) :: self
            integer,intent(in) :: iopt,npts
            integer,intent(out) :: nret
            !
            real,intent(in) :: fill
            real,intent(inout) :: rlon(npts),rlat(npts)
            real,intent(inout) :: xpts(npts),ypts(npts)
            real,optional,intent(out) :: crot(npts),srot(npts)
            real,optional,intent(out) :: xlon(npts),xlat(npts)
            real,optional,intent(out) :: ylon(npts),ylat(npts),area(npts)
        endsubroutine gdswzd_interface

        !> @fn ip_grid_mod::init_grib1_interface::init_grib1_interface(self, g1_desc)
        !> Init GRIB1 interface.
        !>
        !> @param[inout] self ip_grid_mod object.
        !> @param[in] g1_desc GRIB1 descriptor.
        !>
        !> @author Kyle Gerheiser
        !> @date July 2021
        subroutine init_grib1_interface(self,g1_desc)
            import
            class(ip_grid),intent(inout) :: self
            type(grib1_descriptor),intent(in) :: g1_desc
        endsubroutine init_grib1_interface

        !> @fn ip_grid_mod::init_grib2_interface::init_grib2_interface(self, g2_desc)
        !> Init GRIB2 interface.
        !>
        !> @param[inout] self ip_grid_mod object.
        !> @param[in] g2_desc GRIB2 descriptor.
        !>
        !> @author Kyle Gerheiser
        !> @date July 2021
        subroutine init_grib2_interface(self,g2_desc)
            import
            class(ip_grid),intent(inout) :: self
            type(grib2_descriptor),intent(in) :: g2_desc
        endsubroutine init_grib2_interface

    endinterface

    !> Check equality.
    !> @author Kyle Gerheiser @date July 2021
    interface operator(==)
        module procedure is_same_grid
    endinterface operator(==)

contains

    !> Compares two grids.
    !>
    !> @param[in] grid1 An ip_grid
    !> @param[in] grid2 Another ip_grid
    !>
    !> @return True if the grids are the same, false if not.
    !>
    !> @author Kyle Gerheiser
    !> @date July 2021
    logical function is_same_grid(grid1,grid2)
        class(ip_grid),intent(in) :: grid1,grid2
        is_same_grid=grid1%descriptor.eq.grid2%descriptor
    endfunction is_same_grid

    !> Returns the field position for a given grid point.
    !>
    !> @param[in] self
    !> @param[in] i
    !> @param[in] j
    !>
    !> @return Integer position in grib field to locate grid point.
    !>
    !> @author Mark Iredell, George Gayno, Kyle Gerheiser
    !> @date April 1996
    function field_pos(self,i,j)
        class(ip_grid),intent(in) :: self
        integer,intent(in) :: i,j
        integer :: field_pos

        integer :: ii,jj,im,jm
        integer :: iif,jjf,is1,iwrap
        integer :: jwrap1,jwrap2,kscan,nscan

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
    endfunction field_pos

endmodule ip_grid_mod

