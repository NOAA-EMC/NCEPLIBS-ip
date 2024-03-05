!> @file
!> @brief Determine earth radius and shape.
!> @author Kyle Gerheiser @date 2021-07-21

!> @brief Determine earth radius and shape.
!> @author Gayno, Kyle Gerheiser
module earth_radius_mod
    implicit none

    private
    public :: earth_radius

contains

    !> Determine earth radius and shape.
    !>
    !> Determine the radius and shape of the earth from
    !> the grib 2 grid definition template array - section 3.
    !>
    !> @param[in] igdtmpl integer (igdtlen) grid definition template
    !> array. Corresponds to the gfld%igdtmpl component of the ncep g2
    !> library gridmod data structure. For all map projections recognized
    !> by iplib, the entries use by this routine are:
    !> - 1 shape of earth, section 3, octet 15
    !> - 2 scale factor of spherical earth radius, octet 16
    !> - 3 scaled value of radius of spherical earth, octets 17-20
    !> - 4 scale factor of major axis of elliptical earth, octet 21
    !> - 5 scaled value of major axis of elliptical earth, octets 22-25
    !> - 6 scale factor of minor axis of elliptical earth, octet 26
    !> - 7 scaled value of minor axis of elliptical earth, octets 27-30
    !> @param[in] igdtlen integer number of elements of the grid
    !> definition template array. Corresponds to the gfld%igdtlen
    !> component of the ncep g2 library gridmod data structure.
    !> @param[out] radius real earth radius in meters. For ellipitical
    !> earths, this is the semi major axis. See "map projectsions - a
    !> working manual" by Snyder (1987) for details.
    !> @param[out] eccen_squared real earth eccentricity squared
    !>
    !> @author Gayno @date 2015-07-14
    subroutine earth_radius(igdtmpl, igdtlen, radius, eccen_squared)
        implicit none

        integer, intent(in) :: igdtlen
        integer, intent(in) :: igdtmpl(igdtlen)

        real, intent(out) :: eccen_squared
        real, intent(out) :: radius

        real                                  :: flat
        real                                  :: major_axis, minor_axis

        select case (igdtmpl(1))
        case (0)
            radius = 6367470.0
            eccen_squared = 0.0
        case (1)  ! USER SPECIFIED SPHERICAL
            radius = float(igdtmpl(3))/float(10**igdtmpl(2))
            eccen_squared = 0.0
        case (2)  ! IAU 1965
            radius = 6378160.0      ! SEMI MAJOR AXIS
            flat = 1.0/297.0      ! FLATTENING
            eccen_squared = (2.0*flat)-(flat**2)
        case (3)  ! USER SPECIFIED ELLIPTICAL (KM)
            major_axis = float(igdtmpl(5))/float(10**igdtmpl(4))
            major_axis = major_axis*1000.0
            minor_axis = float(igdtmpl(7))/float(10**igdtmpl(6))
            minor_axis = minor_axis*1000.0
            eccen_squared = 1.0-(minor_axis**2/major_axis**2)
            radius = major_axis
        case (4)  ! IAG-GRS80 MODEL
            radius = 6378137.0      ! SEMI MAJOR AXIS
            flat = 1.0/298.2572   ! FLATTENING
            eccen_squared = (2.0*flat)-(flat**2)
        case (5)  ! WGS84 DATUM
            radius = 6378137.0      ! SEMI MAJOR AXIS
            eccen_squared = 0.00669437999013
        case (6)
            radius = 6371229.0
            eccen_squared = 0.0
        case (7)  ! USER SPECIFIED ELLIPTICAL (M)
            major_axis = float(igdtmpl(5))/float(10**igdtmpl(4))
            minor_axis = float(igdtmpl(7))/float(10**igdtmpl(6))
            eccen_squared = 1.0-(minor_axis**2/major_axis**2)
            radius = major_axis
        case (8)
            radius = 6371200.0
            eccen_squared = 0.0
        case default
            radius = -9999.
            eccen_squared = -9999.
        end select
        !
        return
        !
    end subroutine earth_radius
end module earth_radius_mod
