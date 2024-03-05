!> @file
!! @brief Routines for creating an ip_grid given a Grib descriptor.
!!
!! @author Mark Iredell, George Gayno, Kyle Gerheiser
!! @date July 2021

!> @brief Routines for creating an ip_grid given a Grib descriptor.
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
!! @date July 2021
module ip_grid_factory_mod
    use ip_grid_descriptor_mod
    use ip_grids_mod
    use ip_grid_mod
    implicit none

    private
    public :: init_grid

    interface init_grid
        module procedure init_grid_generic
    end interface init_grid

contains

    !> Initializes a polymorphic ip_grid object from an ip_grid_descriptor.
  !!
  !! @param[out] grid Grid to initialize
  !! @param[in] grid_desc Grid descriptor created from a grib1/grib2 template.
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
    subroutine init_grid_generic(grid, grid_desc)
        class(ip_grid_descriptor), intent(in) :: grid_desc
        class(ip_grid), allocatable, intent(out) :: grid

        select type (grid_desc)
        type is (grib1_descriptor)
            call init_grid_grib1(grid, grid_desc)
        type is (grib2_descriptor)
            call init_grid_grib2(grid, grid_desc)
        end select
    end subroutine init_grid_generic

    !> Initializes a polymorphic ip_grid from a grib1_descriptor.
  !! The concrete grid type is chosen based on the grid number in the descriptor.
  !!
  !! @param[out] grid Grid to initialize
  !! @param[in] g1_desc
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
    subroutine init_grid_grib1(grid, g1_desc)
        type(grib1_descriptor), intent(in) :: g1_desc
        class(ip_grid), allocatable, intent(out) :: grid

        select case (g1_desc%grid_num)
        case (:-1)
            allocate (ip_station_points_grid::grid)
        case (equid_cylind_grid_id_grib1)
            allocate (ip_equid_cylind_grid::grid)
        case (mercator_grid_id_grib1)
            allocate (ip_mercator_grid::grid)
        case (lambert_conf_grid_id_grib1)
            allocate (ip_lambert_conf_grid::grid)
        case (gaussian_grid_id_grib1)
            allocate (ip_gaussian_grid::grid)
        case (polar_stereo_grid_id_grib1)
            allocate (ip_polar_stereo_grid::grid)
        case (rot_equid_cylind_e_grid_id_grib1)
            allocate (ip_rot_equid_cylind_egrid::grid)
        case (rot_equid_cylind_b_grid_id_grib1)
            allocate (ip_rot_equid_cylind_grid::grid)
        end select

        call grid%init(g1_desc)
        allocate (grid%descriptor, source=g1_desc)
    end subroutine init_grid_grib1

    !> Initializes a polymorphic ip_grid from a grib2_descriptor.
  !! The concrete grid type is chosen based on the grid number in the descriptor.
  !!
  !! @param[out] grid Grid to initialize
  !! @param[in] g2_desc Grib2 descriptor
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
    subroutine init_grid_grib2(grid, g2_desc)
        type(grib2_descriptor), intent(in) :: g2_desc
        class(ip_grid), allocatable, intent(out) :: grid

        integer :: i_offset_odd, i_offset_even

        select case (g2_desc%grid_num)
        case (:-1)
            allocate (ip_station_points_grid::grid)
        case (equid_cylind_grid_id_grib2)
            allocate (ip_equid_cylind_grid::grid)
        case (rot_equid_cylind_grid_id_grib2)
            i_offset_odd = mod(g2_desc%gdt_tmpl(19)/8, 2)
            i_offset_even = mod(g2_desc%gdt_tmpl(19)/4, 2)
            if (i_offset_odd .ne. i_offset_even) then
                allocate (ip_rot_equid_cylind_egrid::grid)
            else
                allocate (ip_rot_equid_cylind_grid::grid)
            end if
        case (mercator_grid_id_grib2)
            allocate (ip_mercator_grid::grid)
        case (polar_stereo_grid_id_grib2)
            allocate (ip_polar_stereo_grid::grid)
        case (lambert_conf_grid_id_grib2)
            allocate (ip_lambert_conf_grid::grid)
        case (gaussian_grid_id_grib2)
            allocate (ip_gaussian_grid::grid)
        case (rot_equid_cylind_e_grid_id_grib2)
            allocate (ip_rot_equid_cylind_egrid::grid)
        case (rot_equid_cylind_b_grid_id_grib2)
            allocate (ip_rot_equid_cylind_grid::grid)
        case default
            print *, "gdt_num: ", g2_desc%gdt_num, " not recognized"
            error stop
        end select

        call grid%init(g2_desc)
        allocate (grid%descriptor, source=g2_desc)
    end subroutine init_grid_grib2

end module ip_grid_factory_mod
