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

  function init_grid_generic(grid_desc) result(grid)
    class(ip_grid_descriptor), intent(in) :: grid_desc
    class(ip_grid), allocatable :: grid

    select type(grid_desc)
    type is(grib1_descriptor)
       grid = init_grid_grib1(grid_desc)
    type is(grib2_descriptor)
       grid = init_grid_grib2(grid_desc)
    end select
  end function init_grid_generic
  

  function init_grid_grib1(g1_desc) result(grid)
    type(grib1_descriptor), intent(in) :: g1_desc
    class(ip_grid), allocatable :: grid
    
    select case(g1_desc%grid_num)
    case(:-1)
       allocate(ip_station_points_grid::grid)
    case(EQUID_CYLIND_GRID_ID_GRIB1)
       allocate(ip_equid_cylind_grid::grid)
    case(MERCATOR_GRID_ID_GRIB1)
       allocate(ip_mercator_grid::grid)
    case(LAMBERT_CONF_GRID_ID_GRIB1)
       allocate(ip_lambert_conf_grid::grid)
    case(GAUSSIAN_GRID_ID_GRIB1)
       allocate(ip_gaussian_grid::grid)
    case(POLAR_STEREO_GRID_ID_GRIB1)
       allocate(ip_polar_stereo_grid::grid)
    case(ROT_EQUID_CYLIND_E_GRID_ID_GRIB1)
       allocate(ip_rot_equid_cylind_egrid::grid)
    case(ROT_EQUID_CYLIND_B_GRID_ID_GRIB1)
       allocate(ip_rot_equid_cylind_grid::grid)
    end select

    call grid%init(g1_desc)
    allocate(grid%descriptor, source = g1_desc)
  end function init_grid_grib1


  function init_grid_grib2(g2_desc) result(grid)
    type(grib2_descriptor), intent(in) :: g2_desc
    class(ip_grid), allocatable :: grid

    integer :: i_offset_odd, i_offset_even

    select case(g2_desc%grid_num)
    case(:-1)
       allocate(ip_station_points_grid::grid)
    case(EQUID_CYLIND_GRID_ID_GRIB2)
       allocate(ip_equid_cylind_grid::grid)
    case(ROT_EQUID_CYLIND_GRID_ID_GRIB2)
       i_offset_odd = mod(g2_desc%gdt_tmpl(19) / 8, 2)
       i_offset_even = mod(g2_desc%gdt_tmpl(19) / 4, 2)
       if (i_offset_odd /= i_offset_even) then
          allocate(ip_rot_equid_cylind_egrid::grid)
       else
          allocate(ip_rot_equid_cylind_grid::grid)
       end if
    case(MERCATOR_GRID_ID_GRIB2)
       allocate(ip_mercator_grid::grid)
    case(POLAR_STEREO_GRID_ID_GRIB2)
       allocate(ip_polar_stereo_grid::grid)
    case(LAMBERT_CONF_GRID_ID_GRIB2)
       allocate(ip_lambert_conf_grid::grid)
    case(GAUSSIAN_GRID_ID_GRIB2)
       allocate(ip_gaussian_grid::grid)
    case default
       print *, "gdt_num: ", g2_desc%gdt_num, " not recognized"
       error stop
    end select

    call grid%init(g2_desc)
    allocate(grid%descriptor, source = g2_desc)
  end function init_grid_grib2
  
end module ip_grid_factory_mod
