!> @file
!> @brief Interpolate gridded data to a series of station points.
!>
!> @author Kyle Gerheiser @date 7/21/21

!> @brief Interpolate gridded data to a series of station points.
!>
!> @author Kyle Gerheiser @date 7/21/21
module ip_station_points_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  implicit none

  ! Not really a grid
  private
  public :: ip_station_points_grid

  type, extends(ip_grid) :: ip_station_points_grid
   contains
     !> Initializes a gaussian grid given a grib1_descriptor object.
     !> @return N/A @copydoc ip_station_points_grid_mod::init_grib1
     procedure :: init_grib1 
     !> Initializes a gaussian grid given a grib2_descriptor object.
     !> @return N/A @copydoc ip_station_points_grid_mod::init_grib2
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
     !> for IP Station Point grids.
     !> @return N/A @copydoc ip_station_points_grid_mod::gdswzd_station_points
     procedure :: gdswzd => GDSWZD_station_points 
  end type ip_station_points_grid

contains

  !> Initializes an IP Station grid given a grib1_descriptor object.
  !>
  !> @param[inout] self The grid to initialize
  !> @param[in] g1_desc A grib1_descriptor
  !>
  !> @author Iredell @date 96-04-10
  subroutine init_grib1(self, g1_desc)
    class(ip_station_points_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc
  end subroutine init_grib1

  !> Initializes an IP Station grid given a grib2_descriptor object.
  !>
  !> @param[inout] self The grid to initialize
  !> @param[in] g2_desc A grib2_descriptor
  !>
  !> @author Iredell @date 96-04-10
  subroutine init_grib2(self, g2_desc)
    class(ip_station_points_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc
  end subroutine init_grib2

  !> Interpolate gridded data to a series of station points.
  !>
  !> @param[in] self The grid.
  !> @param[in] IOPT must be minus 1 (return grid coordinates for selected earth coordinates).
  !> @param[in] NPTS  Maximum number of points.
  !> @param[in] FILL Fill value to set invalid output data.
  !> Must be impossible value; suggested value: -9999.
  !> @param[inout] XPTS X point coordinates. Always output.
  !> @param[inout] YPTS Y point coordinates. Always output.
  !> @param[inout] RLON Point longitudes. Always input.
  !> @param[inout] RLAT Point latitudes. Always input.
  !> @param[out] NRET Number of valid points computed.
  !> @param[out] CROT Not used.
  !> @param[out] SROT Not used.
  !> @param[out] XLON Not used.
  !> @param[out] XLAT Not used.
  !> @param[out] YLON Not used.
  !> @param[out] YLAT Not used.
  !> @param[out] AREA Not used.
  !>
  !> @author Kyle Gerheiser @date 7/21/21
  !> @author Eric Engle @date 5/4/23
  SUBROUTINE GDSWZD_station_points(self,IOPT,NPTS, &
       FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    class(ip_station_points_grid), intent(in) :: self
    INTEGER,                 INTENT(IN   ) :: IOPT, NPTS
    INTEGER,                 INTENT(  OUT) :: NRET
    !
    REAL,                    INTENT(IN   ) :: FILL
    REAL,                    INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,                    INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)

    ! This is all that needs to be done for GDSWZD for station points.
    NRET = NPTS

  end subroutine GDSWZD_STATION_POINTS

end module ip_station_points_grid_mod
