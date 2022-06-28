!> @file
!> @brief ???
!>
!> @author Iredell @date 96-04-10

!> @brief ???
!>
!> @author Iredell @date 96-04-10
module ip_station_points_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  implicit none

  ! Not really a grid
  private
  public :: ip_station_points_grid

  type, extends(ip_grid) :: ip_station_points_grid
   contains
     !> Initializes a gaussian grid given a grib1_descriptor object. @return N/A
     procedure :: init_grib1 
     !> Initializes a gaussian grid given a grib2_descriptor object. @return N/A
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
     !> for IP Station Point grids. @return N/A
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

  !> ???
  !>
  !> @param[in] self ???
  !> @param[in] IOPT ???
  !> @param[in] NPTS ???
  !> @param[in] FILL ???
  !> @param[inout] XPTS ???
  !> @param[inout] YPTS ???
  !> @param[inout] RLON ???
  !> @param[inout] RLAT ???
  !> @param[out] NRET ???
  !> @param[out] CROT ???
  !> @param[out] SROT ???
  !> @param[out] XLON ???
  !> @param[out] XLAT ???
  !> @param[out] YLON ???
  !> @param[out] YLAT ???
  !> @param[out] AREA ???
  !>
  !> @author Iredell @date 96-04-10
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

  end subroutine GDSWZD_STATION_POINTS

end module ip_station_points_grid_mod
