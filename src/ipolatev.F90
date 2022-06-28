!> @file
!> @brief Top-level driver for vector interpolation routine ipolates.
!> @author Mark Iredell, Kyle Gerheiser
!> @date July 2021

!> @brief Top-level driver for vector interpolation interpolation
!> routine ipolates.  Ipolatev is overloaded with interfaces for grib1
!> and grib2 descriptors
!>
!> @author George Gayno, Mark Iredell, Kyle Gerheiser
module ipolatev_mod
  use ip_interpolators_mod
  use ip_grid_factory_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod

  implicit none

  private
  public :: ipolatev, ipolatev_grib2, ipolatev_grib1_single_field, ipolatev_grib1

  interface ipolatev
     module procedure ipolatev_grib1
     module procedure ipolatev_grib1_single_field
     module procedure ipolatev_grib2
  end interface ipolatev
  
contains

  !> Interpolates scalar fields between grids given ip_grid objects.
  !>
  !> Calls the specific interpolation routines on the generic ip_grids
  !> created from a GRIB1 or GRIB2 descriptor.
  !>
  !> @param[in] ip Interpolation method.
  !> @param[in] ipopt Interpolation options.
  !> @param[in] grid_in Input grid.
  !> @param[in] grid_out Output grid object created.
  !> @param[in] mi Skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo Skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1.
  !> @param[in] km Number of fields to interpolate.
  !> @param[in] ibi Input bitmap flags.
  !> @param[in] li Input bitmaps (if respective ibi(k)=1).
  !> @param[in] ui Input u-component fields to interpolate.
  !> @param[in] vi Input v-component fields to interpolate.
  !> @param[out] no Number of output points (only if kgdso(1)<0).
  !> @param[out] rlat Output latitudes in degrees (if kgdso(1)<0).
  !> @param[out] rlon Output longitudes in degrees (if kgdso(1)<0).
  !> @param[inout] crot Vector rotation cosines (if igdtnumo>=0).
  !> @param[inout] srot Vector rotation sines (if igdtnumo>=0).
  !> @param[out] ibo Output bitmap flags.
  !> @param[out] lo Output bitmaps (always output).
  !> @param[out] uo Output u-component fields interpolated.
  !> @param[out] vo Output v-component fields interpolated.
  !> @param[out] iret Return code.
  !> - 0 Successful interpolation.
  !> - 1 Unrecognized interpolation method.
  !> - 2 Unrecognized input grid or no grid overlap.
  !> - 3 Unrecognized output grid.
  !> - 1x Invalid bicubic method parameters.
  !> - 3x Invalid budget method parameters.
  !> - 4x Invalid spectral method parameters.
  !> @date July 2021
  !> @author Kyle Gerheiser
  SUBROUTINE ipolatev_grid(IP,IPOPT,grid_in,grid_out, &
       MI,MO,KM,IBI,LI,UI,VI, &
       NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,               INTENT(IN   ) :: IP, IPOPT(20), IBI(KM)
    INTEGER,               INTENT(IN   ) :: KM, MI, MO
    INTEGER,               INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,             INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,             INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,                  INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,                  INTENT(INOUT) :: CROT(MO),SROT(MO)
    REAL,                  INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)

    select case(ip)
    case(BILINEAR_INTERP_ID)
       CALL interpolate_bilinear(IPOPT,grid_in,grid_out, &
            MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    case(BICUBIC_INTERP_ID)
       CALL interpolate_bicubic(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    case(NEIGHBOR_INTERP_ID)
       CALL interpolate_neighbor(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    case(BUDGET_INTERP_ID)
       CALL interpolate_budget(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    case(SPECTRAL_INTERP_ID)
       CALL interpolate_spectral(IPOPT,grid_in,grid_out, &
            MI,MO,KM,IBI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    case(NEIGHBOR_BUDGET_INTERP_ID)
       CALL interpolate_neighbor_budget(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,UI,VI,&
            NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)
    case default
       print *, "unrecognized interpolation option: ", ip
       error stop
       ! IF(IGDTNUMO.GE.0) NO=0
       ! DO K=1,KM
       !    IBO(K)=1
       !    DO N=1,NO
       !       LO(N,K)=.FALSE.
       !       UO(N,K)=0.
       !       VO(N,K)=0.
       !    ENDDO
       ! ENDDO
       ! IRET=1
    end select

  end subroutine ipolatev_grid
  
  !> This subprogram interpolates vector fields from any grid to any
  !> grid given a grib2 descriptor.
  !>
  !> This is a wrapper for ipolates_grid which converts a grib1
  !> descriptor into an ip_grid_descriptor, which is used to create
  !> an ip_grid. Only horizontal interpolation is performed.
  !>
  !> The following interpolation methods are possible:
  !> - (ip=0) bilinear
  !> - (ip=1) bicubic
  !> - (ip=2) neighbor
  !> - (ip=3) budget
  !> - (ip=4) spectral
  !> - (ip=6) neighbor-budget
  !>
  !> Some of these methods have interpolation options and/or
  !> restrictions on the input or output grids, both of which
  !> are documented more fully in their respective subprograms.
  !>
  !> Input and output grids are defined by their grib 2 grid
  !> definition template as decoded by the ncep g2 library. The
  !> current code recognizes the following projections, where
  !> "igdtnumi/o" is the grib 2 grid defintion template number
  !> for the input and output grids, respectively:
  !> - (igdtnumi/o=00) equidistant cylindrical
  !> - (igdtnumi/o=01) rotated equidistant cylindrical. "e" and non-"e" staggered
  !> - (igdtnumi/o=10) mercator cylindrical
  !> - (igdtnumi/o=20) polar stereographic azimuthal
  !> - (igdtnumi/o=30) lambert conformal conical
  !> - (igdtnumi/o=40) gaussian cylindrical
  !>
  !> As an added bonus the number of output grid points
  !> and their latitudes and longitudes are also returned.
  !>
  !> On the other hand, data may be interpolated to a set of station
  !> points if "igdtnumo"<0 (or subtracted from 255 for the budget 
  !> option), in which case the number of points and
  !> their latitudes and longitudes must be input.
  !>
  !> Input bitmaps will be interpolated to output bitmaps.
  !> Output bitmaps will also be created when the output grid
  !> extends outside of the domain of the input grid.
  !>
  !> The output field is set to 0 where the output bitmap is off.
  !>
  !> @param[in] ip Interpolation method
  !> - ip=0 for bilinear
  !> - ip=1 for bicubic
  !> - ip=2 for neighbor;
  !> - ip=3 for budget;
  !> - ip=4 for spectral;
  !> - ip=6 for neighbor-budget
  !>
  !> @param[in] ipopt Interpolation options
  !> - ip=0: (No options)
  !> - ip=1: Constraint option
  !> - ip=2: (No options)
  !> - ip=3: Number in radius, radius weights, search radius
  !> - ip=4: Spectral shape, spectral truncation
  !> - ip=6: Number in radius, radius weights ...)
  !>
  !> @param[in] igdtnumi Grid definition template number for the input grid.
  !> Corresponds to the gfld%igdtnum component of the ncep g2 library
  !> gridmod data structure:
  !> - 00 - Equidistant Cylindrical
  !> - 01 - Rotated Equidistant cylindrical. "e" and non-"e" staggered
  !> - 10 - Mercator Cyclindrical
  !> - 20 - Polar Stereographic azimuthal
  !> - 30 - Lambert Conformal Conical
  !> - 40 - Gaussian Equidistant Cyclindrical
  !>
  !> @param[in] igdtmpli Grid definition template array input grid.
  !> Corresponds to the gfld%igdtmpl component of the NCEPLIBS-g2
  !> gridmod data structure
  !>
  !> Section 3 Info:
  !>
  !> All map projections:
  !> - (1): Shape of earth, octet 15.
  !> - (2): Scale factor of spherical earth radius, octet 16.
  !> - (3): Scaled value of radius of spherical earth, octets 17-20.
  !> - (4): Scale factor of major axis of elliptical earth, octet 21.
  !> - (5): Scaled value of major axis of elliptical earth, octets 22-25.
  !> - (6): Scale factor of minor axis of elliptical earth, octet 26.
  !> - (7): Scaled value of minor axis of elliptical earth, octets 27-30.
  !>
  !> Equidistant Cyclindrical:
  !> - (8):  Number of points along a parallel, octs 31-34.
  !> - (9):  Number of points along a meridian, octs 35-38.
  !> - (10): Basic angle of initial production domain, octets 39-42.
  !> - (11): Subdivisions of basic angle, octets 43-46.
  !> - (12): Latitude of first grid point, octets 47-50.
  !> - (13): Longitude of first grid point, octets 51-54.
  !> - (14): Resolution and component flags, octet 55.
  !> - (15): Latitude of last grid point, octets 56-59.
  !> - (16): Longitude of last grid point, octets 60-63.
  !> - (17): i-direction increment, octets 64-67.
  !> - (18): j-direction increment, octets 68-71.
  !> - (19): Scanning mode, octet 72.
  !>
  !> Mercator Cyclindrical:
  !> - (8):  Number of points along a parallel, octs 31-34.
  !> - (9):  Number of points along a meridian, octs 35-38.
  !> - (10): Latitude of first point, octets 39-42.
  !> - (11): Longitude of first point, octets 43-46.
  !> - (12): Resolution and component flags, octet 47.
  !> - (13): Tangent latitude, octets 48-51.
  !> - (14): Latitude of last point, octets 52-55.
  !> - (15): Longitude of last point, octets 56-59.
  !> - (16): Scanning mode flags, octet 60.
  !> - (17): Orientation of grid, octets 61-64.
  !> - (18): Longitudinal grid length, octets 65-68.
  !> - (19): Latitudinal grid length, octets 69-72.
  !>
  !> Lambert Conformal Conical:
  !> - (8):  Number of points along x-axis, octs 31-34.
  !> - (9):  Number of points along y-axis, octs 35-38.
  !> - (10): Latitude of first point, octets 39-42.
  !> - (11): Longitude of first point, octets 43-46.
  !> - (12): Resolution of component flag, octet 47.
  !> - (13): Latitude where grid lengths specified,octets 48-51.
  !> - (14): Longitude of meridian that is parallel to y-axis, octets 52-55.
  !> - (15): x-direction grid length, octets 56-59.
  !> - (16): y-direction grid length, octets 60-63.
  !> - (17): Projection center flag, octet 64.
  !> - (18): Scanning mode, octet 65.
  !> - (19): First tangent latitude from pole, octets 66-69.
  !> - (20): Second tangent latitude from pole, octets 70-73.
  !> - (21): Latitude of south pole of projection, octets 74-77.
  !> - (22): Longitude of south pole of projection, octets 78-81.
  !>
  !> Gaussian Cylindrical:
  !> - (8):  Number of points along a parallel, octs 31-34.
  !> - (9):  Number of points along a meridian, octs 35-38.
  !> - (10): Basic angle of initial production domain, octets 39-42.
  !> - (11): Subdivisions of basic angle, octets 43-46.
  !> - (12): Latitude of first grid point, octets 47-50.
  !> - (13): Longitude of first grid point, octets 51-54.
  !> - (14): Resolution and component flags, octet 55.
  !> - (15): Latitude of last grid point, octets 56-59.
  !> - (16): Longitude of last grid point, octets 60-63.
  !> - (17): i-direction increment, octets 64-67.
  !> - (18): Number of parallels between pole and equator, octets 68-71.
  !> - (19): Scanning mode, octet 72.
  !>
  !> Polar Stereographic Azimuthal:
  !> - (8):  Number of points along x-axis, octets 31-34.
  !> - (9):  Number of points along y-axis, octets 35-38.
  !> - (10): Latitude of first grid point, octets 39-42.
  !> - (11): Longitude of first grid point, octets 43-46.
  !> - (12): Resolution and component flags, octet 47.
  !> - (13): True latitude, octets 48-51.
  !> - (14): Orientation longitude, octets 52-55.
  !> - (15): x-direction grid length, octets 56-59.
  !> - (16): y-direction grid length, octets 60-63.
  !> - (17): Projection center flag, octet 64.
  !> - (18): Scanning mode flags, octet 65.
  !>
  !> Rotated Equidistant Cyclindrical:
  !> - (8):  Number of points along a parallel, octs 31-34.
  !> - (9):  Number of points along a meridian, octs 35-38.
  !> - (10): Basic angle of initial production domain, octets 39-42.
  !> - (11): Subdivisions of basic angle, octets 43-46.
  !> - (12): Latitude of first grid point, octets 47-50.
  !> - (13): Longitude of first grid point, octets 51-54.
  !> - (14): Resolution and component flags, octet 55.
  !> - (15): Latitude of last grid point, octets 56-59.
  !> - (16): Longitude of last grid point, octets 60-63.
  !> - (17): i-direction increment, octets 64-67.
  !> - (18): j-direction increment, octets 68-71.
  !> - (19): Scanning mode, octet 72.
  !> - (20): Latitude of southern pole of projection, octets 73-76.
  !> - (21): Longitude of southern pole of projection, octets 77-80.
  !> - (22): Angle of rotation of projection, octs 81-84.
  !>
  !> @param[in] igdtleni Number of elements of the grid definition
  !> template array for the input grid. Corresponds to the gfld%igdtlen
  !> component of the ncep g2 library gridmod data structure.
  !>
  !> @param[in] igdtnumo Grid definition template number for the output grid.
  !> Corresponds to the gfld%igdtnum component of the
  !> ncep g2 library gridmod data structure.
  !> See "igdtnumi" for specific template definitions.
  !> Note: igdtnumo<0 means interpolate to random station points.
  !>
  !> @param[in] igdtmplo Grid definition template array for the output grid.
  !> Corresponds to the gfld%igdtmpl component of the ncep g2 library
  !> gridmod data structure.
  !> See "igdtmpli" for definition of array elements.
  !>
  !> @param[in] igdtleno Number of elements of the grid definition
  !> template array for the output grid. Corresponds to the gfld%igdtlen
  !> component of the ncep g2 library gridmod data structure.
  !>
  !> @param[in] mi    Skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo    Skip number between output grid fields if km>1
  !> or dimension of output grid fields if km=1.
  !> @param[in] km    Number of fields to interpolate.
  !> @param[in] ibi   Input bitmap flags.
  !> @param[in] li    Input bitmaps (if respective ibi(k)=1).
  !> @param[in] ui    Input u-component fields to interpolate.
  !> @param[in] vi    Input v-component fields to interpolate.
  !> @param[out] no Number of output points (only if kgdso(1)<0).
  !> @param[inout] rlat Output latitudes in degrees (if kgdso(1)<0).
  !> @param[inout] rlon Output longitudes in degrees (if kgdso(1)<0).
  !> @param[inout] crot Vector rotation cosines (if igdtnumo>=0).
  !> @param[inout] srot Vector rotation sines (if igdtnumo>=0).
  !> @param[out] ibo Output bitmap flags.
  !> @param[out] lo  Output bitmaps (always output).
  !> @param[out] uo  Output u-component fields interpolated.
  !> @param[out] vo  Output v-component fields interpolated.
  !> @param[out] iret Return code.
  !> - 0 Successful interpolation.
  !> - 1 Unrecognized interpolation method.
  !> - 2 Unrecognized input grid or no grid overlap.
  !> - 3 Unrecognized output grid.
  !> - 1x Invalid bicubic method parameters.
  !> - 3x Invalid budget method parameters.
  !> - 4x Invalid spectral method parameters.
  !>
  !> @note Examples demonstrating relative cpu costs.
  !> This example is interpolating 12 levels of winds
  !> from the 360 x 181 global grid (ncep grid 3)
  !> to the 93 x 68 hawaiian mercator grid (ncep grid 204).
  !>
  !> The example times are for the c90. As a reference, the cp time
  !> for unpacking the global 12 wind fields is 0.07 seconds.
  !>
  !>   METHOD  | IP| IPOPT       |  CP SECONDS
  !>   --------| --|-------------| ----------
  !>   BILINEAR| 0 |             |   0.05
  !>   BICUBIC | 1 | 0           |   0.16
  !>   BICUBIC | 1 | 1           |   0.17
  !>   NEIGHBOR| 2 |             |   0.02
  !>   BUDGET  | 3 | -1,-1       |   0.94
  !>   SPECTRAL| 4 | 0,40        |   0.31
  !>   SPECTRAL| 4 | 1,40        |   0.33
  !>   SPECTRAL| 4 | 0,-1        |   0.59
  !>   N-BUDGET| 6 | -1,-1       |   0.31
  !>
  !>   The spectral interpolation is fast for the mercator grid.
  !>   However, for some grids the spectral interpolation is slow.
  !>
  !>   The following example is interpolating 12 levels of winds
  !>   from the 360 x 181 global grid (ncep grid 3)
  !>   to the 93 x 65 conus lambert conformal grid (ncep grid 211).
  !>
  !>   METHOD  | IP| IPOPT        |CP SECONDS
  !>   --------| --| -------------|----------
  !>   BILINEAR| 0 |              | 0.05
  !>   BICUBIC | 1 | 0            | 0.15
  !>   BICUBIC | 1 | 1            | 0.16
  !>   NEIGHBOR| 2 |              | 0.02
  !>   BUDGET  | 3 | -1,-1        | 0.92
  !>   SPECTRAL| 4 | 0,40         | 4.51
  !>   SPECTRAL| 4 | 1,40         | 5.77
  !>   SPECTRAL| 4 | 0,-1         | 12.60
  !>   N-BUDGET| 6 | -1,-1        | 0.33
  !>
  !> @date July 2021
  !> @author Kyle Gerheiser
  subroutine ipolatev_grib2(ip,ipopt,igdtnumi,igdtmpli,igdtleni, &
       igdtnumo,igdtmplo,igdtleno, &
       mi,mo,km,ibi,li,ui,vi, &
       no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret) bind(c)

    INTEGER,               INTENT(IN   ) :: IP, IPOPT(20), IBI(KM)
    INTEGER,               INTENT(IN   ) :: KM, MI, MO
    INTEGER,               INTENT(IN   ) :: IGDTNUMI, IGDTLENI
    INTEGER,               INTENT(IN   ) :: IGDTMPLI(IGDTLENI)
    INTEGER,               INTENT(IN   ) :: IGDTNUMO, IGDTLENO
    INTEGER,               INTENT(IN   ) :: IGDTMPLO(IGDTLENO)
    INTEGER,               INTENT(  OUT) :: IBO(KM), IRET, NO
    !
    LOGICAL*1,             INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,             INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,                  INTENT(IN   ) :: UI(MI,KM),VI(MI,KM)
    REAL,                  INTENT(INOUT) :: CROT(MO),SROT(MO)
    REAL,                  INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT) :: UO(MO,KM),VO(MO,KM)
    !
    INTEGER                              :: K, N

    type(grib2_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(igdtnumi, igdtleni, igdtmpli)
    desc_out = init_descriptor(igdtnumo, igdtleno, igdtmplo)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)

    CALL ipolatev_grid(ip,IPOPT,grid_in,grid_out, &
         MI,MO,KM,IBI,LI,UI,VI,&
         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)

  end subroutine ipolatev_grib2

  !> @brief This subprogram interpolates vector field from any grid
  !> to any grid given a grib1 Grid Descriptor Section.
  !>
  !> @details Only horizontal interpolation is performed.
  !> The following interpolation methods are possible:
  !> - (ip=0) bilinear
  !> - (ip=1) bicubic
  !> - (ip=2) neighbor
  !> - (ip=3) budget
  !> - (ip=4) spectral
  !> - (ip=6) neighbor-budget
  !>
  !> Some of these methods have interpolation options and/or
  !> restrictions on the input or output grids, both of which
  !> are documented more fully in their respective subprograms.
  !>
  !> The grids are defined by their grid description sections
  !> (passed in integer form as decoded by subprogram w3fi63).
  !>
  !> The current code recognizes the following projections:
  !> - (kgds(1)=000) equidistant cylindrical
  !> - (kgds(1)=001) mercator cylindrical
  !> - (kgds(1)=003) lambert conformal conical
  !> - (kgds(1)=004) gaussian cylindrical
  !> - (kgds(1)=005) polar stereographic azimuthal
  !> - (kgds(1)=203) rotated equidistant cylindrical - e-stagger
  !> - (kgds(1)=205) rotated equidistant cylindrical - b-stagger
  !>
  !> Where kgds could be either input kgdsi or output kgdso.
  !>
  !> As an added bonus the number of output grid points
  !> and their latitudes and longitudes are also returned.
  !>
  !> On the other hand, the output can be a set of station points
  !> if kgdso(1)<0, in which case the number of points
  !> and their latitudes and longitudes must be input.
  !> for the budget approach, a subsection of the grid may
  !> be output by subtracting kgdso(1) from 255 and passing
  !> in the latitudes and longitudes of the points.
  !> Input bitmaps will be interpolated to output bitmaps.
  !>
  !> Output bitmaps will also be created when the output grid
  !> extends outside of the domain of the input grid.
  !> the output field is set to 0 where the output bitmap is off.
  !>
  !> @param ip Interpolation method
  !> - ip = BILINEAR_INTERP_ID = 0 for bilinear
  !> - ip = BICUBIC_INTERP_ID = 1 for bicubic
  !> - ip = NEIGHBOR_INTERP_ID = 2 for neighbor;
  !> - ip = BUDGET_INTERP_ID = 3 for budget;
  !> - ip = SPECTRAL_INTERP_ID = 4 for spectral;
  !> - ip = NEIGHBOR_BUDGET_INTERP_ID = 6 for neighbor-budget
  !> @param ipopt Interpolation options
  !> - ip=0 (bilinear): (No options)
  !> - ip=1 Cbicubic): constraint option
  !> - ip=2 (neighbor): (No options)
  !> - ip=3 (budget): Number in radius, radius weights, search radius
  !> - ip=4 (spectral): Spectral shape, spectral truncation
  !> - ip=6 (neighbor-budget): Number in radius, radius weights ...)
  !> @param[in] kgdsi Input gds parameters as decoded by w3fi63.
  !> @param[in] kgdso Output gds parameters.
  !> @param[in] mi Skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo Skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1.
  !> @param[in] km    Number of fields to interpolate.
  !> @param[in] ibi   Input bitmap flags.
  !> @param[in] li    Input bitmaps (if respective ibi(k)=1).
  !> @param[in] ui    Input u-component fields to interpolate.
  !> @param[in] vi    Input v-component fields to interpolate.
  !> @param[out] no Number of output points (only if kgdso(1)<0).
  !> @param[out] rlat Output latitudes in degrees (if kgdso(1)<0).
  !> @param[out] rlon Output longitudes in degrees (if kgdso(1)<0).
  !> @param[inout] crot Vector rotation cosines (if igdtnumo>=0).
  !> @param[inout] srot Vector rotation sines (if igdtnumo>=0).
  !> @param[out] ibo Output bitmap flags.
  !> @param[out] lo  Output bitmaps (always output).
  !> @param[out] uo  Output u-component fields interpolated.
  !> @param[out] vo  Output v-component fields interpolated.
  !> @param[out] iret Return code.
  !> - 0 Successful interpolation.
  !> - 1 Unrecognized interpolation method.
  !> - 2 Unrecognized input grid or no grid overlap.
  !> - 3 Unrecognized output grid.
  !> - 1x Invalid bicubic method parameters.
  !> - 3x Invalid budget method parameters.
  !> - 4x Invalid spectral method parameters.
  !>
  !> @note Examples demonstrating relative cpu costs.
  !> This example is interpolating 12 levels of winds
  !> from the 360 x 181 global grid (ncep grid 3)
  !> to the 93 x 68 hawaiian mercator grid (ncep grid 204).
  !>
  !> The example times are for the c90. As a reference, the cp time
  !> for unpacking the global 12 temperature fields is 0.07 seconds.
  !>
  !>   METHOD  | IP| IPOPT       |  CP SECONDS
  !>   --------| --|-------------| ----------
  !>   BILINEAR| 0 |             |   0.05
  !>   BICUBIC | 1 | 0           |   0.16
  !>   BICUBIC | 1 | 1           |   0.17
  !>   NEIGHBOR| 2 |             |   0.02
  !>   BUDGET  | 3 | -1,-1       |   0.94
  !>   SPECTRAL| 4 | 0,40        |   0.31
  !>   SPECTRAL| 4 | 1,40        |   0.33
  !>   SPECTRAL| 4 | 0,-1        |   0.59
  !>   N-BUDGET| 6 | -1,-1       |   0.31
  !>
  !>   The spectral interpolation is fast for the mercator grid.
  !>   However, for some grids the spectral interpolation is slow.
  !>
  !>   The following example is interpolating 12 levels of temperatures
  !>   from the 360 x 181 global grid (ncep grid 3)
  !>   to the 93 x 65 conus lambert conformal grid (ncep grid 211).
  !>
  !>   METHOD  | IP| IPOPT        |CP SECONDS
  !>   --------| --| -------------|----------
  !>   BILINEAR| 0 |              | 0.05
  !>   BICUBIC | 1 | 0            | 0.15
  !>   BICUBIC | 1 | 1            | 0.16
  !>   NEIGHBOR| 2 |              | 0.02
  !>   BUDGET  | 3 | -1,-1        | 0.92
  !>   SPECTRAL| 4 | 0,40         | 4.51
  !>   SPECTRAL| 4 | 1,40         | 5.77
  !>   SPECTRAL| 4 | 0,-1         | 12.60
  !>   N-BUDGET| 6 | -1,-1        | 0.33
  !>
  !> @date July 2021
  !> @author Kyle Gerheiser
  subroutine ipolatev_grib1(ip,ipopt,kgdsi,kgdso,mi,mo,km,ibi,li,ui,vi, &
       no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret) bind(c)
    IMPLICIT NONE
    !
    INTEGER,               INTENT(IN   ):: IP, IPOPT(20), IBI(KM)
    INTEGER,               INTENT(IN   ):: KM, MI, MO
    INTEGER,               INTENT(INOUT):: KGDSI(200), KGDSO(200)
    INTEGER,               INTENT(  OUT):: IBO(KM), IRET, NO
    !
    LOGICAL*1,             INTENT(IN   ):: LI(MI,KM)
    LOGICAL*1,             INTENT(  OUT):: LO(MO,KM)
    !
    REAL,                  INTENT(IN   ):: UI(MI,KM),VI(MI,KM)
    REAL,                  INTENT(INOUT):: CROT(MO),SROT(MO)
    REAL,                  INTENT(INOUT):: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT):: UO(MO,KM),VO(MO,KM)
    !
    INTEGER                             :: K, N, KGDSI11, KGDSO11

    type(grib1_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    IF(KGDSI(1).EQ.203) THEN
       KGDSI11=KGDSI(11)
       KGDSI(11)=IOR(KGDSI(11),256)
    ENDIF
    IF(KGDSO(1).EQ.203) THEN
       KGDSO11=KGDSO(11)
       KGDSO(11)=IOR(KGDSO(11),256)
    ENDIF

    desc_in = init_descriptor(kgdsi)
    desc_out = init_descriptor(kgdso)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)

    CALL ipolatev_grid(ip,IPOPT,grid_in,grid_out, &
         MI,MO,KM,IBI,LI,UI,VI,&
         NO,RLAT,RLON,CROT,SROT,IBO,LO,UO,VO,IRET)

    IF(KGDSI(1).EQ.203) THEN
       KGDSI(11)=KGDSI11
    ENDIF
    IF(KGDSO(1).EQ.203) THEN
       KGDSO(11)=KGDSO11
    ENDIF

  END SUBROUTINE IPOLATEV_grib1

  !> Special case of ipolatev_grib1 when interpolating a single field.
  !> Removes the km dimension of input arrays so scalars can be passed to ibi/ibo.
  !>
  !> @param ip Interpolation method
  !> - ip = BILINEAR_INTERP_ID = 0 for bilinear
  !> - ip = BICUBIC_INTERP_ID = 1 for bicubic
  !> - ip = NEIGHBOR_INTERP_ID = 2 for neighbor;
  !> - ip = BUDGET_INTERP_ID = 3 for budget;
  !> - ip = SPECTRAL_INTERP_ID = 4 for spectral;
  !> - ip = NEIGHBOR_BUDGET_INTERP_ID = 6 for neighbor-budget
  !> @param ipopt Interpolation options
  !> - ip=0 (bilinear): (No options)
  !> - ip=1 Cbicubic): constraint option
  !> - ip=2 (neighbor): (No options)
  !> - ip=3 (budget): Number in radius, radius weights, search radius
  !> - ip=4 (spectral): Spectral shape, spectral truncation
  !> - ip=6 (neighbor-budget): Number in radius, radius weights ...)
  !> @param[in] kgdsi Input gds parameters as decoded by w3fi63.
  !> @param[in] kgdso Output gds parameters.
  !> @param[in] mi Skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo Skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1.
  !> @param[in] km    Number of fields to interpolate.
  !> @param[in] ibi   Input bitmap flags.
  !> @param[in] li    Input bitmaps (if respective ibi(k)=1).
  !> @param[in] ui    Input u-component fields to interpolate.
  !> @param[in] vi    Input v-component fields to interpolate.
  !> @param[out] no Number of output points (only if kgdso(1)<0).
  !> @param[out] rlat Output latitudes in degrees (if kgdso(1)<0).
  !> @param[out] rlon Output longitudes in degrees (if kgdso(1)<0).
  !> @param[inout] crot Vector rotation cosines (if igdtnumo>=0).
  !> @param[inout] srot Vector rotation sines (if igdtnumo>=0).
  !> @param[out] ibo Output bitmap flags.
  !> @param[out] lo  Output bitmaps (always output).
  !> @param[out] uo  Output u-component fields interpolated.
  !> @param[out] vo  Output v-component fields interpolated.
  !> @param[out] iret Return code.
  !> - 0 Successful interpolation.
  !> - 1 Unrecognized interpolation method.
  !> - 2 Unrecognized input grid or no grid overlap.
  !> - 3 Unrecognized output grid.
  !> - 1x Invalid bicubic method parameters.
  !> - 3x Invalid budget method parameters.
  !> - 4x Invalid spectral method parameters.
  !>
  !> @date Jan 2022
  !> @author Kyle Gerheiser
  subroutine ipolatev_grib1_single_field(ip,ipopt,kgdsi,kgdso,mi,mo,km,ibi,li,ui,vi, &
       no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret) bind(c)
    IMPLICIT NONE
    !
    INTEGER,               INTENT(IN   ):: IP, IPOPT(20), IBI
    INTEGER,               INTENT(IN   ):: KM, MI, MO
    INTEGER,               INTENT(INOUT):: KGDSI(200), KGDSO(200)
    INTEGER,               INTENT(  OUT):: IBO, IRET, NO
    !
    LOGICAL*1,             INTENT(IN   ):: LI(MI)
    LOGICAL*1,             INTENT(  OUT):: LO(MO)
    !
    REAL,                  INTENT(IN   ):: UI(MI),VI(MI)
    REAL,                  INTENT(INOUT):: CROT(MO),SROT(MO)
    REAL,                  INTENT(INOUT):: RLAT(MO),RLON(MO)
    REAL,                  INTENT(  OUT):: UO(MO),VO(MO)
    !
    INTEGER                             :: K, N, KGDSI11, KGDSO11

    type(grib1_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out
    integer :: ibo_array(1)

    ! Can't pass expression (e.g. [ibo]) to intent(out) argument.
    ! Initialize placeholder array of size 1 to make rank match.
    ibo_array(1) = ibo

    IF(KGDSI(1).EQ.203) THEN
       KGDSI11=KGDSI(11)
       KGDSI(11)=IOR(KGDSI(11),256)
    ENDIF
    IF(KGDSO(1).EQ.203) THEN
       KGDSO11=KGDSO(11)
       KGDSO(11)=IOR(KGDSO(11),256)
    ENDIF

    desc_in = init_descriptor(kgdsi)
    desc_out = init_descriptor(kgdso)

    call init_grid(grid_in, desc_in)
    call init_grid(grid_out, desc_out)

    CALL ipolatev_grid(ip,IPOPT,grid_in,grid_out, &
         MI,MO,KM,[IBI],LI,UI,VI,&
         NO,RLAT,RLON,CROT,SROT,IBO_array,LO,UO,VO,IRET)

    ibo = ibo_array(1)

    IF(KGDSI(1).EQ.203) THEN
       KGDSI(11)=KGDSI11
    ENDIF
    IF(KGDSO(1).EQ.203) THEN
       KGDSO(11)=KGDSO11
    ENDIF

  END SUBROUTINE ipolatev_grib1_single_field

end module ipolatev_mod

