!> @file
!! @brief Top-level driver for scalar interpolation routine ipolates.
!! @author Mark Iredell, Kyle Gerheiser

!> Top-level driver for scalar interpolation interpolation routine ipolates.
!! Ipolates is overloaded with interfaces for grib1 and grib2 descriptors
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
module ipolates_mod
  use ip_interpolators_mod
  use ip_grid_descriptor_mod
  use ip_grid_factory_mod
  use ip_interpolators_mod
  use ip_grid_mod
  implicit none

  private
  public :: ipolates, ipolates_grib1, ipolates_grib2

  interface ipolates
     module procedure ipolates_grib1
     module procedure ipolates_grib2
  end interface ipolates

contains

  !> @brief Interpolates scalar fields between grids given ip_grid objects.
  !! @details Calls the specific interpolation routines on the generic ip_grids created from a grib1/grib2 descriptor.
  !! @param[in]  ip Interpolation method.
  !! @param[in]  ipopt Interpolation options.
  !! @param[in]  grid_in Input grid.
  !! @param[in]  grid_out Output grid object created.
  !! @param[in]  mi   Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !! @param[in]  mo   Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !! @param[in]  km   Number of fields to interpolate.
  !! @param[in]  ibi  Input bitmap flags.
  !! @param[in]  li   Input bitmaps (if respective ibi(k)=1).
  !! @param[in]  gi   Input fields to interpolate.
  !! @param[out] no   Number of output points (only if kgdso(1)<0).
  !! @param[out] rlat Output latitudes in degrees (if kgdso(1)<0).
  !! @param[out] rlon Output longitudes in degrees (if kgdso(1)<0).
  !! @param[out] ibo  Output bitmap flags.
  !! @param[out] lo   Output bitmaps (always output).
  !! @param[out] go   Output fields interpolated.
  !! @param[out] iret Return code.
  !! - 0 Successful interpolation.
  !! - 1 Unrecognized interpolation method.
  !! - 2 Unrecognized input grid or no grid overlap.
  !! - 3 Unrecognized output grid.
  !! - 1x Invalid bicubic method parameters.
  !! - 3x Invalid budget method parameters.
  !! - 4x Invalid spectral method parameters.
  subroutine ipolates_grid(ip, ipopt, grid_in, grid_out, mi, mo, km,&
       & ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    class(ip_grid), intent(in) :: grid_in, grid_out
    INTEGER,    INTENT(IN   ) :: IP, IPOPT(20), KM, MI, MO
    INTEGER,    INTENT(IN   ) :: IBI(KM)
    INTEGER,    INTENT(INOUT) :: NO
    INTEGER,    INTENT(  OUT) :: IRET, IBO(KM)
    !
    LOGICAL*1,  INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,  INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,       INTENT(IN   ) :: GI(MI,KM)
    REAL,       INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,       INTENT(  OUT) :: GO(MO,KM)
    !

    select case(ip)
    case(BILINEAR_INTERP_ID)
       CALL interpolate_bilinear(IPOPT,grid_in,grid_out,MI,MO,KM,IBI&
            &,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
    case(BICUBIC_INTERP_ID)
       CALL interpolate_bicubic(IPOPT,grid_in,grid_out,MI,MO,KM,IBI&
            &,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
    case(NEIGHBOR_INTERP_ID)
       CALL interpolate_neighbor(IPOPT,grid_in,grid_out,MI,MO,KM,IBI&
            &,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
    case(BUDGET_INTERP_ID)
       CALL interpolate_budget(IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI&
            &,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
    case(SPECTRAL_INTERP_ID)
       CALL interpolate_spectral(IPOPT,grid_in,grid_out,MI,MO,KM,IBI&
            &,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
    case(NEIGHBOR_BUDGET_INTERP_ID)
       CALL interpolate_neighbor_budget(IPOPT,grid_in,grid_out,MI,MO&
            &,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
    case default
       ! IF(KGDSO(1).GE.0) NO=0
       ! DO K=1,KM
       !    IBO(K)=1
       !    DO N=1,NO
       !       LO(N,K)=.FALSE.
       !       GO(N,K)=0.
       !    ENDDO
       ! ENDDO
       IRET=1
       print *, "Unrecognized interp option: ", ip
       error stop
    end select

  end subroutine ipolates_grid

  !> @brief This subprogram interpolates scalar field from any grid
  !! to any grid given a grib1 Grid Descriptor Section.
  !!
  !! @details Only horizontal interpolation is performed.
  !! The following interpolation methods are possible:
  !! - (ip=0) bilinear
  !! - (ip=1) bicubic
  !! - (ip=2) neighbor
  !! - (ip=3) budget
  !! - (ip=4) spectral
  !! - (ip=6) neighbor-budget
  !!
  !! Some of these methods have interpolation options and/or
  !! restrictions on the input or output grids, both of which
  !! are documented more fully in their respective subprograms.
  !!
  !! The grids are defined by their grid description sections
  !! (passed in integer form as decoded by subprogram w3fi63).
  !!
  !! The current code recognizes the following projections:
  !! - (kgds(1)=000) equidistant cylindrical
  !! - (kgds(1)=001) mercator cylindrical
  !! - (kgds(1)=003) lambert conformal conical
  !! - (kgds(1)=004) gaussian cylindrical
  !! - (kgds(1)=005) polar stereographic azimuthal
  !! - (kgds(1)=203) rotated equidistant cylindrical - e-stagger
  !! - (kgds(1)=205) rotated equidistant cylindrical - b-stagger
  !!
  !! Where kgds could be either input kgdsi or output kgdso.
  !!
  !! As an added bonus the number of output grid points
  !! and their latitudes and longitudes are also returned.
  !!
  !! On the other hand, the output can be a set of station points
  !! if kgdso(1)<0, in which case the number of points
  !! and their latitudes and longitudes must be input.
  !! for the budget approach, a subsection of the grid may
  !! be output by subtracting kgdso(1) from 255 and passing
  !! in the latitudes and longitudes of the points.
  !! Input bitmaps will be interpolated to output bitmaps.
  !!
  !! Output bitmaps will also be created when the output grid
  !! extends outside of the domain of the input grid.
  !! the output field is set to 0 where the output bitmap is off.
  !!        
  !! @param ip Interpolation method
  !! - ip = BILINEAR_INTERP_ID = 0 for bilinear
  !! - ip = BICUBIC_INTERP_ID = 1 for bicubic
  !! - ip = NEIGHBOR_INTERP_ID = 2 for neighbor;
  !! - ip = BUDGET_INTERP_ID = 3 for budget;
  !! - ip = SPECTRAL_INTERP_ID = 4 for spectral;
  !! - ip = NEIGHBOR_BUDGET_INTERP_ID = 6 for neighbor-budget
  !!
  !! @param ipopt Interpolation options
  !! - ip=0 (bilinear): (No options)
  !! - ip=1 Cbicubic): constraint option
  !! - ip=2 (neighbor): (No options)
  !! - ip=3 (budget): Number in radius, radius weights, search radius
  !! - ip=4 (spectral): Spectral shape, spectral truncation
  !! - ip=6 (neighbor-budget): Number in radius, radius weights ...)
  !!
  !! @param[in] kgdsi Input gds parameters as decoded by w3fi63.
  !! @param[in] kgdso Output gds parameters.
  !! @param[in] mi    Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !! @param[in] mo    Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !! @param[in] km    Number of fields to interpolate.
  !! @param[in] ibi   Input bitmap flags.
  !! @param[in] li    Input bitmaps (if respective ibi(k)=1).
  !! @param[in] gi    Input fields to interpolate.
  !! @param[out] no Number of output points (only if kgdso(1)<0).
  !! @param[out] rlat Output latitudes in degrees (if kgdso(1)<0).
  !! @param[out] rlon Output longitudes in degrees (if kgdso(1)<0).
  !! @param[out] ibo Output bitmap flags.
  !! @param[out] lo  Output bitmaps (always output).
  !! @param[out] go  Output fields interpolated.
  !! @param[out] iret Return code.
  !! - 0 Successful interpolation.
  !! - 1 Unrecognized interpolation method.
  !! - 2 Unrecognized input grid or no grid overlap.
  !! - 3 Unrecognized output grid.
  !! - 1x Invalid bicubic method parameters.
  !! - 3x Invalid budget method parameters.
  !! - 4x Invalid spectral method parameters.
  subroutine ipolates_grib1(ip,ipopt,kgdsi,kgdso,mi,mo,km,ibi,li,gi, &
       no,rlat,rlon,ibo,lo,go,iret) bind(c)
    !
    INTEGER,    INTENT(IN   ) :: IP, IPOPT(20), KM, MI, MO
    INTEGER,    INTENT(IN   ) :: IBI(KM), KGDSI(200), KGDSO(200)
    INTEGER,    INTENT(INOUT) :: NO
    INTEGER,    INTENT(  OUT) :: IRET, IBO(KM)
    !
    LOGICAL*1,  INTENT(IN   ) :: LI(MI,KM)
    LOGICAL*1,  INTENT(  OUT) :: LO(MO,KM)
    !
    REAL,       INTENT(IN   ) :: GI(MI,KM)
    REAL,       INTENT(INOUT) :: RLAT(MO),RLON(MO)
    REAL,       INTENT(  OUT) :: GO(MO,KM)
    !
    INTEGER                   :: K, N

    type(grib1_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(kgdsi)
    desc_out = init_descriptor(kgdso)

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)

    call ipolates_grid(ip, ipopt, grid_in, grid_out, mi, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)

  END SUBROUTINE IPOLATES_GRIB1


  !> @brief This subprogram interpolates scalar field from any grid to any grid given a grib2 descriptor.
  !! @details Wrapper for ipolates_grid which converts a grib1 descriptor into an ip_grid_descriptor,
  !! which is used to create an ip_grid.
  !! Only horizontal interpolation is performed.
  !!
  !! The following interpolation methods are possible:
  !! - (ip=0) bilinear
  !! - (ip=1) bicubic
  !! - (ip=2) neighbor
  !! - (ip=3) budget
  !! - (ip=4) spectral
  !! - (ip=6) neighbor-budget
  !!
  !! Some of these methods have interpolation options and/or
  !! restrictions on the input or output grids, both of which
  !! are documented more fully in their respective subprograms.
  !!
  !! Input and output grids are defined by their grib 2 grid
  !! definition template as decoded by the ncep g2 library. The
  !! current code recognizes the following projections, where
  !! "igdtnumi/o" is the grib 2 grid defintion template number
  !! for the input and output grids, respectively:
  !! - (igdtnumi/o=00) equidistant cylindrical
  !! - (igdtnumi/o=01) rotated equidistant cylindrical. "e" and non-"e" staggered
  !! - (igdtnumi/o=10) mercator cylindrical
  !! - (igdtnumi/o=20) polar stereographic azimuthal
  !! - (igdtnumi/o=30) lambert conformal conical
  !! - (igdtnumi/o=40) gaussian cylindrical
  !!
  !! As an added bonus the number of output grid points
  !! and their latitudes and longitudes are also returned.
  !!
  !! On the other hand, data may be interpolated to a set of station
  !! points if "igdtnumo"<0 (or subtracted from 255 for the budget 
  !! option), in which case the number of points and
  !! their latitudes and longitudes must be input.
  !!
  !! Input bitmaps will be interpolated to output bitmaps.
  !! Output bitmaps will also be created when the output grid
  !! extends outside of the domain of the input grid.
  !!
  !! The output field is set to 0 where the output bitmap is off.
  !!
  !! @param[in] ip Interpolation method
  !! - ip=0 for bilinear
  !! - ip=1 for bicubic
  !! - ip=2 for neighbor;
  !! - ip=3 for budget;
  !! - ip=4 for spectral;
  !! - ip=6 for neighbor-budget
  !!
  !! @param[in] ipopt Interpolation options
  !! - ip=0: (No options)
  !! - ip=1: Constraint option
  !! - ip=2: (No options)
  !! - ip=3: Number in radius, radius weights, search radius
  !! - ip=4: Spectral shape, spectral truncation
  !! - ip=6: Number in radius, radius weights ...)
  !!
  !! @param[in] igdtnumi Grid definition template number for the input grid.
  !! Corresponds to the gfld%igdtnum component of the ncep g2 library gridmod data structure:
  !! - 00 - EQUIDISTANT CYLINDRICAL
  !! - 01 - Rotated equidistant cylindrical. "e" and non-"e" staggered
  !! - 10 - MERCATOR CYCLINDRICAL
  !! - 20 - POLAR STEREOGRAPHIC AZIMUTHAL
  !! - 30 - LAMBERT CONFORMAL CONICAL
  !! - 40 - GAUSSIAN EQUIDISTANT CYCLINDRICAL
  !!
  !! @param[in] igdtmpli Grid definition template array input grid.
  !! Corresponds to the gfld%igdtmpl component of the NCEPLIBS-g2 gridmod data structure
  !!
  !! Section 3 Info:
  !!
  !! All map projections:
  !! - (1): SHAPE OF EARTH, OCTET 15
  !! - (2): SCALE FACTOR OF SPHERICAL EARTH RADIUS, OCTET 16
  !! - (3): SCALED VALUE OF RADIUS OF SPHERICAL EARTH, OCTETS 17-20
  !! - (4): SCALE FACTOR OF MAJOR AXIS OF ELLIPTICAL EARTH, OCTET 21
  !! - (5): SCALED VALUE OF MAJOR AXIS OF ELLIPTICAL EARTH, OCTETS 22-25
  !! - (6): SCALE FACTOR OF MINOR AXIS OF ELLIPTICAL EARTH, OCTET 26
  !! - (7): SCALED VALUE OF MINOR AXIS OF ELLIPTICAL EARTH, OCTETS 27-30
  !!
  !! Equidistant Cyclindrical:
  !! - (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
  !! - (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
  !! - (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN, OCTETS 39-42.
  !! - (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
  !! - (12): LATITUDE OF FIRST GRID POINT, OCTETS 47-50
  !! - (13): LONGITUDE OF FIRST GRID POINT, OCTETS 51-54
  !! - (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
  !! - (15): LATITUDE OF LAST GRID POINT, OCTETS 56-59
  !! - (16): LONGITUDE OF LAST GRID POINT, OCTETS 60-63
  !! - (17): I-DIRECTION INCREMENT, OCTETS 64-67
  !! - (18): J-DIRECTION INCREMENT, OCTETS 68-71
  !! - (19): SCANNING MODE, OCTET 72
  !!
  !! Mercator Cyclindrical:
  !! - (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
  !! - (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
  !! - (10): LATITUDE OF FIRST POINT, OCTETS 39-42
  !! - (11): LONGITUDE OF FIRST POINT, OCTETS 43-46
  !! - (12): RESOLUTION AND COMPONENT FLAGS, OCTET 47
  !! - (13): TANGENT LATITUDE, OCTETS 48-51
  !! - (14): LATITUDE OF LAST POINT, OCTETS 52-55
  !! - (15): LONGITUDE OF LAST POINT, OCTETS 56-59
  !! - (16): SCANNING MODE FLAGS, OCTET 60
  !! - (17): ORIENTATION OF GRID, OCTETS 61-64
  !! - (18): LONGITUDINAL GRID LENGTH, OCTETS 65-68
  !! - (19): LATITUDINAL GRID LENGTH, OCTETS 69-72
  !!
  !! Lambert Conformal Conical:
  !! - (8):  NUMBER OF POINTS ALONG X-AXIS, OCTS 31-34
  !! - (9):  NUMBER OF POINTS ALONG Y-AXIS, OCTS 35-38
  !! - (10): LATITUDE OF FIRST POINT, OCTETS 39-42
  !! - (11): LONGITUDE OF FIRST POINT, OCTETS 43-46
  !! - (12): RESOLUTION OF COMPONENT FLAG, OCTET 47
  !! - (13): LATITUDE WHERE GRID LENGTHS SPECIFIED,OCTETS 48-51
  !! - (14): LONGITUDE OF MERIDIAN THAT IS PARALLEL TO Y-AXIS, OCTETS 52-55
  !! - (15): X-DIRECTION GRID LENGTH, OCTETS 56-59
  !! - (16): Y-DIRECTION GRID LENGTH, OCTETS 60-63
  !! - (17): PROJECTION CENTER FLAG, OCTET 64
  !! - (18): SCANNING MODE, OCTET 65
  !! - (19): FIRST TANGENT LATITUDE FROM POLE, OCTETS 66-69
  !! - (20): SECOND TANGENT LATITUDE FROM POLE, OCTETS 70-73
  !! - (21): LATITUDE OF SOUTH POLE OF PROJECTION, OCTETS 74-77
  !! - (22): LONGITUDE OF SOUTH POLE OF PROJECTION, OCTETS 78-81
  !!
  !! Gaussian Cylindrical:
  !! - (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
  !! - (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
  !! - (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN, OCTETS 39-42
  !! - (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
  !! - (12): LATITUDE OF FIRST GRID POINT, OCTETS 47-50
  !! - (13): LONGITUDE OF FIRST GRID POINT, OCTETS 51-54
  !! - (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
  !! - (15): LATITUDE OF LAST GRID POINT, OCTETS 56-59
  !! - (16): LONGITUDE OF LAST GRID POINT, OCTETS 60-63
  !! - (17): I-DIRECTION INCREMENT, OCTETS 64-67
  !! - (18): NUMBER OF PARALLELS BETWEEN POLE AND EQUATOR, OCTETS 68-71
  !! - (19): SCANNING MODE, OCTET 72
  !!
  !! Polar Stereographic Azimuthal:
  !! - (8):  NUMBER OF POINTS ALONG X-AXIS, OCTETS 31-34
  !! - (9):  NUMBER OF POINTS ALONG Y-AXIS, OCTETS 35-38
  !! - (10): LATITUDE OF FIRST GRID POINT, OCTETS 39-42
  !! - (11): LONGITUDE OF FIRST GRID POINT, OCTETS 43-46
  !! - (12): RESOLUTION AND COMPONENT FLAGS, OCTET 47
  !! - (13): TRUE LATITUDE, OCTETS 48-51
  !! - (14): ORIENTATION LONGITUDE, OCTETS 52-55
  !! - (15): X-DIRECTION GRID LENGTH, OCTETS 56-59
  !! - (16): Y-DIRECTION GRID LENGTH, OCTETS 60-63
  !! - (17): PROJECTION CENTER FLAG, OCTET 64
  !! - (18): SCANNING MODE FLAGS, OCTET 65
  !!
  !! Rotated Equidistant Cyclindrical:
  !! - (8):  NUMBER OF POINTS ALONG A PARALLEL, OCTS 31-34
  !! - (9):  NUMBER OF POINTS ALONG A MERIDIAN, OCTS 35-38
  !! - (10): BASIC ANGLE OF INITIAL PRODUCTION DOMAIN, OCTETS 39-42
  !! - (11): SUBDIVISIONS OF BASIC ANGLE, OCTETS 43-46
  !! - (12): LATITUDE OF FIRST GRID POINT, OCTETS 47-50
  !! - (13): LONGITUDE OF FIRST GRID POINT, OCTETS 51-54
  !! - (14): RESOLUTION AND COMPONENT FLAGS, OCTET 55
  !! - (15): LATITUDE OF LAST GRID POINT, OCTETS 56-59
  !! - (16): LONGITUDE OF LAST GRID POINT, OCTETS 60-63
  !! - (17): I-DIRECTION INCREMENT, OCTETS 64-67
  !! - (18): J-DIRECTION INCREMENT, OCTETS 68-71
  !! - (19): SCANNING MODE, OCTET 72
  !! - (20): LATITUDE OF SOUTHERN POLE OF PROJECTION, OCTETS 73-76
  !! - (21): LONGITUDE OF SOUTHERN POLE OF PROJECTION, OCTETS 77-80
  !! - (22): ANGLE OF ROTATION OF PROJECTION, OCTS 81-84
  !!
  !! @param[in] igdtleni Number of elements of the grid definition
  !! template array for the input grid. Corresponds to the gfld%igdtlen
  !! component of the ncep g2 library gridmod data structure.
  !!
  !! @param[in] igdtnumo Grid definition template number for the output grid.
  !! Corresponds to the gfld%igdtnum component of the
  !! ncep g2 library gridmod data structure.
  !! See "igdtnumi" for specific template definitions.
  !! Note: igdtnumo<0 means interpolate to random station points.
  !!
  !! @param[in] igdtmplo Grid definition template array for the output grid.
  !! Corresponds to the gfld%igdtmpl component of the ncep g2 library gridmod data structure.
  !! See "igdtmpli" for definition of array elements.
  !!
  !! @param[in] igdtleno Number of elements of the grid definition template array for the output grid.  c
  !! Corresponds to the gfld%igdtlen component of the ncep g2 library gridmod data structure.
  !!
  !! @param[in] mi    Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !! @param[in] mo    Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !! @param[in] km    Number of fields to interpolate.
  !! @param[in] ibi   Input bitmap flags.
  !! @param[in] li    Input bitmaps (if respective ibi(k)=1).
  !! @para[in] gi    Input fields to interpolate.
  !! @param[out] no Number of output points (only if kgdso(1)<0).
  !! @param[out] rlat Output latitudes in degrees (if kgdso(1)<0).
  !! @param[out] rlon Output longitudes in degrees (if kgdso(1)<0).
  !! @param[out] ibo Output bitmap flags.
  !! @param[out] lo  Output bitmaps (always output).
  !! @param[out] go  Output fields interpolated.
  !! @param[out] iret Return code.
  !! - 0 Successful interpolation.
  !! - 1 Unrecognized interpolation method.
  !! - 2 Unrecognized input grid or no grid overlap.
  !! - 3 Unrecognized output grid.
  !! - 1x Invalid bicubic method parameters.
  !! - 3x Invalid budget method parameters.
  !! - 4x Invalid spectral method parameters.
  !!
  !! @note Examples demonstrating relative cpu costs.
  !! This example is interpolating 12 levels of temperatures
  !! from the 360 x 181 global grid (ncep grid 3)
  !! to the 93 x 68 hawaiian mercator grid (ncep grid 204).
  !!
  !! The example times are for the c90. As a reference, the cp time
  !! for unpacking the global 12 temperature fields is 0.04 seconds.
  !!
  !!   METHOD  | IP| IPOPT       |  CP SECONDS
  !!   --------| --|-------------| ----------
  !!   BILINEAR| 0 |             |   0.03
  !!   BICUBIC | 1 | 0           |   0.07
  !!   BICUBIC | 1 | 1           |   0.07
  !!   NEIGHBOR| 2 |             |   0.01
  !!   BUDGET  | 3 | -1,-1       |   0.48
  !!   SPECTRAL| 4 | 0,40        |   0.22
  !!   SPECTRAL| 4 | 1,40        |   0.24
  !!   SPECTRAL| 4 | 0,-1        |   0.42
  !!   N-BUDGET| 6 | -1,-1       |   0.15
  !!
  !!   The spectral interpolation is fast for the mercator grid.
  !!   However, for some grids the spectral interpolation is slow.
  !!
  !!   The following example is interpolating 12 levels of temperatures
  !!   from the 360 x 181 global grid (ncep grid 3)
  !!   to the 93 x 65 conus lambert conformal grid (ncep grid 211).
  !!
  !!   METHOD  | IP| IPOPT        |CP SECONDS
  !!   --------| --| -------------|----------
  !!   BILINEAR| 0 |              | 0.03
  !!   BICUBIC | 1 | 0            | 0.07
  !!   BICUBIC | 1 | 1            | 0.07
  !!   NEIGHBOR| 2 |              | 0.01
  !!   BUDGET  | 3 | -1,-1        | 0.51
  !!   SPECTRAL| 4 | 0,40         | 3.94
  !!   SPECTRAL| 4 | 1,40         | 5.02
  !!   SPECTRAL| 4 | 0,-1         | 11.36
  !!   N-BUDGET| 6 | -1,-1        | 0.18
  !!
  SUBROUTINE IPOLATES_grib2(IP,IPOPT,IGDTNUMI,IGDTMPLI,IGDTLENI, &
       IGDTNUMO,IGDTMPLO,IGDTLENO, &
       MI,MO,KM,IBI,LI,GI, &
       NO,RLAT,RLON,IBO,LO,GO,IRET) bind(C)
    INTEGER,        INTENT(IN   )     :: IP, IPOPT(20), KM, MI, MO
    INTEGER,        INTENT(IN   )     :: IBI(KM)
    INTEGER,        INTENT(IN   )     :: IGDTNUMI, IGDTLENI
    INTEGER,        INTENT(IN   )     :: IGDTMPLI(IGDTLENI)
    INTEGER,        INTENT(IN   )     :: IGDTNUMO, IGDTLENO
    INTEGER,        INTENT(IN   )     :: IGDTMPLO(IGDTLENO)
    INTEGER,        INTENT(  OUT)     :: NO
    INTEGER,        INTENT(  OUT)     :: IRET, IBO(KM)
    !
    LOGICAL*1,      INTENT(IN   )     :: LI(MI,KM)
    LOGICAL*1,      INTENT(  OUT)     :: LO(MO,KM)
    !
    REAL,           INTENT(IN   )     :: GI(MI,KM)
    REAL,           INTENT(INOUT)     :: RLAT(MO),RLON(MO)
    REAL,           INTENT(  OUT)     :: GO(MO,KM)

    type(grib2_descriptor) :: desc_in, desc_out
    class(ip_grid), allocatable :: grid_in, grid_out

    desc_in = init_descriptor(igdtnumi, igdtleni, igdtmpli)
    desc_out = init_descriptor(igdtnumo, igdtleno, igdtmplo)

    grid_in = init_grid(desc_in)
    grid_out = init_grid(desc_out)

    CALL ipolates_grid(ip,IPOPT,grid_in,grid_out,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)

  END SUBROUTINE IPOLATES_GRIB2



end module ipolates_mod

