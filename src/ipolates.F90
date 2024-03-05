!> @file
!! @brief Top-level driver for scalar interpolation routine ipolates().
!! @author Mark Iredell, Kyle Gerheiser

!> @brief Top-level driver for scalar interpolation interpolation
!! routine ipolates().
!!
!! ipolates() is overloaded with interfaces for GRIB1 and GRIB2
!! descriptors.
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
module ipolates_mod
    use ip_interpolators_mod
    use ip_grid_descriptor_mod
    use ip_grid_factory_mod,only:init_grid
    use ip_interpolators_mod
    use ip_grid_mod
    implicit none

    private
    public :: ipolates,ipolates_grib2,ipolates_grib1_single_field,ipolates_grib1,ipolates_grib2_single_field

    interface ipolates
        module procedure ipolates_grib1
        module procedure ipolates_grib1_single_field
        module procedure ipolates_grib2
        module procedure ipolates_grib2_single_field
    endinterface ipolates

contains

    !> Interpolates scalar fields between grids given ip_grid objects.
  !!
  !! Calls the specific interpolation routines on the generic ip_grids
  !! created from a grib1/grib2 descriptor.
  !!
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
  !!
  !! @author Mark Iredell, Kyle Gerheiser
    subroutine ipolates_grid(ip,ipopt,grid_in,grid_out,mi,mo,km,&
         & ibi,li,gi,no,rlat,rlon,ibo,lo,go,iret)
        class(ip_grid),intent(in) :: grid_in,grid_out
        integer,intent(in) :: ip,ipopt(20),km,mi,mo
        integer,intent(in) :: ibi(km)
        integer,intent(inout) :: no
        integer,intent(out) :: iret,ibo(km)
        !
        logical*1,intent(in) :: li(mi,km)
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: gi(mi,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: go(mo,km)
        !

        select case(ip)
        case(bilinear_interp_id)
            call interpolate_bilinear(ipopt,grid_in,grid_out,mi,mo,km,ibi&
                 &,li,gi,no,rlat,rlon,ibo,lo,go,iret)
        case(bicubic_interp_id)
            call interpolate_bicubic(ipopt,grid_in,grid_out,mi,mo,km,ibi&
                 &,li,gi,no,rlat,rlon,ibo,lo,go,iret)
        case(neighbor_interp_id)
            call interpolate_neighbor(ipopt,grid_in,grid_out,mi,mo,km,ibi&
                 &,li,gi,no,rlat,rlon,ibo,lo,go,iret)
        case(budget_interp_id)
            call interpolate_budget(ipopt,grid_in,grid_out,mi,mo,km,ibi,li&
                 &,gi,no,rlat,rlon,ibo,lo,go,iret)
        case(spectral_interp_id)
            call interpolate_spectral(ipopt,grid_in,grid_out,mi,mo,km,ibi&
                 &,gi,no,rlat,rlon,ibo,lo,go,iret)
        case(neighbor_budget_interp_id)
            call interpolate_neighbor_budget(ipopt,grid_in,grid_out,mi,mo&
                 &,km,ibi,li,gi,no,rlat,rlon,ibo,lo,go,iret)
        case default
            ! IF(KGDSO(1).GE.0) NO=0
            ! DO K=1,KM
            !    IBO(K)=1
            !    DO N=1,NO
            !       LO(N,K)=.FALSE.
            !       GO(N,K)=0.
            !    ENDDO
            ! ENDDO
            iret=1
            print*,"Unrecognized interp option: ",ip
            error stop
        endselect

    endsubroutine ipolates_grid

    !> Special case of ipolates_grib1 when interpolating a single field.
  !! Removes the km dimension of input arrays so scalars can be passed to ibi/ibo.
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
  !!
  !! @date Jan 2022
  !! @author Kyle Gerheiser
    subroutine ipolates_grib1_single_field(ip,ipopt,kgdsi,kgdso,mi,mo,km,ibi,li,gi, &
                                           no,rlat,rlon,ibo,lo,go,iret) bind(c)
        !
        use iso_c_binding,only:c_int,c_float,c_double,c_bool,c_long
#if (LSIZE==8)
        integer(c_long),intent(in) :: ip,ipopt(20),km,mi,mo
        integer(c_long),intent(in) :: ibi,kgdsi(200),kgdso(200)
        integer(c_long),intent(inout) :: no
        integer(c_long),intent(out) :: iret,ibo
#else
        integer(c_int),intent(in) :: ip,ipopt(20),km,mi,mo
        integer(c_int),intent(in) :: ibi,kgdsi(200),kgdso(200)
        integer(c_int),intent(inout) :: no
        integer(c_int),intent(out) :: iret,ibo
#endif
        !
        logical(c_bool),intent(in) :: li(mi)
        logical(c_bool),intent(out) :: lo(mo)
        !
#if (LSIZE==4)
        real(c_float),intent(in) :: gi(mi)
        real(c_float),intent(inout) :: rlat(mo),rlon(mo)
        real(c_float),intent(out) :: go(mo)
#else
        real(c_double),intent(in) :: gi(mi)
        real(c_double),intent(inout) :: rlat(mo),rlon(mo)
        real(c_double),intent(out) :: go(mo)
#endif
        !

        type(grib1_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out
        integer :: ibo_array(1)

        desc_in=init_descriptor(kgdsi)
        desc_out=init_descriptor(kgdso)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)

        ! Can't pass expression (e.g. [ibo]) to intent(out) argument.
        ! Initialize placeholder array of size 1 to make rank match.
        ibo_array(1)=ibo

        call ipolates_grid(ip,ipopt,grid_in,grid_out,mi,mo,km,[ibi],li,gi,no,rlat,rlon,ibo_array,lo,go,iret)

        ibo=ibo_array(1)

    endsubroutine ipolates_grib1_single_field

    !> This subprogram interpolates scalar field from any grid
  !! to any grid given a grib1 Grid Descriptor Section.
  !!
  !! Only horizontal interpolation is performed.
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
  !!
  !! @author Mark Iredell, Kyle Gerheiser
    subroutine ipolates_grib1(ip,ipopt,kgdsi,kgdso,mi,mo,km,ibi,li,gi, &
                              no,rlat,rlon,ibo,lo,go,iret) bind(c)
        !
        use iso_c_binding,only:c_int,c_float,c_double,c_bool,c_long
#if (LSIZE==8)
        integer(c_long),intent(in) :: ip,ipopt(20),km,mi,mo
        integer(c_long),intent(in) :: ibi(km),kgdsi(200),kgdso(200)
        integer(c_long),intent(inout) :: no
        integer(c_long),intent(out) :: iret,ibo(km)
#else
        integer(c_int),intent(in) :: ip,ipopt(20),km,mi,mo
        integer(c_int),intent(in) :: ibi(km),kgdsi(200),kgdso(200)
        integer(c_int),intent(inout) :: no
        integer(c_int),intent(out) :: iret,ibo(km)
#endif
        !
        logical(c_bool),intent(in) :: li(mi,km)
        logical(c_bool),intent(out) :: lo(mo,km)
        !
#if (LSIZE==4)
        real(c_float),intent(in) :: gi(mi,km)
        real(c_float),intent(inout) :: rlat(mo),rlon(mo)
        real(c_float),intent(out) :: go(mo,km)
#else
        real(c_double),intent(in) :: gi(mi,km)
        real(c_double),intent(inout) :: rlat(mo),rlon(mo)
        real(c_double),intent(out) :: go(mo,km)
#endif
        !

        type(grib1_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out

        desc_in=init_descriptor(kgdsi)
        desc_out=init_descriptor(kgdso)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)

        call ipolates_grid(ip,ipopt,grid_in,grid_out,mi,mo,km,ibi,li,gi,no,rlat,rlon,ibo,lo,go,iret)

    endsubroutine ipolates_grib1

    !> This subprogram interpolates scalar field from any grid to any
  !! grid given a grib2 descriptor.
  !!
  !! Wrapper for ipolates_grid which converts a grib1 descriptor into
  !! an ip_grid_descriptor, which is used to create an ip_grid. Only
  !! horizontal interpolation is performed.
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
  !! @author Mark Iredell, Kyle Gerheiser
    subroutine ipolates_grib2(ip,ipopt,igdtnumi,igdtmpli,igdtleni, &
                              igdtnumo,igdtmplo,igdtleno, &
                              mi,mo,km,ibi,li,gi, &
                              no,rlat,rlon,ibo,lo,go,iret) bind(c)
        use iso_c_binding,only:c_int,c_float,c_double,c_bool,c_long
#if (LSIZE==8)
        integer(c_long),intent(in)     :: ip,ipopt(20),km,mi,mo
        integer(c_long),intent(in)     :: ibi(km)
        integer(c_long),intent(in)     :: igdtnumi,igdtleni
        integer(c_long),intent(in)     :: igdtmpli(igdtleni)
        integer(c_long),intent(in)     :: igdtnumo,igdtleno
        integer(c_long),intent(in)     :: igdtmplo(igdtleno)
        integer(c_long),intent(out)     :: no
        integer(c_long),intent(out)     :: iret,ibo(km)
#else
        integer(c_int),intent(in)     :: ip,ipopt(20),km,mi,mo
        integer(c_int),intent(in)     :: ibi(km)
        integer(c_int),intent(in)     :: igdtnumi,igdtleni
        integer(c_int),intent(in)     :: igdtmpli(igdtleni)
        integer(c_int),intent(in)     :: igdtnumo,igdtleno
        integer(c_int),intent(in)     :: igdtmplo(igdtleno)
        integer(c_int),intent(out)     :: no
        integer(c_int),intent(out)     :: iret,ibo(km)
#endif
        !
        logical(c_bool),intent(in)     :: li(mi,km)
        logical(c_bool),intent(out)     :: lo(mo,km)
        !
#if (LSIZE==4)
        real(c_float),intent(in) :: gi(mi,km)
        real(c_float),intent(inout) :: rlat(mo),rlon(mo)
        real(c_float),intent(out) :: go(mo,km)
#else
        real(c_double),intent(in) :: gi(mi,km)
        real(c_double),intent(inout) :: rlat(mo),rlon(mo)
        real(c_double),intent(out) :: go(mo,km)
#endif

        type(grib2_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out

        desc_in=init_descriptor(igdtnumi,igdtleni,igdtmpli)
        desc_out=init_descriptor(igdtnumo,igdtleno,igdtmplo)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)

        call ipolates_grid(ip,ipopt,grid_in,grid_out,mi,mo,km,ibi,li,gi,no,rlat,rlon,ibo,lo,go,iret)

    endsubroutine ipolates_grib2

    !> Special case of ipolates_grib2 when interpolating a single field.
  !! Removes the km dimension of input arrays so scalars can be passed to ibi/ibo.
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
  !!
  !! @author Eric Engle @date November 2022
    subroutine ipolates_grib2_single_field(ip,ipopt,igdtnumi,igdtmpli,igdtleni, &
                                           igdtnumo,igdtmplo,igdtleno, &
                                           mi,mo,km,ibi,li,gi, &
                                           no,rlat,rlon,ibo,lo,go,iret) bind(c)
        use iso_c_binding,only:c_int,c_float,c_double,c_bool,c_long
#if (LSIZE==8)
        integer(c_long),intent(in)     :: ip,ipopt(20),km,mi,mo
        integer(c_long),intent(in)     :: ibi
        integer(c_long),intent(in)     :: igdtnumi,igdtleni
        integer(c_long),intent(in)     :: igdtmpli(igdtleni)
        integer(c_long),intent(in)     :: igdtnumo,igdtleno
        integer(c_long),intent(in)     :: igdtmplo(igdtleno)
        integer(c_long),intent(out)     :: no
        integer(c_long),intent(out)     :: iret,ibo
#else
        integer(c_int),intent(in)     :: ip,ipopt(20),km,mi,mo
        integer(c_int),intent(in)     :: ibi
        integer(c_int),intent(in)     :: igdtnumi,igdtleni
        integer(c_int),intent(in)     :: igdtmpli(igdtleni)
        integer(c_int),intent(in)     :: igdtnumo,igdtleno
        integer(c_int),intent(in)     :: igdtmplo(igdtleno)
        integer(c_int),intent(out)     :: no
        integer(c_int),intent(out)     :: iret,ibo
#endif
        !
        logical(c_bool),intent(in)     :: li(mi)
        logical(c_bool),intent(out)     :: lo(mo)
        !
#if (LSIZE==4)
        real(c_float),intent(in) :: gi(mi)
        real(c_float),intent(inout) :: rlat(mo),rlon(mo)
        real(c_float),intent(out) :: go(mo)
#else
        real(c_double),intent(in) :: gi(mi)
        real(c_double),intent(inout) :: rlat(mo),rlon(mo)
        real(c_double),intent(out) :: go(mo)
#endif

        type(grib2_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out
        integer :: ibo_array(1)

        desc_in=init_descriptor(igdtnumi,igdtleni,igdtmpli)
        desc_out=init_descriptor(igdtnumo,igdtleno,igdtmplo)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)

        ! Can't pass expression (e.g. [ibo]) to intent(out) argument.
        ! Initialize placeholder array of size 1 to make rank match.
        ibo_array(1)=ibo

        call ipolates_grid(ip,ipopt,grid_in,grid_out,mi,mo,km,[ibi],li,gi,no,rlat,rlon,ibo_array,lo,go,iret)

        ibo=ibo_array(1)

    endsubroutine ipolates_grib2_single_field

endmodule ipolates_mod

