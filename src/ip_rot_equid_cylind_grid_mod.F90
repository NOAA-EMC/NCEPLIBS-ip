!> @file
!> @brief Rotated equidistant cylindrical GRIB decoder and grid
!> coordinate transformations for Arakawa grids A through D.
!>
!> @author Mark Iredell, George Gayno, Kyle Gerheiser
!> @date July 2021

!> Rotated equidistant cylindrical GRIB decoder and grid coordinate
!> transformations for Arakawa grids A through D. (To handle the E
!> grid, see ip_rot_equid_cylind_egrid_mod).
!>
!> See more info about [Awakawa
!> grids](https://en.wikipedia.org/wiki/Arakawa_grids).
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.1 Rotate
!> Latitude/Longitude](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-1.shtml).
!>
!> @author Gayno @date 2007-NOV-15
module ip_rot_equid_cylind_grid_mod
    use iso_fortran_env, only: real64
    use ip_grid_descriptor_mod
    use ip_grid_mod
    use ip_constants_mod, only: dpr, pi
    use earth_radius_mod
    implicit none

    private
    public :: ip_rot_equid_cylind_grid

    integer, parameter :: kd = real64 !< Fortran kind for reals.

    type, extends(ip_grid) :: ip_rot_equid_cylind_grid
        real(kd) :: clat0 !< Cosine of the latitude of the southern pole of projection.
        real(kd) :: dlats !< 'J'-direction grid increment.
        real(kd) :: dlons !< 'I'-direction grid increment.
        real(kd) :: rlon0 !< Longitude of southern pole of projection.
        real(kd) :: slat0 !< Sine of the latitude of the southern pole of projection.
        real(kd) :: wbd !< Longitude of the western boundary of the grid before rotation.
        real(kd) :: sbd !<  Latitude of the southern boundary of the grid before rotation.
        !> Rotation flag. When '0' the u/v vector components are relative
        !> to north/east. When '1' the u/v vector components are grid
        !> relative.
        integer :: irot
    contains
        !> Initializes a Rotated equidistant cylindrical grid given a
        !> grib1_descriptor object. @return N/A
        procedure :: init_grib1
        !> Initializes a Rotated equidistant cylindrical given a
        !> grib2_descriptor object. @return N/A
        procedure :: init_grib2
        !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
        !> for Gaussian grids. @return N/A
        procedure :: gdswzd => gdswzd_rot_equid_cylind
    end type ip_rot_equid_cylind_grid

    integer :: irot !< Local copy of irot.
    real(KIND=kd) :: rerth !< Radius of the Earth.
    real(KIND=kd) :: clat0 !< Local copy of clat0.
    real(KIND=kd) :: dlats !< Local copy of dlats.
    real(KIND=kd) :: dlons !< Local copy of dlons.
    real(KIND=kd) :: rlon0 !< Local copy of rlon0.
    real(KIND=kd) :: slat0 !< Local copy of slat0.

contains

    !> Initializes a Rotated equidistant cylindrical grid given a
    !> grib1_descriptor object.
    !>
    !> @param[inout] self The grid to initialize
    !> @param[in] g1_desc A grib1_descriptor
    !>
    !> @author Gayno @date 2007-NOV-15
    subroutine init_grib1(self, g1_desc)
        class(ip_rot_equid_cylind_grid), intent(inout) :: self
        type(grib1_descriptor), intent(in) :: g1_desc

        real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
        real(kd) :: hs, hs2, slat1, slat2, slatr, clon1, clon2, clat1, clat2, clatr, clonr, rlonr, rlatr

        associate (kgds => g1_desc%gds)
            self%rerth = 6.3712e6_kd
            self%eccen_squared = 0d0

            rlat1 = kgds(4)*1.e-3_kd
            rlon1 = kgds(5)*1.e-3_kd
            rlat0 = kgds(7)*1.e-3_kd
            self%rlon0 = kgds(8)*1.e-3_kd
            rlat2 = kgds(12)*1.e-3_kd
            rlon2 = kgds(13)*1.e-3_kd

            self%irot = mod(kgds(6)/8, 2)
            self%im = kgds(2)
            self%jm = kgds(3)

            slat1 = sin(rlat1/dpr)
            clat1 = cos(rlat1/dpr)
            self%slat0 = sin(rlat0/dpr)
            self%clat0 = cos(rlat0/dpr)

            hs = sign(1._kd, mod(rlon1-self%rlon0+180+3600, 360._kd)-180)
            clon1 = cos((rlon1-self%rlon0)/dpr)
            slatr = self%clat0*slat1-self%slat0*clat1*clon1
            clatr = sqrt(1-slatr**2)
            clonr = (self%clat0*clat1*clon1+self%slat0*slat1)/clatr
            rlatr = dpr*asin(slatr)
            rlonr = hs*dpr*acos(clonr)

            self%wbd = rlonr
            self%sbd = rlatr
            slat2 = sin(rlat2/dpr)
            clat2 = cos(rlat2/dpr)
            hs2 = sign(1._kd, mod(rlon2-self%rlon0+180+3600, 360._kd)-180)
            clon2 = cos((rlon2-self%rlon0)/dpr)
            slatr = self%clat0*slat2-self%slat0*clat2*clon2
            clatr = sqrt(1-slatr**2)
            clonr = (self%clat0*clat2*clon2+self%slat0*slat2)/clatr
            nbd = dpr*asin(slatr)
            ebd = hs2*dpr*acos(clonr)
            self%dlats = (nbd-self%sbd)/float(self%jm-1)
            self%dlons = (ebd-self%wbd)/float(self%im-1)

            self%iwrap = 0
            self%jwrap1 = 0
            self%jwrap2 = 0
            self%nscan = mod(kgds(11)/32, 2)
            self%nscan_field_pos = self%nscan
            self%kscan = 0
        end associate

    end subroutine init_grib1

    !> Initializes a Rotated equidistant cylindrical grid given a
    !> grib2_descriptor object.
    !>
    !> @param[inout] self The grid to initialize
    !> @param[in] g2_desc A grib2_descriptor
    !>
    !> @author Gayno @date 2007-NOV-15
    subroutine init_grib2(self, g2_desc)
        class(ip_rot_equid_cylind_grid), intent(inout) :: self
        type(grib2_descriptor), intent(in) :: g2_desc

        real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
        integer :: iscale
        integer :: i_offset_odd, i_offset_even, j_offset

        associate (igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

            call earth_radius(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

            i_offset_odd = mod(igdtmpl(19)/8, 2)
            i_offset_even = mod(igdtmpl(19)/4, 2)
            j_offset = mod(igdtmpl(19)/2, 2)

            iscale = igdtmpl(10)*igdtmpl(11)
            if (iscale .eq. 0) iscale = 10**6

            rlat1 = float(igdtmpl(12))/float(iscale)
            rlon1 = float(igdtmpl(13))/float(iscale)
            rlat0 = float(igdtmpl(20))/float(iscale)
            rlat0 = rlat0+90.0_kd

            self%rlon0 = float(igdtmpl(21))/float(iscale)

            rlat2 = float(igdtmpl(15))/float(iscale)
            rlon2 = float(igdtmpl(16))/float(iscale)

            self%irot = mod(igdtmpl(14)/8, 2)
            self%im = igdtmpl(8)
            self%jm = igdtmpl(9)

            self%slat0 = sin(rlat0/dpr)
            self%clat0 = cos(rlat0/dpr)

            self%wbd = rlon1
            if (self%wbd .gt. 180.0) self%wbd = self%wbd-360.0
            self%sbd = rlat1

            nbd = rlat2
            ebd = rlon2

            self%dlats = (nbd-self%sbd)/float(self%jm-1)
            self%dlons = (ebd-self%wbd)/float(self%im-1)

            if (i_offset_odd .eq. 1) self%wbd = self%wbd+(0.5_kd*self%dlons)
            if (j_offset .eq. 1) self%sbd = self%sbd+(0.5_kd*self%dlats)

            self%iwrap = 0
            self%jwrap1 = 0
            self%jwrap2 = 0
            self%kscan = 0
            self%nscan = mod(igdtmpl(19)/32, 2)
            self%nscan_field_pos = self%nscan
        end associate
    end subroutine init_grib2

    !> GDS wizard for rotated equidistant cylindrical.
    !>
    !> This subprogram decodes the grib 2 grid definition template
    !> (passed in integer form as decoded by the ncep g2 library) and
    !> returns one of the following:
    !> - (iopt=+1) earth coordinates of selected grid coordinates
    !> - (iopt=-1) grid coordinates of selected earth coordinates
    !>
    !> Works for non-"e" staggered rotated equidistant cylindrical
    !> projections. the scan mode (section 3, octet 72, bits 5-6)
    !> determine whether this is an "h" or "v" grid.
    !>
    !> If the selected coordinates are more than one gridpoint beyond
    !> the the edges of the grid domain, then the relevant output
    !> elements are set to fill values. The actual number of valid
    !> points computed is returned too.
    !>
    !> Optionally, the vector rotations, the map jacobians and the grid
    !> box areas may be returned as well.
    !>
    !> To compute the vector rotations, the optional arguments 'srot'
    !> and 'crot' must be present. To compute the map jacobians, the
    !> optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be
    !> present. To compute the grid box areas, the optional argument
    !> 'area' must be present.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2010-jan-15 | gayno | based on routines gdswzdcb and gdswzdca
    !> 2015-jan-21 | gayno | merger of gdswizcd and gdswzdcd. make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. make part of a module. move vector rotation, map jacobian and grid box area computations to separate subroutines.
    !> 2015-jul-13 | gayno | convert to grib 2. replace grib 1 kgds array with grib 2 grid definition template array. rename as "gdswzd_rot_equid_cylind."
    !> 2018-07-20 | wesley | add threads.
    !>
    !> @param[in] self Module reference.
    !> @param[in] iopt integer option flag
    !> - 1 to compute earth coords of selected grid coords
    !> - -1 to compute grid coords of selected earth coords
    !> @param[in] npts integer maximum number of coordinates
    !> @param[in] fill real fill value to set invalid output data
    !> (must be impossible value; suggested value: -9999.)
    !> @param[inout] xpts real (npts) grid x point coordinates if iopt>0
    !> @param[inout] ypts real (npts) grid y point coordinates if iopt>0
    !> @param[inout] rlon real (npts) earth longitudes in degrees e if iopt<0
    !> (acceptable range: -360. to 360.)
    !> @param[inout] rlat real (npts) earth latitudes in degrees n if iopt<0
    !> (acceptable range: -90. to 90.)
    !> @param[out] nret integer number of valid points computed
    !> @param[out] crot real, optional (npts) clockwise vector rotation cosines
    !> @param[out] srot real, optional (npts) clockwise vector rotation sines
    !> (ugrid=crot*uearth-srot*vearth;
    !> vgrid=srot*uearth+crot*vearth)
    !> @param[out] xlon real, optional (npts) dx/dlon in 1/degrees
    !> @param[out] xlat real, optional (npts) dx/dlat in 1/degrees
    !> @param[out] ylon real, optional (npts) dy/dlon in 1/degrees
    !> @param[out] ylat real, optional (npts) dy/dlat in 1/degrees
    !> @param[out] area real, optional (npts) area weights in m**2
    !>
    !> @author Gayno @date 2007-NOV-15
    subroutine gdswzd_rot_equid_cylind(self, iopt, npts, &
                                       fill, xpts, ypts, rlon, rlat, nret, &
                                       crot, srot, xlon, xlat, ylon, ylat, area)
        implicit none

        class(ip_rot_equid_cylind_grid), intent(in) :: self
        integer, intent(in) :: iopt, npts
        integer, intent(out) :: nret
        !
        real, intent(in) :: fill
        real, intent(inout) :: rlon(npts), rlat(npts)
        real, intent(inout) :: xpts(npts), ypts(npts)
        real, optional, intent(out) :: crot(npts), srot(npts)
        real, optional, intent(out) :: xlon(npts), xlat(npts)
        real, optional, intent(out) :: ylon(npts), ylat(npts), area(npts)
        !
        integer                                :: im, jm, n
        !
        logical                                :: lrot, lmap, larea
        !
        real(KIND=kd)                          :: hs
        real(KIND=kd)                          :: clonr, clatr, slatr
        real(KIND=kd)                          :: clat, slat, clon
        real(KIND=kd)                          :: rlatr, rlonr
        real(KIND=kd)                          :: wbd, sbd
        real                                   :: xmin, xmax, ymin, ymax
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (present(crot)) crot = fill
        if (present(srot)) srot = fill
        if (present(xlon)) xlon = fill
        if (present(xlat)) xlat = fill
        if (present(ylon)) ylon = fill
        if (present(ylat)) ylat = fill
        if (present(area)) area = fill
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! IS THE EARTH RADIUS DEFINED?
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  IS THIS AN "E"-STAGGER GRID?  ROUTINE CAN'T PROCESS THOSE.
        ! I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
        ! I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
        ! J_OFFSET=MOD(IGDTMPL(19)/2,2)
        ! IF(I_OFFSET_ODD/=I_OFFSET_EVEN) THEN
        !    CALL ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
        !    RETURN
        ! ENDIF
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        rlon0 = self%rlon0
        irot = self%irot

        im = self%im
        jm = self%jm

        slat0 = self%slat0
        clat0 = self%clat0

        wbd = self%wbd
        sbd = self%sbd

        dlats = self%dlats
        dlons = self%dlons

        xmin = 0
        xmax = im+1
        ymin = 0
        ymax = jm+1
        nret = 0

        rerth = self%rerth
        if (rerth .lt. 0.) then
            call rot_equid_cylind_error(iopt, fill, rlat, rlon, xpts, ypts, npts)
            return
        end if

        if (present(crot) .and. present(srot)) then
            lrot = .true.
        else
            lrot = .false.
        end if
        if (present(xlon) .and. present(xlat) .and. present(ylon) .and. present(ylat)) then
            lmap = .true.
        else
            lmap = .false.
        end if
        if (present(area)) then
            larea = .true.
        else
            larea = .false.
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
        if (iopt .eq. 0 .or. iopt .eq. 1) then
            !$omp parallel do private(n, rlonr, rlatr, hs, clonr, slatr, clatr, slat, clat, clon) &
                !$omp&reduction(+:nret) schedule(static)
            do n = 1, npts
                if (xpts(n) .ge. xmin .and. xpts(n) .le. xmax .and. &
                    ypts(n) .ge. ymin .and. ypts(n) .le. ymax) then
                    rlonr = wbd+(xpts(n)-1._kd)*dlons
                    rlatr = sbd+(ypts(n)-1._kd)*dlats
                    if (rlonr .le. 0._kd) then
                        hs = -1.0_kd
                    else
                        hs = 1.0_kd
                    end if
                    clonr = cos(rlonr/dpr)
                    slatr = sin(rlatr/dpr)
                    clatr = cos(rlatr/dpr)
                    slat = clat0*slatr+slat0*clatr*clonr
                    if (slat .le. -1) then
                        clat = 0.
                        clon = cos(rlon0/dpr)
                        rlon(n) = 0.
                        rlat(n) = -90.
                    elseif (slat .ge. 1) then
                        clat = 0.
                        clon = cos(rlon0/dpr)
                        rlon(n) = 0.
                        rlat(n) = 90.
                    else
                        clat = sqrt(1-slat**2)
                        clon = (clat0*clatr*clonr-slat0*slatr)/clat
                        clon = min(max(clon, -1._kd), 1._kd)
                        rlon(n) = real(mod(rlon0+hs*dpr*acos(clon)+3600, 360._kd))
                        rlat(n) = real(dpr*asin(slat))
                    end if
                    nret = nret+1
                    if (lrot) call rot_equid_cylind_vect_rot(rlon(n), clatr, slatr, &
                                                             clat, slat, clon, crot(n), srot(n))
                    if (lmap) call rot_equid_cylind_map_jacob(fill, rlon(n), clatr, &
                                                              clat, slat, clon, xlon(n), xlat(n), ylon(n), ylat(n))
                    if (larea) call rot_equid_cylind_grid_area(clatr, fill, area(n))
                else
                    rlon(n) = fill
                    rlat(n) = fill
                end if
            end do
            !$omp end parallel do
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
        elseif (iopt .eq. -1) then
            !$omp parallel do private(n, hs, clon, slat, clat, slatr, clatr, clonr, rlonr, rlatr) &
                !$omp&reduction(+:nret) schedule(static)
            do n = 1, npts
                if (abs(rlon(n)) .le. 360 .and. abs(rlat(n)) .le. 90) then
                    hs = sign(1._kd, mod(rlon(n)-rlon0+180+3600, 360._kd)-180)
                    clon = cos((rlon(n)-rlon0)/dpr)
                    slat = sin(rlat(n)/dpr)
                    clat = cos(rlat(n)/dpr)
                    slatr = clat0*slat-slat0*clat*clon
                    if (slatr .le. -1) then
                        clatr = 0._kd
                        rlonr = 0.
                        rlatr = -90.
                    elseif (slatr .ge. 1) then
                        clatr = 0._kd
                        rlonr = 0.
                        rlatr = 90.
                    else
                        clatr = sqrt(1-slatr**2)
                        clonr = (clat0*clat*clon+slat0*slat)/clatr
                        clonr = min(max(clonr, -1._kd), 1._kd)
                        rlonr = hs*dpr*acos(clonr)
                        rlatr = dpr*asin(slatr)
                    end if
                    xpts(n) = real((rlonr-wbd)/dlons+1._kd)
                    ypts(n) = real((rlatr-sbd)/dlats+1._kd)
                    if (xpts(n) .ge. xmin .and. xpts(n) .le. xmax .and. &
                        ypts(n) .ge. ymin .and. ypts(n) .le. ymax) then
                        nret = nret+1
                        if (lrot) call rot_equid_cylind_vect_rot(rlon(n), clatr, slatr, &
                                                                 clat, slat, clon, crot(n), srot(n))
                        if (lmap) call rot_equid_cylind_map_jacob(fill, rlon(n), clatr, &
                                                                  clat, slat, clon, xlon(n), xlat(n), ylon(n), ylat(n))
                        if (larea) call rot_equid_cylind_grid_area(clatr, fill, area(n))
                    else
                        xpts(n) = fill
                        ypts(n) = fill
                    end if
                else
                    xpts(n) = fill
                    ypts(n) = fill
                end if
            end do
            !$omp end parallel do
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine gdswzd_rot_equid_cylind

    !> Error handler.
    !>
    !> Upon an error, this subprogram assigns a "fill" value to the
    !> output fields.
    !>
    !> @param[in] iopt integer option flag
    !> - +1 to compute earth coords of selected grid coords
    !> - -1 to compute grid coords of selected earth coords
    !> @param[in] fill real fill value to set invalid output data
    !> (must be impossible value; suggested value: -9999.)
    !> @param[out] rlat real (npts) earth latitudes in degrees n if iopt<0
    !> @param[out] rlon real (npts) earth longitudes in degrees e if iopt<0
    !> @param[out] xpts real (npts) grid x point coordinates if iopt>0
    !> @param[out] ypts real (npts) grid y point coordinates if iopt>0
    !> @param[in] npts integer maximum number of coordinates
    !>
    !> @author Gayno @date 2015-07-13
    subroutine rot_equid_cylind_error(iopt, fill, rlat, rlon, xpts, ypts, npts)
        implicit none
        !
        integer, intent(in) :: iopt, npts
        !
        real, intent(in) :: fill
        real, intent(out) :: rlat(npts), rlon(npts)
        real, intent(out) :: xpts(npts), ypts(npts)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (iopt .ge. 0) then
            rlon = fill
            rlat = fill
        end if
        if (iopt .le. 0) then
            xpts = fill
            ypts = fill
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine rot_equid_cylind_error

    !> Vector rotation fields for rotated equidistant cylindrical grids -
    !> non "e" stagger.
    !>
    !> This subprogram computes the vector rotation sines and cosines
    !> for a rotated equidistant cylindrical grid - non "e" stagger.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-07-19 | gayno | rename as "rot_equid_cylind_vect_rot."
    !> 2018-07-20 | wesley | pass in clatr, slatr, clat, slat, clon for threading.
    !>
    !> @param[in] rlon longitude in degrees (real)
    !> @param[in] clatr cosine of rotated latitude (real)
    !> @param[in] slatr sine of rotated latitude (real)
    !> @param[in] clat cosine of latitude (real)
    !> @param[in] slat sine of latitude (real)
    !> @param[in] clon cosine of longitude (real)
    !> @param[out] crot clockwise vector rotation cosines (real)
    !> @param[out] srot clockwise vector rotation sines (real)
    !> (ugrid=crot*uearth-srot*vearth;
    !> vgrid=srot*uearth+crot*vearth)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine rot_equid_cylind_vect_rot(rlon, clatr, slatr, clat, slat, &
                                         clon, crot, srot)
        implicit none

        real(KIND=kd), intent(in) :: clat, clatr, clon, slat, slatr
        real, intent(in) :: rlon
        real, intent(out) :: crot, srot

        real(KIND=kd)                   :: slon

        if (irot .eq. 1) then
            if (clatr .le. 0._kd) then
                crot = real(-sign(1._kd, slatr*slat0))
                srot = 0.
            else
                slon = sin((rlon-rlon0)/dpr)
                crot = real((clat0*clat+slat0*slat*clon)/clatr)
                srot = real(slat0*slon/clatr)
            end if
        else
            crot = 1.
            srot = 0.
        end if

    end subroutine rot_equid_cylind_vect_rot

    !> Map jacobians for rotated equidistant cylindrical
    !> grids - non "e" stagger.
    !>
    !> This subprogram computes the map jacobians for a rotated
    !> equidistant cylindrical grid - non "e" stagger.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "rot_equid_cylind_map_jacob".
    !> 2018-07-20 | wesley | pass in clatr, clat, slat, clon to allow threading.
    !>
    !> @param[in] fill fill value for undefined points (real)
    !> @param[in] rlon longitude in degrees (real)
    !> @param[in] clatr cosine of unrotated latitude (real)
    !> @param[in] clat cosine of latitude (real)
    !> @param[in] slat sine of latitude (real)
    !> @param[in] clon cosine of latitude (real)
    !> @param[out] xlon dx/dlon in 1/degrees (real)
    !> @param[out] xlat dx/dlat in 1/degrees (real)
    !> @param[out] ylon dy/dlon in 1/degrees (real)
    !> @param[out] ylat dy/dlat in 1/degrees (real)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine rot_equid_cylind_map_jacob(fill, rlon, clatr, clat, &
                                          slat, clon, xlon, xlat, ylon, ylat)
        implicit none

        real(KIND=kd), intent(in) :: clatr, clat, slat, clon
        real, intent(in) :: fill, rlon
        real, intent(out) :: xlon, xlat, ylon, ylat

        real(KIND=kd)                   :: slon, term1, term2

        if (clatr .le. 0._kd) then
            xlon = fill
            xlat = fill
            ylon = fill
            ylat = fill
        else
            slon = sin((rlon-rlon0)/dpr)
            term1 = (clat0*clat+slat0*slat*clon)/clatr
            term2 = slat0*slon/clatr
            xlon = real(term1*clat/(dlons*clatr))
            xlat = real(-term2/(dlons*clatr))
            ylon = real(term2*clat/dlats)
            ylat = real(term1/dlats)
        end if

    end subroutine rot_equid_cylind_map_jacob

    !> Grid box area for rotated equidistant cylindrical grids - non "e"
    !> stagger.
    !>
    !> This subprogram computes the grid box area for a rotated
    !> equidistant cylindrical grid - non "e" stagger.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | Gayno | initial version
    !> 2015-07-19 | gayno | rename as "rot_equid_cylind_grid_area."
    !> 2018-07-20 | wesley | pass in clatr for threading
    !>
    !> @param[in] clatr cosine of unrotated latitude (real)
    !> @param[in] fill fill value for undefined points (real)
    !> @param[out] area area weights in m**2 (real)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine rot_equid_cylind_grid_area(clatr, fill, area)
        implicit none

        real(KIND=kd), intent(in) :: clatr
        real, intent(in) :: fill
        real, intent(out) :: area

        if (clatr .le. 0._kd) then
            area = fill
        else
            area = real(2._kd*(rerth**2)*clatr*(dlons/dpr)*sin(0.5_kd*dlats/dpr))
        end if

    end subroutine rot_equid_cylind_grid_area

end module ip_rot_equid_cylind_grid_mod

