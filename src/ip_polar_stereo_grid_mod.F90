!> @file
!> @brief GDS wizard for polar stereographic azimuthal.
!>
!> @author Iredell @date 96-04-10

!> @brief GDS wizard for polar stereographic azimuthal.
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.20 Polar
!> stereographic
!> projection](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-20.shtml).
!>
!> @author Iredell @date 96-04-10
module ip_polar_stereo_grid_mod
    use ip_grid_descriptor_mod
    use ip_grid_mod
    use ip_constants_mod, only: dpr, pi, pi2, pi4, rerth_wgs84, e2_wgs84
    use earth_radius_mod
    implicit none

    private
    public :: ip_polar_stereo_grid

    type, extends(ip_grid) :: ip_polar_stereo_grid
        logical :: elliptical !< When true/false, computations are based on an elliptical/spherical earth.
        real :: rlat1 !< Latitude of the first grid point.
        real :: rlon1 !< Longitude of the first grid point.
        real :: orient !< Orientation longitude.
        real :: h !< Hemisphere flag. 0 - NH; 1 - SH.
        real :: dxs !< 'x'-direction grid length, adjusted by the scanning mode.
        real :: dys !< 'y'-direction grid length, adjusted by the scanning mode.
        real :: slatr !< Standard latitude of grid in radians.
        !> Rotation flag. When '0' the u/v vector components are relative
        !> to north/east. When '1' the u/v vector components are grid
        !> relative.
        integer :: irot
    contains
        procedure :: init_grib1 !< Initializes a grid given a grib1_descriptor object. @return N/A
        procedure :: init_grib2 !< Initializes a grid given a grib2_descriptor object. @return N/A
        !> Calculates Earth coordinates (iopt = 1) or grid coorindates
        !> (iopt = -1). @return N/A
        procedure :: gdswzd => gdswzd_polar_stereo
    end type ip_polar_stereo_grid

    integer :: irot !< Local copy of irot.
    real :: de2 !< Square of DE.
    real :: dxs !< Local copy of dxs.
    real :: dys !< Local copy of dys.
    real :: e2 !< Eccentricity squared.
    real :: rerth !< Radius of the Earth.
    real :: h !< Local copy of h.
    real :: orient !< Local copy of orient.
    real :: tinyreal = tiny(1.0) !< Smallest positive real value (use for equality comparisons)

contains

    !> Initializes a polar stereographic grid given a grib1_descriptor
  !! object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Iredell @date 96-04-10
    subroutine init_grib1(self, g1_desc)
        class(ip_polar_stereo_grid), intent(inout) :: self
        type(grib1_descriptor), intent(in) :: g1_desc

        real, parameter :: slat = 60.0  ! standard latitude according grib1 standard

        real :: dx, dy, hi, hj
        integer :: iproj, iscan, jscan

        associate (kgds => g1_desc%gds)
            self%elliptical = mod(kgds(6)/64, 2) .eq. 1

            if (.not. self%elliptical) then
                self%rerth = 6.3712e6
                self%eccen_squared = 0d0
            else
                self%rerth = rerth_wgs84
                self%eccen_squared = e2_wgs84 !wgs84 datum
            end if

            self%im = kgds(2)
            self%jm = kgds(3)

            self%rlat1 = kgds(4)*1.e-3
            self%rlon1 = kgds(5)*1.e-3

            self%irot = mod(kgds(6)/8, 2)

            self%slatr = slat/dpr

            self%orient = kgds(7)*1.e-3

            dx = kgds(8)
            dy = kgds(9)

            iproj = mod(kgds(10)/128, 2)
            iscan = mod(kgds(11)/128, 2)
            jscan = mod(kgds(11)/64, 2)

            self%h = (-1.)**iproj
            hi = (-1.)**iscan
            hj = (-1.)**(1-jscan)

            if (abs(self%h+1.) .lt. tinyreal) self%orient = self%orient+180.

            self%dxs = dx*hi
            self%dys = dy*hj

            self%iwrap = 0
            self%jwrap1 = 0
            self%jwrap2 = 0
            self%nscan = mod(kgds(11)/32, 2)
            self%nscan_field_pos = self%nscan
            self%kscan = 0
        end associate

    end subroutine init_grib1

    !> Initializes a polar stereographic grid given a grib2_descriptor
  !! object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g2_desc A grib2_descriptor
  !!
  !! @author Iredell @date 96-04-10
    subroutine init_grib2(self, g2_desc)
        class(ip_polar_stereo_grid), intent(inout) :: self
        type(grib2_descriptor), intent(in) :: g2_desc

        real :: slat, dx, dy, hi, hj
        integer :: iproj, iscan, jscan

        associate (igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
            call earth_radius(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

            self%elliptical = self%eccen_squared .gt. 0.0

            self%im = igdtmpl(8)
            self%jm = igdtmpl(9)

            self%rlat1 = float(igdtmpl(10))*1.e-6
            self%rlon1 = float(igdtmpl(11))*1.e-6

            self%irot = mod(igdtmpl(12)/8, 2)

            slat = float(abs(igdtmpl(13)))*1.e-6
            self%slatr = slat/dpr

            self%orient = float(igdtmpl(14))*1.e-6

            dx = float(igdtmpl(15))*1.e-3
            dy = float(igdtmpl(16))*1.e-3

            iproj = mod(igdtmpl(17)/128, 2)
            iscan = mod(igdtmpl(18)/128, 2)
            jscan = mod(igdtmpl(18)/64, 2)

            self%h = (-1.)**iproj
            hi = (-1.)**iscan
            hj = (-1.)**(1-jscan)

            self%dxs = dx*hi
            self%dys = dy*hj

            self%nscan = mod(igdtmpl(18)/32, 2)
            self%nscan_field_pos = self%nscan
            self%iwrap = 0
            self%jwrap1 = 0
            self%jwrap2 = 0
            self%kscan = 0
        end associate
    end subroutine init_grib2

    !> GDS wizard for polar stereographic azimuthal
    !>
    !> This subprogram decodes the grib 2 grid definition template
    !> (passed in integer form as decoded by the ncep g2 library) and
    !> returns one of the following:
    !> - (iopt=+1) earth coordinates of selected grid coordinates
    !> - (iopt=-1) grid coordinates of selected earth coordinates
    !>
    !> Works for polar stereographic azimuthal projections.
    !>
    !> If the selected coordinates are more than one gridpoint beyond
    !> the the edges of the grid domain, then the relevant output
    !> elements are set to fill values.
    !>
    !> The actual number of valid points computed is returned too.
    !>
    !> Optionally, the vector rotations, map jacobians, and grid box
    !> areas may be returned as well. Routine works for both spherical
    !> and elliptical earths with the exception of the map jacobians and
    !> grid box areas, which are only computed for spherical earths.
    !>
    !> To compute the vector rotations, the optional arguments 'srot'
    !> and 'crot' must be present. To compute the map jacobians, the
    !> optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be
    !> present.  to compute the grid box areas, the optional argument
    !> 'area' must be present.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | iredell | Initial
    !> 97-10-20 | iredell | include map options
    !> 09-05-13 | gayno | ensure area always positive
    !> 2015-01-21 | gayno | merger of gdswiz05 and gdswzd05. make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. make part of a module. move vector rotation, map jacobian and grid box area computations to separate subroutines. include option for elliptical earths.
    !> 2015-07-13 | gayno | convert to grib 2. replace grib 1 kgds array with grib 2 grid definition template array. rename routine.
    !> 2018-07-20 | wesley | add threading.
    !>
    !> @param[in] self grid
    !> @param[in] iopt option flag
    !> - 1 to compute earth coords of selected grid coords
    !> - -1 to compute grid coords of selected earth coords
    !> @param[in] npts maximum number of coordinates
    !> @param[in] fill fill value to set invalid output data
    !> (must be impossible value; suggested value: -9999.)
    !> @param[inout] xpts (npts) grid x point coordinates if iopt>0
    !> @param[inout] ypts (npts) grid y point coordinates if iopt>0
    !> @param[inout] rlon (npts) earth longitudes in degrees e if iopt<0
    !> (acceptable range: -360. to 360.)
    !> @param[inout] rlat (npts) earth latitudes in degrees n if iopt<0
    !> (acceptable range: -90. to 90.)
    !> @param[out] nret number of valid points computed
    !> @param[out] crot optional (npts) clockwise vector rotation cosines
    !> @param[out] srot optional (npts) clockwise vector rotation sines
    !> (ugrid=crot*uearth-srot*vearth;
    !> vgrid=srot*uearth+crot*vearth)
    !> @param[out] xlon optional (npts) dx/dlon in 1/degrees
    !> @param[out] xlat optional (npts) dx/dlat in 1/degrees
    !> @param[out] ylon optional (npts) dy/dlon in 1/degrees
    !> @param[out] ylat optional (npts) dy/dlat in 1/degrees
    !> @param[out] area optional (npts) area weights in m**2
    !> (proportional to the square of the map factor)
    !>
    !> @author Iredell @date 96-04-10
    subroutine gdswzd_polar_stereo(self, iopt, npts, &
                                   fill, xpts, ypts, rlon, rlat, nret, &
                                   crot, srot, xlon, xlat, ylon, ylat, area)
        implicit none
        !

        class(ip_polar_stereo_grid), intent(in) :: self
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
        integer                         :: im, jm
        integer                         :: iter, n
        !
        logical                         :: elliptical, lrot, lmap, larea
        !
        real                            :: alat, alat1, along, diff
        real                            :: di, dj, de
        real                            :: dr, e, e_over_2
        real                            :: mc, slatr
        real                            :: rlat1, rlon1, rho, t, tc
        real                            :: xmax, xmin, ymax, ymin
        real                            :: xp, yp, dr2
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (present(crot)) crot = fill
        if (present(srot)) srot = fill
        if (present(xlon)) xlon = fill
        if (present(xlat)) xlat = fill
        if (present(ylon)) ylon = fill
        if (present(ylat)) ylat = fill
        if (present(area)) area = fill

        elliptical = self%elliptical
        im = self%im
        jm = self%jm

        rlat1 = self%rlat1
        rlon1 = self%rlon1

        irot = self%irot
        slatr = self%slatr
        orient = self%orient

        h = self%h
        dxs = self%dxs
        dys = self%dys

        rerth = self%rerth
        e2 = self%eccen_squared
        !
        ! FIND X/Y OF POLE
        if (.not. elliptical) then
            de = (1.+sin(slatr))*rerth
            dr = de*cos(rlat1/dpr)/(1+h*sin(rlat1/dpr))
            xp = 1-h*sin((rlon1-orient)/dpr)*dr/dxs
            yp = 1+cos((rlon1-orient)/dpr)*dr/dys
            de2 = de**2
        else
            e = sqrt(e2)
            e_over_2 = e*0.5
            alat = h*rlat1/dpr
            along = (rlon1-orient)/dpr
            t = tan(pi4-alat/2.)/((1.-e*sin(alat))/ &
                                  (1.+e*sin(alat)))**(e_over_2)
            tc = tan(pi4-slatr/2.)/((1.-e*sin(slatr))/ &
                                    (1.+e*sin(slatr)))**(e_over_2)
            mc = cos(slatr)/sqrt(1.0-e2*(sin(slatr)**2))
            rho = rerth*mc*t/tc
            yp = 1.0+rho*cos(h*along)/dys
            xp = 1.0-rho*sin(h*along)/dxs
        end if ! ELLIPTICAL
        xmin = 0
        xmax = im+1
        ymin = 0
        ymax = jm+1
        nret = 0
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
            if (.not. elliptical) then
                !$omp parallel do private(n, di, dj, dr2) reduction(+:nret) schedule(static)
                do n = 1, npts
                    if (xpts(n) .ge. xmin .and. xpts(n) .le. xmax .and. &
                        ypts(n) .ge. ymin .and. ypts(n) .le. ymax) then
                        di = (xpts(n)-xp)*dxs
                        dj = (ypts(n)-yp)*dys
                        dr2 = di**2+dj**2
                        if (dr2 .lt. de2*1.e-6) then
                            rlon(n) = 0.
                            rlat(n) = h*90.
                        else
                            rlon(n) = mod(orient+h*dpr*atan2(di, -dj)+3600, 360.)
                            rlat(n) = h*dpr*asin((de2-dr2)/(de2+dr2))
                        end if
                        nret = nret+1
                        if (lrot) call polar_stereo_vect_rot(rlon(n), crot(n), srot(n))
                        if (lmap) call polar_stereo_map_jacob(rlon(n), rlat(n), dr2, &
                                                              xlon(n), xlat(n), ylon(n), ylat(n))
                        if (larea) call polar_stereo_grid_area(rlat(n), dr2, area(n))
                    else
                        rlon(n) = fill
                        rlat(n) = fill
                    end if
                end do
                !$omp end parallel do
            else ! ELLIPTICAL
                !$omp parallel do private(n, di, dj, rho, t, along, alat1, alat, diff) &
                    !$omp&reduction(+:nret) schedule(static)
                do n = 1, npts
                    if (xpts(n) .ge. xmin .and. xpts(n) .le. xmax .and. &
                        ypts(n) .ge. ymin .and. ypts(n) .le. ymax) then
                        di = (xpts(n)-xp)*dxs
                        dj = (ypts(n)-yp)*dys
                        rho = sqrt(di*di+dj*dj)
                        t = (rho*tc)/(rerth*mc)
                        if (abs(ypts(n)-yp) .lt. 0.01) then
                            if (di .gt. 0.0) along = orient+h*90.0
                            if (di .le. 0.0) along = orient-h*90.0
                        else
                            along = orient+h*atan(di/(-dj))*dpr
                            if (dj .gt. 0) along = along+180.
                        end if
                        alat1 = pi2-2.0*atan(t)
                        do iter = 1, 10
                            alat = pi2-2.0*atan(t*(((1.0-e*sin(alat1))/ &
                                                    (1.0+e*sin(alat1)))**(e_over_2)))
                            diff = abs(alat-alat1)*dpr
                            if (diff .lt. 0.000001) exit
                            alat1 = alat
                        end do
                        rlat(n) = h*alat*dpr
                        rlon(n) = along
                        if (rlon(n) .lt. 0.0) rlon(n) = rlon(n)+360.
                        if (rlon(n) .gt. 360.0) rlon(n) = rlon(n)-360.0
                        nret = nret+1
                        if (lrot) call polar_stereo_vect_rot(rlon(n), crot(n), srot(n))
                    else
                        rlon(n) = fill
                        rlat(n) = fill
                    end if
                end do
                !$omp end parallel do
            end if ! ELLIPTICAL
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
        elseif (iopt .eq. -1) then
            if (.not. elliptical) then
                !$omp parallel do private(n, dr, dr2) reduction(+:nret) schedule(static)
                do n = 1, npts
                    if (abs(rlon(n)) .lt. (360.+tinyreal) .and. abs(rlat(n)) .lt. (90.+tinyreal) .and. &
                        abs(h*rlat(n)+90) .gt. tinyreal) then
                        dr = de*tan((90-h*rlat(n))/2/dpr)
                        dr2 = dr**2
                        xpts(n) = xp+h*sin((rlon(n)-orient)/dpr)*dr/dxs
                        ypts(n) = yp-cos((rlon(n)-orient)/dpr)*dr/dys
                        if (xpts(n) .ge. xmin .and. xpts(n) .le. xmax .and. &
                            ypts(n) .ge. ymin .and. ypts(n) .le. ymax) then
                            nret = nret+1
                            if (lrot) call polar_stereo_vect_rot(rlon(n), crot(n), srot(n))
                            if (lmap) call polar_stereo_map_jacob(rlon(n), rlat(n), dr2, &
                                                                  xlon(n), xlat(n), ylon(n), ylat(n))
                            if (larea) call polar_stereo_grid_area(rlat(n), dr2, area(n))
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
            else  ! ELLIPTICAL CASE
                !$omp parallel do private(n, alat, along, t, rho) reduction(+:nret) schedule(static)
                do n = 1, npts
                    if (abs(rlon(n)) .lt. (360+tinyreal) .and. abs(rlat(n)) .lt. (90+tinyreal) .and. &
                        abs(h*rlat(n)+90) .gt. tinyreal) then
                        alat = h*rlat(n)/dpr
                        along = (rlon(n)-orient)/dpr
                        t = tan(pi4-alat*0.5)/((1.-e*sin(alat))/ &
                                               (1.+e*sin(alat)))**(e_over_2)
                        rho = rerth*mc*t/tc
                        xpts(n) = xp+rho*sin(h*along)/dxs
                        ypts(n) = yp-rho*cos(h*along)/dys
                        if (xpts(n) .ge. xmin .and. xpts(n) .le. xmax .and. &
                            ypts(n) .ge. ymin .and. ypts(n) .le. ymax) then
                            nret = nret+1
                            if (lrot) call polar_stereo_vect_rot(rlon(n), crot(n), srot(n))
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
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine gdswzd_polar_stereo

    !> Vector rotation fields for polar stereographic grids.
    !>
    !>
    !> This subprogram computes the vector rotation sines and
    !> cosines for a polar stereographic azimuthal grid.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "polar_stereo_vect_rot"
    !>
    !> @param[in] rlon grid point longitude in degrees (real)
    !> @param[in] crot clockwise vector rotation cosines (real)
    !> @param[in] srot clockwise vector rotation sines (real)
    !> (ugrid=crot*uearth-srot*vearth;
    !> vgrid=srot*uearth+crot*vearth)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine polar_stereo_vect_rot(rlon, crot, srot)
        implicit none

        real, intent(in) :: rlon
        real, intent(out) :: crot, srot

        if (irot .eq. 1) then
            crot = h*cos((rlon-orient)/dpr)
            srot = sin((rlon-orient)/dpr)
        else
            crot = 1.
            srot = 0.
        end if

    end subroutine polar_stereo_vect_rot

    !> Map jacobians for polar stereographic grids.
    !>
    !> This subprogram computes the map jacobians for
    !> a polar stereographic azimuthal grid (spherical
    !> earth).
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "polar_stereo_map_jacob"
    !> 2018-07-20 | wesley | pass in dr2 for threading.
    !>
    !> @param[in] rlon longitude in degrees (real)
    !> @param[in] rlat latitude in degrees (real)
    !> @param[in] dr2 squared distance from pole (real)
    !> @param[out] xlon dx/dlon in 1/degrees (real)
    !> @param[out] xlat dx/dlat in 1/degrees (real)
    !> @param[out] ylon dy/dlon in 1/degrees (real)
    !> @param[out] ylat dy/dlat in 1/degrees (real)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine polar_stereo_map_jacob(rlon, rlat, dr2, xlon, xlat, ylon, ylat)
        implicit none

        real, intent(in) :: rlon, rlat, dr2
        real, intent(out) :: xlon, xlat, ylon, ylat

        real                            :: clat, de, dr

        if (dr2 .lt. de2*1.e-6) then
            de = sqrt(de2)
            xlon = 0.
            xlat = -sin((rlon-orient)/dpr)/dpr*de/dxs/2
            ylon = 0.
            ylat = h*cos((rlon-orient)/dpr)/dpr*de/dys/2
        else
            dr = sqrt(dr2)
            clat = cos(rlat/dpr)
            xlon = h*cos((rlon-orient)/dpr)/dpr*dr/dxs
            xlat = -sin((rlon-orient)/dpr)/dpr*dr/dxs/clat
            ylon = sin((rlon-orient)/dpr)/dpr*dr/dys
            ylat = h*cos((rlon-orient)/dpr)/dpr*dr/dys/clat
        end if

    end subroutine polar_stereo_map_jacob

    !> Grid box area for polar stereographic grids.
    !>
    !> This subprogram computes the grid box area for
    !> a polar stereographic azimuthal grid (spherical
    !> earth).
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "polar_stereo_grid_area".
    !> 2018-07-20 | wesley | pass in dr2 for threading.
    !>
    !> @param[in] rlat latitude of grid point in degrees (real)
    !> @param[in] dr2 squared distance from pole (real)
    !> @param[out] area area weights in m**2 (real)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine polar_stereo_grid_area(rlat, dr2, area)
        implicit none

        real, intent(in) :: rlat, dr2
        real, intent(out) :: area

        real                            :: clat

        if (dr2 .lt. de2*1.e-6) then
            area = rerth**2*abs(dxs)*abs(dys)*4/de2
        else
            clat = cos(rlat/dpr)
            area = rerth**2*clat**2*abs(dxs)*abs(dys)/dr2
        end if

    end subroutine polar_stereo_grid_area

end module ip_polar_stereo_grid_mod
