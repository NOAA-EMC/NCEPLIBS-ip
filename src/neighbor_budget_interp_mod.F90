!> @file
!> @brief Interpolate scalar and vector fields with neighbor budget interpolation.
!> @author Mark Iredell @date 96-04-10

!> @brief Interpolate scalar fields (neighbor).
!>
!> @author Mark Iredell @date 96-04-10
module neighbor_budget_interp_mod
    use gdswzd_mod
    use polfix_mod
    use ip_grids_mod
    implicit none

    private
    public :: interpolate_neighbor_budget

    interface interpolate_neighbor_budget
        module procedure interpolate_neighbor_budget_scalar
        module procedure interpolate_neighbor_budget_vector
    end interface interpolate_neighbor_budget

    ! Smallest positive real value (use for equality comparisons)
    real :: tinyreal = tiny(1.0)

contains

    !> Interpolate scalar fields (budget).
    !>
    !> This subprogram performs budget interpolation from any grid to
    !> any grid for scalar fields.
    !>
    !> The algorithm simply computes (weighted) averages of neighbor
    !> points arranged in a square box centered around each output grid
    !> point and stretching nearly halfway to each of the neighboring
    !> grid points.
    !>
    !> Options allow choices of number of points in each radius from the
    !> center point (ipopt(1)) which defaults to 2 (if ipopt(1)=-1)
    !> meaning that 25 points will be averaged; further options are the
    !> respective weights for the radius points starting at the center
    !> point (ipopt(2:2+ipopt(1)) which defaults to all 1 (if
    !> ipopt(1)=-1 or ipopt(2)=-1).
    !>
    !> Another option is the minimum percentage for mask, i.e. percent
    !> valid input data required to make output data, (ipopt(3+ipopt(1))
    !> which defaults to 50 (if -1).
    !>
    !> Only horizontal interpolation is performed.
    !>
    !> The code recognizes the following projections, where "igdtnumi/o"
    !> is the grib 2 grid defintion template number for the input and
    !> output grids, respectively:
    !> - (igdtnumi/o=00) equidistant cylindrical
    !> - (igdtnumi/o=01) rotated equidistant cylindrical. "e" and non-"e" staggered
    !> - (igdtnumi/o=10) mercator cylindrical
    !> - (igdtnumi/o=20) polar stereographic azimuthal
    !> - (igdtnumi/o=30) lambert conformal conical
    !> - (igdtnumi/o=40) gaussian cylindrical
    !>
    !> As an added bonus the number of output grid points and their
    !> latitudes and longitudes are also returned. Input bitmaps will
    !> be interpolated to output bitmaps. Output bitmaps will also be
    !> created when the output grid extends outside of the domain of the
    !> input grid. The output field is set to 0 where the output bitmap
    !> is off.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | Iredell | Initial
    !> 96-10-04 | Iredell | neighbor points not bilinear interpolation
    !> 1999-04-08 | Iredell | split ijkgds into two pieces
    !> 2001-06-18 | Iredell | include minimum mask percentage option
    !> 2015-01-27 | Gayno | replace calls to gdswiz with new merged version of gdswzd.
    !> 2015-07-13 | Gayno | replace grib 1 kgds arrays with grib 2 grid definition template arrays.
    !> 2023-05-04 | Engle | allow calls to GDSWZD for station points
    !>
    !> @param[in] ipopt (20) interpolation options ipopt(1) is number of
    !> radius points (defaults to 2 if ipopt(1)=-1); ipopt(2:2+ipopt(1))
    !> are respective weights (defaults to all 1 if ipopt(1)=-1 or
    !> ipopt(2)=-1).  ipopt(3+ipopt(1)) is minimum percentage for mask
    !> (defaults to 50 if ipopt(3+ipopt(1)=-1)
    !> @param[in] grid_in The input grid.
    !> @param[in] grid_out The output grid.
    !> @param[in] mi skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1
    !> @param[in] mo skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1
    !> @param[in] km number of fields to interpolate
    !> @param[in] ibi (km) input bitmap flags
    !> @param[in] li (mi,km) input bitmaps (if some ibi(k)=1)
    !> @param[in] gi (mi,km) input fields to interpolate
    !> @param[out] no number of output points
    !> @param[out] rlat (mo) output latitudes in degrees
    !> @param[out] rlon (mo) output longitudes in degrees
    !> @param[out] ibo (km) output bitmap flags
    !> @param[out] lo (mo,km) output bitmaps (always output)
    !> @param[out] go (mo,km) output fields interpolated
    !> @param[out] iret return code
    !> - 0    successful interpolation
    !> - 2    unrecognized input grid or no grid overlap
    !> - 3    unrecognized output grid
    !> - 31   invalid undefined output grid
    !> - 32   invalid budget method parameters
    !>
    !> @author Mark Iredell @date 96-04-10
    subroutine interpolate_neighbor_budget_scalar(ipopt, grid_in, grid_out, &
                                                  mi, mo, km, ibi, li, gi, &
                                                  no, rlat, rlon, ibo, lo, go, iret)
        class(ip_grid), intent(in) :: grid_in, grid_out

        integer, intent(in) :: ibi(km), ipopt(20), km, mi, mo
        integer, intent(out) :: ibo(km), iret, no
        !
        logical*1, intent(in) :: li(mi, km)
        logical*1, intent(out) :: lo(mo, km)
        !
        real, intent(in) :: gi(mi, km)
        real, intent(out) :: go(mo, km), rlat(mo), rlon(mo)
        !
        real, parameter         :: fill = -9999.
        !
        integer                       :: ib, i1
        integer                       :: jb, j1, k, lb, lsw, mp, n
        integer                       :: n11(mo), nb, nb1, nb2, nb3, nb4, nv
        !
        real                          :: pmp, rlob(mo), rlab(mo)
        real                          :: wb, wo(mo, km), xi, yi
        real                          :: xptb(mo), yptb(mo), xpts(mo), ypts(mo)

        logical :: to_station_points

        select type (grid_out)
        type is (ip_station_points_grid)
            to_station_points = .true.
        class default
            to_station_points = .false.
        end select

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        iret = 0
        if (to_station_points) then
            call gdswzd(grid_out, 0, mo, fill, xpts, ypts, rlon, rlat, no)
            if (no .eq. 0) iret = 3
            call gdswzd(grid_in, -1, no, fill, xpts, ypts, rlon, rlat, nv)
            if (nv .eq. 0) iret = 2
        else
            call gdswzd(grid_out, 0, mo, fill, xpts, ypts, rlon, rlat, no)
            if (no .eq. 0) iret = 3
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        nb1 = ipopt(1)
        if (nb1 .eq. -1) nb1 = 2
        if (iret .eq. 0 .and. nb1 .lt. 0) iret = 32
        lsw = 1
        if (ipopt(1) .eq. -1 .or. ipopt(2) .eq. -1) lsw = 0
        if (iret .eq. 0 .and. lsw .eq. 1 .and. nb1 .gt. 15) iret = 32
        mp = ipopt(3+ipopt(1))
        if (mp .eq. -1 .or. mp .eq. 0) mp = 50
        if (mp .lt. 0 .or. mp .gt. 100) iret = 32
        pmp = mp*0.01
        if (iret .eq. 0) then
            nb2 = 2*nb1+1
            nb3 = nb2*nb2
            nb4 = nb3
            if (lsw .eq. 1) then
                nb4 = ipopt(2)
                do ib = 1, nb1
                    nb4 = nb4+8*ib*ipopt(2+ib)
                end do
            end if
        else
            nb2 = 0
            nb3 = 0
            nb4 = 0
        end if
        do k = 1, km
            do n = 1, no
                go(n, k) = 0.
                wo(n, k) = 0.
            end do
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX

        do nb = 1, nb3
            !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS
            jb = (nb-1)/nb2-nb1
            ib = nb-(jb+nb1)*nb2-nb1-1
            lb = max(abs(ib), abs(jb))
            wb = 1
            if (lsw .eq. 1) wb = ipopt(2+lb)
            if (abs(wb) .gt. tinyreal) then
                do n = 1, no
                    xptb(n) = xpts(n)+ib/real(nb2)
                    yptb(n) = ypts(n)+jb/real(nb2)
                end do
                if (to_station_points) then
                    call gdswzd(grid_in, 1, no, fill, xptb, yptb, rlob, rlab, nv)
                    call gdswzd(grid_in, -1, no, fill, xptb, yptb, rlob, rlab, nv)
                else
                    call gdswzd(grid_out, 1, no, fill, xptb, yptb, rlob, rlab, nv)
                    call gdswzd(grid_in, -1, no, fill, xptb, yptb, rlob, rlab, nv)
                end if
                if (iret .eq. 0 .and. nv .eq. 0 .and. lb .eq. 0) iret = 2
                do n = 1, no
                    xi = xptb(n)
                    yi = yptb(n)
                    if (abs(xi-fill) .gt. tinyreal .and. abs(yi-fill) .gt. tinyreal) then
                        i1 = nint(xi)
                        j1 = nint(yi)
                        n11(n) = grid_in%field_pos(i1, j1)
                    else
                        n11(n) = 0
                    end if
                end do
                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                !  INTERPOLATE WITH OR WITHOUT BITMAPS
                do k = 1, km
                    do n = 1, no
                        if (n11(n) .gt. 0) then
                            if (ibi(k) .eq. 0 .or. li(n11(n), k)) then
                                go(n, k) = go(n, k)+wb*gi(n11(n), k)
                                wo(n, k) = wo(n, k)+wb
                            end if
                        end if
                    end do
                end do
            end if
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE OUTPUT BITMAPS AND FIELDS
        do k = 1, km
            ibo(k) = ibi(k)
            do n = 1, no
                lo(n, k) = wo(n, k) .ge. pmp*nb4
                if (lo(n, k)) then
                    go(n, k) = go(n, k)/wo(n, k)
                else
                    ibo(k) = 1
                    go(n, k) = 0.
                end if
            end do
        end do

        select type (grid_out)
        type is (ip_equid_cylind_grid)
            call polfixs(no, mo, km, rlat, ibo, lo, go)
        end select

    end subroutine interpolate_neighbor_budget_scalar

    !> Interpolate vector fields (budget).
    !>
    !> This subprogram performs budget interpolation from any grid to
    !> any grid for vector fields.
    !>
    !> The algorithm simply computes (weighted) averages of neighbor
    !> points arranged in a square box centered around each output grid
    !> point and stretching nearly halfway to each of the neighboring
    !> grid points.
    !>
    !> Options allow choices of number of points in each radius from the
    !> center point (ipopt(1)) which defaults to 2 (if ipopt(1)=-1)
    !> meaning that 25 points will be averaged; further options are the
    !> respective weights for the radius points starting at the center
    !> point (ipopt(2:2+ipopt(1)) which defaults to all 1 (if
    !> ipopt(1)=-1 or ipopt(2)=-1).
    !>
    !> Another option is the minimum percentage for mask, i.e. percent
    !> valid input data required to make output data, (ipopt(3+ipopt(1))
    !> which defaults to 50 (if -1).
    !>
    !> Only horizontal interpolation is performed.
    !>
    !> The input and output grids are defined by their grib 2 grid
    !> definition template as decoded by the ncep g2 library. the code
    !> recognizes the following projections, where "igdtnumi/o" is the
    !> grib 2 grid defintion template number for the input and output
    !> grids, respectively:
    !> - (igdtnumi/o=00) equidistant cylindrical
    !> - (igdtnumi/o=01) rotated equidistant cylindrical. "e" and non-"e" staggered
    !> - (igdtnumi/o=10) mercator cylindrical
    !> - (igdtnumi/o=20) polar stereographic azimuthal
    !> - (igdtnumi/o=30) lambert conformal conical
    !> - (igdtnumi/o=40) gaussian cylindrical
    !>
    !> The input and output vectors are rotated so that they are either
    !> resolved relative to the defined grid in the direction of
    !> increasing x and y coordinates or resolved relative to easterly
    !> and northerly directions, as designated by their respective grid
    !> description sections.
    !>
    !> As an added bonus the number of output grid points and their
    !> latitudes and longitudes are also returned along with their
    !> vector rotation parameters. Input bitmaps will be interpolated
    !> to output bitmaps.
    !>
    !> Output bitmaps will also be created when the output grid extends
    !> outside of the domain of the input grid. The output field is set
    !> to 0 where the output bitmap is off.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | Iredell | Initial
    !> 1999-04-08 | Iredell | split ijkgds into two pieces
    !> 2001-06-18 | Iredell | include minimum mask percentage option
    !> 2002-01-17 | Iredell | save data from last call for optimization
    !> 2015-01-27 | Gayno | replace calls to gdswiz with new merged routine gdswzd.
    !> 2015-07-13 | Gayno | replace grib 1 kgds arrays with grib 2 grid definition template arrays.
    !> 2023-05-04 | Engle | allow calls to GDSWZD for station points
    !>
    !> @param[in] ipopt (20) interpolation options ipopt(1) is number of
    !> radius points (defaults to 2 if ipopt(1)=-1); ipopt(2:2+ipopt(1))
    !> are respective weights (defaults to all 1 if ipopt(1)=-1 or
    !> ipopt(2)=-1).  ipopt(3+ipopt(1)) is minimum percentage for mask
    !> (defaults to 50 if ipopt(3+ipopt(1)=-1)
    !> @param[in] grid_in The input grid.
    !> @param[in] grid_out The output grid.
    !> @param[in] mi skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1
    !> @param[in] mo skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1
    !> @param[in] km number of fields to interpolate
    !> @param[in] ibi (km) input bitmap flags
    !> @param[in] li (mi,km) input bitmaps (if some ibi(k)=1)
    !> @param[in] ui (mi,km) input u-component fields to interpolate
    !> @param[in] vi (mi,km) input v-component fields to interpolate
    !> @param[out] no number of output points
    !> @param[out] rlat (mo) output latitudes in degrees
    !> @param[out] rlon (mo) output longitudes in degrees
    !> @param[out] crot (mo) vector rotation cosines
    !> @param[out] srot (mo) vector rotation sines
    !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
    !> @param[out] ibo (km) output bitmap flags
    !> @param[out] lo (mo,km) output bitmaps (always output)
    !> @param[out] uo (mo,km) output u-component fields interpolated
    !> @param[out] vo (mo,km) output v-component fields interpolated
    !> @param[out] iret return code
    !> - 0    successful interpolation
    !> - 2    unrecognized input grid or no grid overlap
    !> - 3    unrecognized output grid
    !> - 31   invalid undefined output grid
    !> - 32   invalid budget method parameters
    !>
    !> @author Mark Iredell @date 96-04-10
    subroutine interpolate_neighbor_budget_vector(ipopt, grid_in, grid_out, &
                                                  mi, mo, km, ibi, li, ui, vi, &
                                                  no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
        class(ip_grid), intent(in) :: grid_in, grid_out

        integer, intent(in) :: ipopt(20), ibi(km)
        integer, intent(in) :: km, mi, mo
        integer, intent(out) :: iret, no, ibo(km)
        !
        logical*1, intent(in) :: li(mi, km)
        logical*1, intent(out) :: lo(mo, km)
        !
        real, intent(in) :: ui(mi, km), vi(mi, km)
        real, intent(inout) :: rlat(mo), rlon(mo)
        real, intent(out) :: uo(mo, km), vo(mo, km)
        real, intent(out) :: crot(mo), srot(mo)
        !
        real, parameter     :: fill = -9999.
        !
        integer                        :: n11(mo)
        integer                        :: ib, jb, i1, j1
        integer                        :: k, lb, lsw, mp, n, nv
        integer                        :: nb, nb1, nb2, nb3, nb4
        !
        logical                        :: same_grid
        !
        real                           :: c11(mo), s11(mo)
        real                           :: cm11, sm11, pmp
        real                           :: u11, v11, urot, vrot
        real                           :: wb, wo(mo, km), xi, yi
        real                           :: rlob(mo), rlab(mo)
        real                           :: xpts(mo), ypts(mo)
        real                           :: xptb(mo), yptb(mo)

        logical :: to_station_points

        ! Save coeffecients between runs and only compute if grid has changed
        integer, save  :: mix = -1
        real, allocatable, save  :: croi(:), sroi(:)
        real, allocatable, save  :: xpti(:), ypti(:)
        real, allocatable, save  :: rloi(:), rlai(:)
        class(ip_grid), allocatable, save :: prev_grid_in

        select type (grid_out)
        type is (ip_station_points_grid)
            to_station_points = .true.
        class default
            to_station_points = .false.
        end select

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        iret = 0
        call gdswzd(grid_out, 0, mo, fill, xpts, ypts, rlon, rlat, no, crot, srot)
        if (no .eq. 0) iret = 3
        if (to_station_points) then
            call gdswzd(grid_in, -1, no, fill, xpts, ypts, rlon, rlat, nv, crot, srot)
            if (nv .eq. 0) iret = 2
        end if

        if (.not. allocated(prev_grid_in)) then
            allocate (prev_grid_in, source=grid_in)

            same_grid = .false.
        else
            same_grid = grid_in .eq. prev_grid_in

            if (.not. same_grid) then
                deallocate (prev_grid_in)
                allocate (prev_grid_in, source=grid_in)
            end if
        end if

        if (.not. same_grid) then
            if (mix .ne. mi) then
                if (mix .ge. 0) deallocate (xpti, ypti, rloi, rlai, croi, sroi)
                allocate (xpti(mi), ypti(mi), rloi(mi), rlai(mi), croi(mi), sroi(mi))
                mix = mi
            end if
            call gdswzd(grid_in, 0, mi, fill, xpti, ypti, &
                        rloi, rlai, nv, croi, sroi)
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        nb1 = ipopt(1)
        if (nb1 .eq. -1) nb1 = 2
        if (iret .eq. 0 .and. nb1 .lt. 0) iret = 32
        lsw = 1
        if (ipopt(1) .eq. -1 .or. ipopt(2) .eq. -1) lsw = 0
        if (iret .eq. 0 .and. lsw .eq. 1 .and. nb1 .gt. 15) iret = 32
        mp = ipopt(3+ipopt(1))
        if (mp .eq. -1 .or. mp .eq. 0) mp = 50
        if (mp .lt. 0 .or. mp .gt. 100) iret = 32
        pmp = mp*0.01
        if (iret .eq. 0) then
            nb2 = 2*nb1+1
            nb3 = nb2*nb2
            nb4 = nb3
            if (lsw .eq. 1) then
                nb4 = ipopt(2)
                do ib = 1, nb1
                    nb4 = nb4+8*ib*ipopt(2+ib)
                end do
            end if
        else
            nb2 = 0
            nb3 = 0
            nb4 = 0
        end if
        do k = 1, km
            do n = 1, no
                uo(n, k) = 0
                vo(n, k) = 0
                wo(n, k) = 0.
            end do
        end do
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
        do nb = 1, nb3
            !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS AND ROTATIONS
            jb = (nb-1)/nb2-nb1
            ib = nb-(jb+nb1)*nb2-nb1-1
            lb = max(abs(ib), abs(jb))
            wb = 1
            if (lsw .eq. 1) wb = ipopt(2+lb)
            if (abs(wb) .gt. tinyreal) then
                do n = 1, no
                    xptb(n) = xpts(n)+ib/real(nb2)
                    yptb(n) = ypts(n)+jb/real(nb2)
                end do
                if (to_station_points) then
                    call gdswzd(grid_in, 1, no, fill, xptb, yptb, rlob, rlab, nv)
                    call gdswzd(grid_in, -1, no, fill, xptb, yptb, rlob, rlab, nv)
                else
                    call gdswzd(grid_out, 1, no, fill, xptb, yptb, rlob, rlab, nv)
                    call gdswzd(grid_in, -1, no, fill, xptb, yptb, rlob, rlab, nv)
                end if
                if (iret .eq. 0 .and. nv .eq. 0 .and. lb .eq. 0) iret = 2
                do n = 1, no
                    xi = xptb(n)
                    yi = yptb(n)
                    if (abs(xi-fill) .gt. tinyreal .and. abs(yi-fill) .gt. tinyreal) then
                        i1 = nint(xi)
                        j1 = nint(yi)
                        n11(n) = grid_in%field_pos(i1, j1)
                        if (n11(n) .gt. 0) then
                            call movect(rlai(n11(n)), rloi(n11(n)), rlat(n), rlon(n), cm11, sm11)
                            c11(n) = cm11*croi(n11(n))+sm11*sroi(n11(n))
                            s11(n) = sm11*croi(n11(n))-cm11*sroi(n11(n))
                        end if
                    else
                        n11(n) = 0
                    end if
                end do
                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                !  INTERPOLATE WITH OR WITHOUT BITMAPS
                do k = 1, km
                    do n = 1, no
                        if (n11(n) .gt. 0) then
                            if (ibi(k) .eq. 0 .or. li(n11(n), k)) then
                                u11 = c11(n)*ui(n11(n), k)-s11(n)*vi(n11(n), k)
                                v11 = s11(n)*ui(n11(n), k)+c11(n)*vi(n11(n), k)
                                uo(n, k) = uo(n, k)+wb*u11
                                vo(n, k) = vo(n, k)+wb*v11
                                wo(n, k) = wo(n, k)+wb
                            end if
                        end if
                    end do
                end do
            end if
        end do  ! NB LOOP
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE OUTPUT BITMAPS AND FIELDS
        do k = 1, km
            ibo(k) = ibi(k)
            do n = 1, no
                lo(n, k) = wo(n, k) .ge. pmp*nb4
                if (lo(n, k)) then
                    uo(n, k) = uo(n, k)/wo(n, k)
                    vo(n, k) = vo(n, k)/wo(n, k)
                    urot = crot(n)*uo(n, k)-srot(n)*vo(n, k)
                    vrot = srot(n)*uo(n, k)+crot(n)*vo(n, k)
                    uo(n, k) = urot
                    vo(n, k) = vrot
                else
                    ibo(k) = 1
                    uo(n, k) = 0.
                    vo(n, k) = 0.
                end if
            end do
        end do

        select type (grid_out)
        type is (ip_equid_cylind_grid)
            call polfixv(no, mo, km, rlat, rlon, ibo, lo, uo, vo)
        end select

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine interpolate_neighbor_budget_vector

end module neighbor_budget_interp_mod
