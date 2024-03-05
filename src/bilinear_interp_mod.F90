!> @file
!> @brief Bilinear interpolation routines for scalars and vectors.
!> @author Mark Iredell, Kyle Gerheiser, Eric Engle

!> @brief Bilinear interpolation routines for scalars and vectors.
!>
!> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
module bilinear_interp_mod
    use gdswzd_mod
    use ip_grids_mod
    use ip_grid_descriptor_mod
    use ip_grid_factory_mod
    use polfix_mod
    implicit none

    private
    public :: interpolate_bilinear

    interface interpolate_bilinear
        module procedure interpolate_bilinear_scalar
        module procedure interpolate_bilinear_vector
    end interface interpolate_bilinear

    ! Smallest positive real value (use for equality comparisons)
    real :: tinyreal = tiny(1.0)

contains

    !> This subprogram performs bilinear interpolation
    !> from any grid to any grid for scalar fields.
    !>
    !> @details Options allow varying the minimum percentage for mask,
    !> i.e. percent valid input data required to make output data,
    !> (ipopt(1)) which defaults to 50 (if ipopt(1)=-1).
    !> only horizontal interpolation is performed.
    !> if no input data is found near the output point, a spiral
    !> search may be invoked by setting ipopt(2)> 0.
    !> no searching is done if output point is outside the input grid.
    !> as an added bonus the number of output grid points
    !> and their latitudes and longitudes are also returned.
    !> on the other hand, the output can be a set of station points
    !> if igdtnumo<0, in which case the number of points
    !> and their latitudes and longitudes must be input.
    !> input bitmaps will be interpolated to output bitmaps.
    !> output bitmaps will also be created when the output grid
    !> extends outside of the domain of the input grid.
    !>
    !> The output field is set to 0 where the output bitmap is off.
    !> @param[in] ipopt interpolation options
    !> - ipopt(1) is minimum percentage for mask (defaults to 50 if ipopt(1)=-1)
    !> - ipopt(2) is width of square to examine in spiral search (defaults to no search if ipopt(2)=-1)
    !> @param[in] grid_in input grid
    !> @param[in] grid_out output grid
    !> @param[in]  mi skip number between input grid fields if km>1 or dimension of input grid fields if km=1
    !> @param[out] mo skip number between output grid fields if km>1 or dimension of output grid fields if km=1
    !> @param[in]  km number of fields to interpolate
    !> @param[in]  ibi input bitmap flags
    !> @param[in]  li input bitmaps (if some ibi(k)=1)
    !> @param[in]  gi input fields to interpolate
    !> @param[in,out] no  number of output points (only if igdtnumo<0)
    !> @param[in,out] rlat output latitudes in degrees (if igdtnumo<0)
    !> @param[in,out] rlon output longitudes in degrees (if igdtnumo<0)
    !> @param[out] ibo output bitmap flags
    !> @param[out] lo output bitmaps (always output)
    !> @param[out] go output fields interpolated
    !> @param[out] iret return code
    !> - 0 successful interpolation
    !> - 2 unrecognized input grid or no grid overlap
    !> - 3 unrecognized output grid
    !>
    !> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
    subroutine interpolate_bilinear_scalar(ipopt, grid_in, grid_out, mi, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
        class(ip_grid), intent(in) :: grid_in, grid_out
        integer, intent(in) :: ipopt(20)
        integer, intent(in) :: mi, mo, km
        integer, intent(in) :: ibi(km)
        integer, intent(inout) :: no
        integer, intent(out) :: iret, ibo(km)
        !
        logical*1, intent(in) :: li(mi, km)
        logical*1, intent(out) :: lo(mo, km)
        !
        real, intent(in) :: gi(mi, km)
        real, intent(inout) :: rlat(mo), rlon(mo)
        real, intent(out) :: go(mo, km)
        !
        real, parameter     :: fill = -9999.
        !
        integer                              :: ijx(2), ijy(2)
        integer                              :: mp, n, i, j, k
        integer                              :: nk, nv
        integer                              :: mspiral, i1, j1, ixs, jxs
        integer                              :: mx, kxs, kxt, ix, jx, nx
        !
        logical                              :: same_gridi, same_grido
        !
        real                                 :: wx(2), wy(2)
        real                                 :: xpts(mo), ypts(mo)
        real                                 :: pmp, xij, yij, xf, yf, g, w

        logical :: to_station_points

        ! Save coeffecients between calls and only compute if grids have changed
        integer, save  :: nox = -1, iretx = -1
        integer, allocatable, save  :: nxy(:, :, :)
        real, allocatable, save  :: rlatx(:), rlonx(:)
        real, allocatable, save  :: wxy(:, :, :)
        class(ip_grid), allocatable, save :: prev_grid_in, prev_grid_out

        iret = 0
        mp = ipopt(1)
        if (mp .eq. -1 .or. mp .eq. 0) mp = 50
        if (mp .lt. 0 .or. mp .gt. 100) iret = 32
        pmp = mp*0.01
        mspiral = max(ipopt(2), 0)

        if (.not. allocated(prev_grid_in) .or. .not. allocated(prev_grid_out)) then
            allocate (prev_grid_in, source=grid_in)
            allocate (prev_grid_out, source=grid_out)

            same_gridi = .false.
            same_grido = .false.
        else
            same_gridi = grid_in .eq. prev_grid_in
            same_grido = grid_out .eq. prev_grid_out

            if (.not. same_gridi .or. .not. same_grido) then
                deallocate (prev_grid_in)
                deallocate (prev_grid_out)

                allocate (prev_grid_in, source=grid_in)
                allocate (prev_grid_out, source=grid_out)
            end if
        end if

        select type (grid_out)
        type is (ip_station_points_grid)
            to_station_points = .true.
        class default
            to_station_points = .false.
        end select

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SAVE OR SKIP WEIGHT COMPUTATION
        if (iret .eq. 0 .and. (to_station_points .or. .not. same_gridi .or. .not. same_grido)) then
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
            call gdswzd(grid_out, 0, mo, fill, xpts, ypts, rlon, rlat, no)
            if (no .eq. 0) iret = 3
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  LOCATE INPUT POINTS
            call gdswzd(grid_in, -1, no, fill, xpts, ypts, rlon, rlat, nv)
            if (iret .eq. 0 .and. nv .eq. 0) iret = 2
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  ALLOCATE AND SAVE GRID DATA
            if (nox .ne. no) then
                if (nox .ge. 0) deallocate (rlatx, rlonx, nxy, wxy)
                allocate (rlatx(no), rlonx(no), nxy(2, 2, no), wxy(2, 2, no))
                nox = no
            end if
            iretx = iret
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  COMPUTE WEIGHTS
            if (iret .eq. 0) then
                !$omp parallel do private(n, xij, yij, ijx, ijy, xf, yf, j, i, wx, wy) schedule(static)
                do n = 1, no
                    rlonx(n) = rlon(n)
                    rlatx(n) = rlat(n)
                    xij = xpts(n)
                    yij = ypts(n)
                    if (abs(xij-fill) .gt. tinyreal .and. abs(yij-fill) .gt. tinyreal) then
                        ijx(1:2) = floor(xij)+(/0, 1/)
                        ijy(1:2) = floor(yij)+(/0, 1/)
                        xf = xij-ijx(1)
                        yf = yij-ijy(1)
                        wx(1) = (1-xf)
                        wx(2) = xf
                        wy(1) = (1-yf)
                        wy(2) = yf
                        do j = 1, 2
                            do i = 1, 2
                                nxy(i, j, n) = grid_in%field_pos(ijx(i), ijy(j))
                                wxy(i, j, n) = wx(i)*wy(j)
                            end do
                        end do
                    else
                        nxy(:, :, n) = 0
                    end if
                end do
            end if
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  INTERPOLATE OVER ALL FIELDS
        if (iret .eq. 0 .and. iretx .eq. 0) then
            if (.not. to_station_points) then
                no = nox
                do n = 1, no
                    rlon(n) = rlonx(n)
                    rlat(n) = rlatx(n)
                end do
            end if
            !$omp parallel do &
                !$omp private(nk, k, n, g, w, j, i) &
                !$omp private(i1, j1, ixs, jxs, mx, kxs, kxt, ix, jx, nx) schedule(static)
            do nk = 1, no*km
                k = (nk-1)/no+1
                n = nk-no*(k-1)
                g = 0
                w = 0
                do j = 1, 2
                    do i = 1, 2
                        if (nxy(i, j, n) .gt. 0) then
                            if (ibi(k) .eq. 0 .or. li(nxy(i, j, n), k)) then
                                g = g+wxy(i, j, n)*gi(nxy(i, j, n), k)
                                w = w+wxy(i, j, n)
                            end if
                        end if
                    end do
                end do
                lo(n, k) = w .ge. pmp
                if (lo(n, k)) then
                    go(n, k) = g/w
                elseif (mspiral .gt. 0 .and. abs(xpts(n)-fill) .gt. tinyreal .and. abs(ypts(n)-fill) .gt. tinyreal) then
                    i1 = nint(xpts(n))
                    j1 = nint(ypts(n))
                    ixs = int(sign(1., xpts(n)-i1))
                    jxs = int(sign(1., ypts(n)-j1))
                    spiral: do mx = 1, mspiral**2
                        kxs = int(sqrt(4*mx-2.5))
                        kxt = mx-(kxs**2/4+1)
                        select case (mod(kxs, 4))
                        case (1)
                            ix = i1-ixs*(kxs/4-kxt)
                            jx = j1-jxs*kxs/4
                        case (2)
                            ix = i1+ixs*(1+kxs/4)
                            jx = j1-jxs*(kxs/4-kxt)
                        case (3)
                            ix = i1+ixs*(1+kxs/4-kxt)
                            jx = j1+jxs*(1+kxs/4)
                        case default
                            ix = i1-ixs*kxs/4
                            jx = j1+jxs*(kxs/4-kxt)
                        end select
                        nx = grid_in%field_pos(ix, jx)
                        if (nx .gt. 0.) then
                            if (li(nx, k) .or. ibi(k) .eq. 0) then
                                go(n, k) = gi(nx, k)
                                lo(n, k) = .true.
                                exit spiral
                            end if
                        end if
                    end do spiral
                    if (.not. lo(n, k)) then
                        ibo(k) = 1
                        go(n, k) = 0.
                    end if
                else
                    go(n, k) = 0.
                end if
            end do
            do k = 1, km
                ibo(k) = ibi(k)
                if (.not. all(lo(1:no, k))) ibo(k) = 1
            end do
            select type (grid_out)
            type is (ip_equid_cylind_grid)
                call polfixs(no, mo, km, rlat, ibo, lo, go)
            end select
        else
            if (iret .eq. 0) iret = iretx
            if (.not. to_station_points) no = 0
        end if

    end subroutine interpolate_bilinear_scalar

    !> This subprogram performs bilinear interpolation from any grid to
    !> any grid for vector fields.
    !>
    !> Options allow varying the minimum percentage for mask,
    !> i.e. percent valid input data required to make output data,
    !> (ipopt(1)) which defaults to 50 (if ipopt(1)=-1).
    !>
    !> Only horizontal interpolation is performed.
    !> the input and output vectors are rotated so that they are
    !> either resolved relative to the defined grid
    !> in the direction of increasing x and y coordinates
    !> or resolved relative to easterly and northerly directions,
    !> as designated by their respective grid description sections.
    !>
    !> As an added bonus the number of output grid points
    !> and their latitudes and longitudes are also returned
    !> along with their vector rotation parameters.
    !> on the other hand, the data may be interpolated to a set of
    !> station points if igdtnumo < 0, in which case the number
    !> of points and their latitudes and longitudes must be
    !> input along with their vector rotation parameters.
    !> input bitmaps will be interpolated to output bitmaps.
    !> output bitmaps will also be created when the output grid
    !> extends outside of the domain of the input grid.
    !> the output field is set to 0 where the output bitmap is off.
    !>
    !> @param[in] ipopt interpolation options
    !> - ipopt(1) is minimum percentage for mask (defaults to 50 if ipopt(1)=-1)
    !> @param[in] grid_in Input grid
    !> @param[in] grid_out Output grid
    !> @param[in]  mi skip number between input grid fields if km>1 or dimension of input grid fields if km=1
    !> @param[out] mo skip number between output grid fields if km>1 or dimension of output grid fields if km=1
    !> @param[in]  km number of fields to interpolate
    !> @param[in]  ibi input bitmap flags
    !> @param[in]  li input bitmaps (if some ibi(k)=1)
    !> @param[in]  ui input u-component fields to interpolate
    !> @param[in]  vi input v-component fields to interpolate
    !> @param[in,out] no  number of output points (only if igdtnumo<0)
    !> @param[in,out] rlat output latitudes in degrees (if igdtnumo<0)
    !> @param[in,out] rlon output longitudes in degrees (if igdtnumo<0)
    !> @param[in,out] crot vector rotation cosines (if igdtnumo<0) ugrid=crot*uearth-srot*vearth;
    !> @param[in,out] srot vector rotation sines (if igdtnumo<0) vgrid=srot*uearth+crot*vearth)
    !> @param[out] ibo output bitmap flags
    !> @param[out] lo output bitmaps (always output)
    !> @param[out] uo output u-component fields interpolated
    !> @param[out] vo output v-component fields interpolated
    !> @param[out] iret return code
    !> - 0 successful interpolation
    !> - 2 unrecognized input grid or no grid overlap
    !> - 3 unrecognized output grid
    !>
    !> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
    subroutine interpolate_bilinear_vector(ipopt, grid_in, grid_out, &
                                           mi, mo, km, ibi, li, ui, vi, &
                                           no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
        class(ip_grid), intent(in) :: grid_in, grid_out
        integer, intent(in) :: ipopt(20), ibi(km), mi, mo, km
        integer, intent(inout) :: no
        integer, intent(out) :: iret, ibo(km)
        !
        logical*1, intent(in) :: li(mi, km)
        logical*1, intent(out) :: lo(mo, km)
        !
        real, intent(in) :: ui(mi, km), vi(mi, km)
        real, intent(inout) :: rlat(mo), rlon(mo), crot(mo), srot(mo)
        real, intent(out) :: uo(mo, km), vo(mo, km)
        !
        real, parameter     :: fill = -9999.
        !
        integer                           :: ijx(2), ijy(2)
        integer                           :: mp, n, i, j, k, nk, nv
        !
        logical                           :: same_gridi, same_grido
        !
        real                              :: cm, sm, urot, vrot
        real                              :: pmp, xij, yij, xf, yf, u, v, w
        real                              :: xpts(mo), ypts(mo)
        real                              :: wx(2), wy(2)
        real                              :: xpti(mi), ypti(mi)
        real                              :: rloi(mi), rlai(mi)
        real                              :: croi(mi), sroi(mi)

        logical :: to_station_points

        ! Save coeffecients between calls and only compute if grids have changed
        integer, save  :: nox = -1, iretx = -1
        integer, allocatable, save  :: nxy(:, :, :)
        real, allocatable, save  :: rlatx(:), rlonx(:)
        real, allocatable, save  :: crotx(:), srotx(:)
        real, allocatable, save  :: wxy(:, :, :), cxy(:, :, :), sxy(:, :, :)
        class(ip_grid), allocatable, save :: prev_grid_in, prev_grid_out

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        iret = 0
        mp = ipopt(1)
        if (mp .eq. -1 .or. mp .eq. 0) mp = 50
        if (mp .lt. 0 .or. mp .gt. 100) iret = 32
        pmp = mp*0.01

        if (.not. allocated(prev_grid_in) .or. .not. allocated(prev_grid_out)) then
            allocate (prev_grid_in, source=grid_in)
            allocate (prev_grid_out, source=grid_out)

            same_gridi = .false.
            same_grido = .false.
        else
            same_gridi = grid_in .eq. prev_grid_in
            same_grido = grid_out .eq. prev_grid_out

            if (.not. same_gridi .or. .not. same_grido) then
                deallocate (prev_grid_in)
                deallocate (prev_grid_out)

                allocate (prev_grid_in, source=grid_in)
                allocate (prev_grid_out, source=grid_out)
            end if
        end if

        select type (grid_out)
        type is (ip_station_points_grid)
            to_station_points = .true.
        class default
            to_station_points = .false.
        end select

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SAVE OR SKIP WEIGHT COMPUTATION
        if (iret .eq. 0 .and. (to_station_points .or. .not. same_gridi .or. .not. same_grido)) then
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
            call gdswzd(grid_out, 0, mo, fill, xpts, ypts, rlon, rlat, no, crot, srot)
            if (no .eq. 0) iret = 3
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  LOCATE INPUT POINTS
            call gdswzd(grid_in, -1, no, fill, xpts, ypts, rlon, rlat, nv)
            if (iret .eq. 0 .and. nv .eq. 0) iret = 2
            call gdswzd(grid_in, 0, mi, fill, xpti, ypti, rloi, rlai, nv, croi, sroi)
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  ALLOCATE AND SAVE GRID DATA
            if (nox .ne. no) then
                if (nox .ge. 0) deallocate (rlatx, rlonx, crotx, srotx, nxy, wxy, cxy, sxy)
                allocate (rlatx(no), rlonx(no), crotx(no), srotx(no), &
                          nxy(2, 2, no), wxy(2, 2, no), cxy(2, 2, no), sxy(2, 2, no))
                nox = no
            end if
            iretx = iret
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  COMPUTE WEIGHTS
            if (iret .eq. 0) then
                !$omp parallel do private(n, xij, yij, ijx, ijy, xf, yf, j, i, wx, wy, cm, sm) schedule(static)
                do n = 1, no
                    rlonx(n) = rlon(n)
                    rlatx(n) = rlat(n)
                    crotx(n) = crot(n)
                    srotx(n) = srot(n)
                    xij = xpts(n)
                    yij = ypts(n)
                    if (abs(xij-fill) .gt. tinyreal .and. abs(yij-fill) .gt. tinyreal) then
                        ijx(1:2) = floor(xij)+(/0, 1/)
                        ijy(1:2) = floor(yij)+(/0, 1/)
                        xf = xij-ijx(1)
                        yf = yij-ijy(1)
                        wx(1) = (1-xf)
                        wx(2) = xf
                        wy(1) = (1-yf)
                        wy(2) = yf
                        do j = 1, 2
                            do i = 1, 2
                                nxy(i, j, n) = grid_in%field_pos(ijx(i), ijy(j))
                                wxy(i, j, n) = wx(i)*wy(j)
                                if (nxy(i, j, n) .gt. 0) then
                                    call movect(rlai(nxy(i, j, n)), rloi(nxy(i, j, n)), &
                                                rlat(n), rlon(n), cm, sm)
                                    cxy(i, j, n) = cm*croi(nxy(i, j, n))+sm*sroi(nxy(i, j, n))
                                    sxy(i, j, n) = sm*croi(nxy(i, j, n))-cm*sroi(nxy(i, j, n))
                                end if
                            end do
                        end do
                    else
                        nxy(:, :, n) = 0
                    end if
                end do
            end if  ! IS IRET 0?
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  INTERPOLATE OVER ALL FIELDS
        if (iret .eq. 0 .and. iretx .eq. 0) then
            if (.not. to_station_points) then
                no = nox
                do n = 1, no
                    rlon(n) = rlonx(n)
                    rlat(n) = rlatx(n)
                    crot(n) = crotx(n)
                    srot(n) = srotx(n)
                end do
            end if
            !$omp parallel do private(nk, k, n, u, v, w, urot, vrot, j, i) schedule(static)
            do nk = 1, no*km
                k = (nk-1)/no+1
                n = nk-no*(k-1)
                u = 0
                v = 0
                w = 0
                do j = 1, 2
                    do i = 1, 2
                        if (nxy(i, j, n) .gt. 0) then
                            if (ibi(k) .eq. 0 .or. li(nxy(i, j, n), k)) then
                                urot = cxy(i, j, n)*ui(nxy(i, j, n), k)-sxy(i, j, n)*vi(nxy(i, j, n), k)
                                vrot = sxy(i, j, n)*ui(nxy(i, j, n), k)+cxy(i, j, n)*vi(nxy(i, j, n), k)
                                u = u+wxy(i, j, n)*urot
                                v = v+wxy(i, j, n)*vrot
                                w = w+wxy(i, j, n)
                            end if
                        end if
                    end do
                end do
                lo(n, k) = w .ge. pmp
                if (lo(n, k)) then
                    urot = crot(n)*u-srot(n)*v
                    vrot = srot(n)*u+crot(n)*v
                    uo(n, k) = urot/w
                    vo(n, k) = vrot/w
                else
                    uo(n, k) = 0.
                    vo(n, k) = 0.
                end if
            end do  ! NK LOOP
            do k = 1, km
                ibo(k) = ibi(k)
                if (.not. all(lo(1:no, k))) ibo(k) = 1
            end do

            select type (grid_out)
            type is (ip_equid_cylind_grid)
                call polfixv(no, mo, km, rlat, rlon, ibo, lo, uo, vo)
            end select

        else
            if (iret .eq. 0) iret = iretx
            if (.not. to_station_points) no = 0
        end if
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine interpolate_bilinear_vector

end module bilinear_interp_mod
