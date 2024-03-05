!> @file
!! @brief Budget interpolation routines for scalars and vectors
!! @author Mark Iredell, Kyle Gerheiser, Eric Engle
!! @date July 2021

!> @brief Budget interpolation routines for scalars and vectors.
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
!! @date July 2021
module budget_interp_mod
    use gdswzd_mod
    use polfix_mod
    use ip_grids_mod
    use ip_grid_descriptor_mod
    use ip_grid_factory_mod
    implicit none

    private
    public :: interpolate_budget

    interface interpolate_budget
        module procedure interpolate_budget_scalar
        module procedure interpolate_budget_vector
    endinterface interpolate_budget

    ! Smallest positive real value (use for equality comparisons)
    real :: tinyreal=tiny(1.0)

contains

    !> Performs budget interpolation from any grid to any grid (or to
    !> random station points) for scalar fields.
    !>
    !> The algorithm simply computes (weighted) averages of bilinearly
    !> interpolated points arranged in a square box centered around each
    !> output grid point and stretching nearly halfway to each of the
    !> neighboring grid points.
    !>
    !> Options allow choices of number of points in each radius from the
    !> center point (ipopt(1)) which defaults to 2 (if ipopt(1)=-1)
    !> meaning that 25 points will be averaged; further options are the
    !> respective weights for the radius points starting at the center
    !> point (ipopt(2:2+ipopt(1)) which defaults to all 1 (if
    !> ipopt(1)=-1 or ipopt(2)=-1).
    !>
    !> A special interpolation is done if ipopt(2)=-2.  in this case,
    !> the boxes stretch nearly all the way to each of the neighboring
    !> grid points and the weights are the adjoint of the bilinear
    !> interpolation weights.  This case gives quasi-second-order budget
    !> interpolation.
    !>
    !> Another option is the minimum percentage for mask, i.e. percent
    !> valid input data required to make output data, (ipopt(3+ipopt(1))
    !> which defaults to 50 (if -1).
    !>
    !> In cases where there is no or insufficient valid input data, the
    !> user may choose to search for the nearest valid data.  this is
    !> invoked by setting ipopt(20) to the width of the search
    !> square. The default is 1 (no search). Squares are searched for
    !> valid data in a spiral pattern starting from the center. No
    !> searching is done where the output grid is outside the input
    !> grid.
    !>
    !> Only horizontal interpolation is performed.
    !>
    !> @param[in] ipopt Interpolation options
    !> - ipopt(1) is number of radius points (defaults to 2 if ipopt(1)=-1).
    !> - ipopt(2:2+ipopt(1)) are respective weights (defaults to all 1 if ipopt(1)=-1 or ipopt(2)=-1).
    !> - ipopt(3+ipopt(1)) is minimum percentage for mask (defaults to 50 if ipopt(3+ipopt(1)=-1).
    !> @param[in] grid_in Input grid
    !> @param[in] grid_out Output grid
    !> @param[in] mi Skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1.
    !> @param[out] mo Skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1.
    !> @param[in] km Number of fields to interpolate.
    !> @param[in] ibi Input bitmap flags.
    !> @param[in] li Input bitmaps (if some ibi(k)=1).
    !> @param[in] gi Input fields to interpolate.
    !> @param[in,out] no  Number of output points (only if igdtnumo<0).
    !> @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0).
    !> @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0).
    !> @param[out] ibo Output bitmap flags.
    !> @param[out] lo Output bitmaps (always output).
    !> @param[out] go Output fields interpolated.
    !> @param[out] iret Return code.
    !> - 0 Successful interpolation.
    !> - 2 Unrecognized input grid or no grid overlap.
    !> - 3 Unrecognized output grid.
    !> - 32 Invalid budget method parameters.
    !>
    !> @author Marke Iredell, George Gayno, Kyle Gerheiser, Eric Engle
    !> @date July 2021
    subroutine interpolate_budget_scalar(ipopt,grid_in,grid_out, &
                                         mi,mo,km,ibi,li,gi, &
                                         no,rlat,rlon,ibo,lo,go,iret)
        class(ip_grid),intent(in) :: grid_in,grid_out
        integer,intent(in)     :: ibi(km),ipopt(20)
        integer,intent(in)     :: km,mi,mo
        integer,intent(out)     :: ibo(km),iret,no
        !
        logical*1,intent(in)     :: li(mi,km)
        logical*1,intent(out)     :: lo(mo,km)
        !
        real,intent(in)     :: gi(mi,km)
        real,intent(inout)     :: rlat(mo),rlon(mo)
        real,intent(out)     :: go(mo,km)
        !
        real,parameter         :: fill=-9999.
        !
        integer                       :: i1,j1,i2,j2,ib,jb
        integer                       :: ix,jx,ixs,jxs
        integer                       :: k,kxs,kxt
        integer                       :: lb,lsw,mp,mspiral,mx
        integer                       :: n,nb,nb1,nb2,nb3,nb4,nv,nx
        integer                       :: n11(mo),n21(mo),n12(mo),n22(mo)
        !
        real                          :: gb,lat(1),lon(1)
        real                          :: pmp,rb2,rlob(mo),rlab(mo),wb
        real                          :: w11(mo),w21(mo),w12(mo),w22(mo)
        real                          :: wo(mo,km),xf,yf,xi,yi,xx,yy
        real                          :: xpts(mo),ypts(mo),xptb(mo),yptb(mo)
        real                          :: xxx(1),yyy(1)

        logical :: to_station_points

        class(ip_grid),allocatable :: grid_out2

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        !  DO SUBSECTION OF GRID IF KGDSO(1) IS SUBTRACTED FROM 255.
        iret=0

        select type(grid_out)
        type is(ip_station_points_grid)
            to_station_points=.true.
            allocate(grid_out2,source=grid_out)
            call gdswzd(grid_out2,0,mo,fill,xpts,ypts,rlon,rlat,no)
            if(no.eq.0) iret=3
            call gdswzd(grid_in,-1,no,fill,xpts,ypts,rlon,rlat,nv)
            if(nv.eq.0) iret=2
        class default
            to_station_points=.false.
            allocate(grid_out2,source=grid_out)
            call gdswzd(grid_out2,0,mo,fill,xpts,ypts,rlon,rlat,no)
            if(no.eq.0) iret=3
        endselect

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        if(ipopt(1).gt.16) iret=32
        mspiral=max(ipopt(20),1)
        nb1=ipopt(1)
        if(nb1.eq.-1) nb1=2
        if(iret.eq.0.and.nb1.lt.0) iret=32
        lsw=1
        if(ipopt(2).eq.-2) lsw=2
        if(ipopt(1).eq.-1.or.ipopt(2).eq.-1) lsw=0
        if(iret.eq.0.and.lsw.eq.1.and.nb1.gt.15) iret=32
        mp=ipopt(3+ipopt(1))
        if(mp.eq.-1.or.mp.eq.0) mp=50
        if(mp.lt.0.or.mp.gt.100) iret=32
        pmp=mp*0.01
        if(iret.eq.0) then
            nb2=2*nb1+1
            rb2=1./nb2
            nb3=nb2*nb2
            nb4=nb3
            if(lsw.eq.2) then
                rb2=1./(nb1+1)
                nb4=(nb1+1)**4
            elseif(lsw.eq.1) then
                nb4=ipopt(2)
                do ib=1,nb1
                    nb4=nb4+8*ib*ipopt(2+ib)
                enddo
            endif
        else
            nb3=0
            nb4=1
        endif
        do k=1,km
            do n=1,no
                go(n,k)=0.
                wo(n,k)=0.
            enddo
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
        do nb=1,nb3
            !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS
            jb=(nb-1)/nb2-nb1
            ib=nb-(jb+nb1)*nb2-nb1-1
            lb=max(abs(ib),abs(jb))
            wb=1
            if(lsw.eq.2) then
                wb=(nb1+1-abs(ib))*(nb1+1-abs(jb))
            elseif(lsw.eq.1) then
                wb=ipopt(2+lb)
            endif
            if(abs(wb).gt.tinyreal) then
                !$omp parallel do private(n) schedule(static)
                do n=1,no
                    xptb(n)=xpts(n)+ib*rb2
                    yptb(n)=ypts(n)+jb*rb2
                enddo
                !$omp end parallel do
                if(to_station_points) then
                    call gdswzd(grid_in,1,no,fill,xptb,yptb,rlob,rlab,nv)
                    call gdswzd(grid_in,-1,no,fill,xptb,yptb,rlob,rlab,nv)
                else
                    call gdswzd(grid_out2,1,no,fill,xptb,yptb,rlob,rlab,nv)
                    call gdswzd(grid_in,-1,no,fill,xptb,yptb,rlob,rlab,nv)
                endif
                if(iret.eq.0.and.nv.eq.0.and.lb.eq.0) iret=2
                !$omp parallel do private(n,xi,yi,i1,i2,j1,j2,xf,yf) schedule(static)
                do n=1,no
                    xi=xptb(n)
                    yi=yptb(n)
                    if(abs(xi-fill).gt.tinyreal.and.abs(yi-fill).gt.tinyreal) then
                        i1=int(xi)
                        i2=i1+1
                        j1=int(yi)
                        j2=j1+1
                        xf=xi-i1
                        yf=yi-j1
                        n11(n)=grid_in%field_pos(i1,j1)
                        n21(n)=grid_in%field_pos(i2,j1)
                        n12(n)=grid_in%field_pos(i1,j2)
                        n22(n)=grid_in%field_pos(i2,j2)
                        if(min(n11(n),n21(n),n12(n),n22(n)).gt.0) then
                            w11(n)=(1-xf)*(1-yf)
                            w21(n)=xf*(1-yf)
                            w12(n)=(1-xf)*yf
                            w22(n)=xf*yf
                        else
                            n11(n)=0
                            n21(n)=0
                            n12(n)=0
                            n22(n)=0
                        endif
                    else
                        n11(n)=0
                        n21(n)=0
                        n12(n)=0
                        n22(n)=0
                    endif
                enddo
                !$omp end parallel do
                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                !  INTERPOLATE WITH OR WITHOUT BITMAPS
                !$omp parallel do private(k,n,gb) schedule(static)
                do k=1,km
                    do n=1,no
                        if(n11(n).gt.0) then
                            if(ibi(k).eq.0) then
                                gb=w11(n)*gi(n11(n),k)+w21(n)*gi(n21(n),k) &
                                    +w12(n)*gi(n12(n),k)+w22(n)*gi(n22(n),k)
                                go(n,k)=go(n,k)+wb*gb
                                wo(n,k)=wo(n,k)+wb
                            else
                                if(li(n11(n),k)) then
                                    go(n,k)=go(n,k)+wb*w11(n)*gi(n11(n),k)
                                    wo(n,k)=wo(n,k)+wb*w11(n)
                                endif
                                if(li(n21(n),k)) then
                                    go(n,k)=go(n,k)+wb*w21(n)*gi(n21(n),k)
                                    wo(n,k)=wo(n,k)+wb*w21(n)
                                endif
                                if(li(n12(n),k)) then
                                    go(n,k)=go(n,k)+wb*w12(n)*gi(n12(n),k)
                                    wo(n,k)=wo(n,k)+wb*w12(n)
                                endif
                                if(li(n22(n),k)) then
                                    go(n,k)=go(n,k)+wb*w22(n)*gi(n22(n),k)
                                    wo(n,k)=wo(n,k)+wb*w22(n)
                                endif
                            endif
                        endif
                    enddo
                enddo
                !$omp end parallel do
            endif
        enddo   ! sub-grid points
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE OUTPUT BITMAPS AND FIELDS
        ! KM is often 1 .. do not do OMP PARALLEL DO here
        km_loop: do k=1,km
            ibo(k)=ibi(k)
            !$omp parallel do private(n,lat,lon,xxx,yyy,nv,xx,yy,ixs,jxs,mx,kxs,kxt,ix,jx,nx) schedule(static)
            n_loop: do n=1,no
                lo(n,k)=wo(n,k).ge.pmp*nb4
                if(lo(n,k)) then
                    go(n,k)=go(n,k)/wo(n,k)
                elseif(mspiral.gt.1) then
                    lat(1)=rlat(n)
                    lon(1)=rlon(n)
                    call gdswzd(grid_in,-1,1,fill,xxx,yyy,lon,lat,nv)
                    xx=xxx(1)
                    yy=yyy(1)
                    if(nv.eq.1) then
                        i1=nint(xx)
                        j1=nint(yy)
                        ixs=int(sign(1.,xx-i1))
                        jxs=int(sign(1.,yy-j1))
                        spiral_loop: do mx=2,mspiral**2
                            kxs=int(sqrt(4*mx-2.5))
                            kxt=mx-(kxs**2/4+1)
                            select case(mod(kxs,4))
                            case(1)
                                ix=i1-ixs*(kxs/4-kxt)
                                jx=j1-jxs*kxs/4
                            case(2)
                                ix=i1+ixs*(1+kxs/4)
                                jx=j1-jxs*(kxs/4-kxt)
                            case(3)
                                ix=i1+ixs*(1+kxs/4-kxt)
                                jx=j1+jxs*(1+kxs/4)
                            case default
                                ix=i1-ixs*kxs/4
                                jx=j1+jxs*(kxs/4-kxt)
                            endselect
                            nx=grid_in%field_pos(ix,jx)
                            if(nx.gt.0.) then
                                if(li(nx,k).or.ibi(k).eq.0) then
                                    go(n,k)=gi(nx,k)
                                    lo(n,k)=.true.
                                    cycle n_loop
                                endif
                            endif
                        enddo spiral_loop
                        ibo(k)=1
                        go(n,k)=0.
                    else
                        ibo(k)=1
                        go(n,k)=0.
                    endif
                else  ! no spiral search option
                    ibo(k)=1
                    go(n,k)=0.
                endif
            enddo n_loop
            !$omp end parallel do
        enddo km_loop

        select type(grid_out2)
        type is(ip_equid_cylind_grid)
            call polfixs(no,mo,km,rlat,ibo,lo,go)
        endselect
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    endsubroutine interpolate_budget_scalar

    !> This subprogram performs budget interpolation from any grid to
    !> any grid (or to random station points) for vector fields.
    !>
    !> The algorithm simply computes (weighted) averages of bilinearly
    !> interpolated points arranged in a square box centered around each
    !> output grid point and stretching nearly halfway to each of the
    !> neighboring grid points.
    !>
    !> Options allow choices of number of points in each radius from the
    !> center point (ipopt(1)) which defaults to 2 (if ipopt(1)=-1)
    !> meaning that 25 points will be averaged; further options are the
    !> respective weights for the radius points starting at the center
    !> point (ipopt(2:2+ipopt(1)) which defaults to all 1 (if
    !> ipopt(1)=-1 or ipopt(2)=-1).
    !>
    !> A special interpolation is done if ipopt(2)=-2.  in this case,
    !> the boxes stretch nearly all the way to each of the neighboring
    !> grid points and the weights are the adjoint of the bilinear
    !> interpolation weights.  This case gives quasi-second-order budget
    !> interpolation.
    !>
    !> Another option is the minimum percentage for mask, i.e. percent
    !> valid input data required to make output data, (ipopt(3+ipopt(1))
    !> which defaults to 50 (if -1).
    !>
    !> In cases where there is no or insufficient valid input data, the
    !> user may choose to search for the nearest valid data.  this is
    !> invoked by setting ipopt(20) to the width of the search
    !> square. The default is 1 (no search). Squares are searched for
    !> valid data in a spiral pattern starting from the center. No
    !> searching is done where the output grid is outside the input
    !> grid.
    !>
    !> Only horizontal interpolation is performed.
    !>
    !> @param[in] ipopt interpolation options ipopt(1) Number of radius
    !> points (defaults to 2 if ipopt(1)=-1); ipopt(2:2+ipopt(1))
    !> Respective weights (defaults to all 1 if ipopt(1)=-1 or
    !> ipopt(2)=-1).  ipopt(3+ipopt(1)) Minimum percentage for mask
    !> (defaults to 50 if ipopt(3+ipopt(1)=-1)
    !> @param[in] grid_in Input grid.
    !> @param[in] grid_out Output grid.
    !> @param[in] mi skip Number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1.
    !> @param[out] mo skip Number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1.
    !> @param[in] km Number of fields to interpolate.
    !> @param[in] ibi Input bitmap flags.
    !> @param[in] li Input bitmaps (if some ibi(k)=1).
    !> @param[in] ui Input u-component fields to interpolate.
    !> @param[in] vi Input v-component fields to interpolate.
    !> @param[in,out] no  Number of output points (only if igdtnumo<0)
    !> @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0)
    !> @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0)
    !> @param[in,out] crot Vector rotation cosines. If interpolating
    !> subgrid ugrid=crot * uearth - srot * vearth.
    !> @param[in,out] srot Vector rotation sines. If interpolating
    !> subgrid vgrid = srot * uearth + crot * vearth.
    !> @param[out] ibo Output bitmap flags.
    !> @param[out] lo Output bitmaps (always output).
    !> @param[out] uo Output u-component fields interpolated.
    !> @param[out] vo Output v-component fields interpolated.
    !> @param[out] iret Return code.
    !> - 0 Successful interpolation.
    !> - 2 Unrecognized input grid or no grid overlap.
    !> - 3 Unrecognized output grid.
    !> - 32 Invalid budget method parameters.
    !>
    !> @author Marke Iredell, George Gayno, Kyle Gerheiser, Eric Engle
    !> @date July 2021
    subroutine interpolate_budget_vector(ipopt,grid_in,grid_out, &
                                         mi,mo,km,ibi,li,ui,vi, &
                                         no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
        class(ip_grid),intent(in) :: grid_in,grid_out
        integer,intent(in) :: ipopt(20),ibi(km)
        integer,intent(in) :: km,mi,mo
        integer,intent(out) :: iret,no,ibo(km)
        !
        logical*1,intent(in) :: li(mi,km)
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: ui(mi,km),vi(mi,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: uo(mo,km),vo(mo,km)
        real,intent(out) :: crot(mo),srot(mo)
        !
        real,parameter     :: fill=-9999.
        !
        integer                         :: i1,i2,j1,j2,ib,jb,lsw,mp
        integer                         :: k,lb,n,nb,nb1,nb2,nb3,nb4,nv
        integer                         :: n11(mo),n21(mo),n12(mo),n22(mo)
        !
        logical                         :: same_grid
        !
        real                            :: cm11,sm11,cm12,sm12
        real                            :: cm21,sm21,cm22,sm22
        real                            :: pmp,rb2
        real                            :: c11(mo),c21(mo),c12(mo),c22(mo)
        real                            :: s11(mo),s21(mo),s12(mo),s22(mo)
        real                            :: w11(mo),w21(mo),w12(mo),w22(mo)
        real                            :: ub,vb,wb,urot,vrot
        real                            :: u11,v11,u21,v21,u12,v12,u22,v22
        real                            :: wi1,wj1,wi2,wj2
        real                            :: wo(mo,km),xi,yi
        real                            :: xpts(mo),ypts(mo)
        real                            :: xptb(mo),yptb(mo),rlob(mo),rlab(mo)

        logical :: to_station_points

        class(ip_grid),allocatable :: grid_out2

        ! Save coeffecients between calls and only compute if grids have changed
        integer,save          :: mix=-1
        real,allocatable,save :: croi(:),sroi(:)
        real,allocatable,save :: xpti(:),ypti(:),rloi(:),rlai(:)

        class(ip_grid),allocatable,save :: prev_grid_in

        iret=0

        ! Negative grid number means interpolate to subgrid
        ! The type of the subgrid is calculated by 255 +
        select type(grid_out)
        type is(ip_station_points_grid)
            to_station_points=.true.
            allocate(grid_out2,source=grid_out)
            call gdswzd(grid_out2,0,mo,fill,xpts,ypts,rlon,rlat,no,crot,srot)
            if(no.eq.0) iret=3
            call gdswzd(grid_in,-1,no,fill,xpts,ypts,rlon,rlat,nv,crot,srot)
            if(nv.eq.0) iret=2
        class default
            to_station_points=.false.
            allocate(grid_out2,source=grid_out)
            call gdswzd(grid_out2,0,mo,fill,xpts,ypts,rlon,rlat,no,crot,srot)
        endselect

        if(.not.allocated(prev_grid_in)) then
            allocate(prev_grid_in,source=grid_in)

            same_grid=.false.
        else
            same_grid=grid_in.eq.prev_grid_in

            if(.not.same_grid) then
                deallocate(prev_grid_in)
                allocate(prev_grid_in,source=grid_in)
            endif
        endif

        if(.not.same_grid) then
            if(mix.ne.mi) then
                if(mix.ge.0) deallocate(xpti,ypti,rloi,rlai,croi,sroi)
                allocate(xpti(mi),ypti(mi),rloi(mi),rlai(mi),croi(mi),sroi(mi))
                mix=mi
            endif
            call gdswzd(grid_in,0,mi,fill,xpti,ypti,rloi,rlai,nv,croi,sroi)
        endif

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        nb1=ipopt(1)
        if(nb1.eq.-1) nb1=2
        if(iret.eq.0.and.nb1.lt.0) iret=32
        lsw=1
        if(ipopt(2).eq.-2) lsw=2
        if(ipopt(1).eq.-1.or.ipopt(2).eq.-1) lsw=0
        if(iret.eq.0.and.lsw.eq.1.and.nb1.gt.15) iret=32
        mp=ipopt(3+ipopt(1))
        if(mp.eq.-1.or.mp.eq.0) mp=50
        if(mp.lt.0.or.mp.gt.100) iret=32
        pmp=mp*0.01
        if(iret.eq.0) then
            nb2=2*nb1+1
            rb2=1./nb2
            nb3=nb2*nb2
            nb4=nb3
            if(lsw.eq.2) then
                rb2=1./(nb1+1)
                nb4=(nb1+1)**4
            elseif(lsw.eq.1) then
                nb4=ipopt(2)
                do ib=1,nb1
                    nb4=nb4+8*ib*ipopt(2+ib)
                enddo
            endif
        else
            nb3=0
            nb4=1
        endif
        do k=1,km
            do n=1,no
                uo(n,k)=0
                vo(n,k)=0
                wo(n,k)=0.
            enddo
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  LOOP OVER SAMPLE POINTS IN OUTPUT GRID BOX
        do nb=1,nb3
            !  LOCATE INPUT POINTS AND COMPUTE THEIR WEIGHTS AND ROTATIONS
            jb=(nb-1)/nb2-nb1
            ib=nb-(jb+nb1)*nb2-nb1-1
            lb=max(abs(ib),abs(jb))
            wb=1
            if(ipopt(2).eq.-2) then
                wb=(nb1+1-abs(ib))*(nb1+1-abs(jb))
            elseif(ipopt(2).ne.-1) then
                wb=ipopt(2+lb)
            endif
            if(abs(wb).gt.tinyreal) then
                !$omp parallel do private(n) schedule(static)
                do n=1,no
                    xptb(n)=xpts(n)+ib*rb2
                    yptb(n)=ypts(n)+jb*rb2
                enddo
                !$omp end parallel do
                if(to_station_points) then
                    call gdswzd(grid_in,1,no,fill,xptb,yptb,rlob,rlab,nv)
                    call gdswzd(grid_in,-1,no,fill,xptb,yptb,rlob,rlab,nv)
                else
                    call gdswzd(grid_out2,1,no,fill,xptb,yptb,rlob,rlab,nv)
                    call gdswzd(grid_in,-1,no,fill,xptb,yptb,rlob,rlab,nv)
                endif
                if(iret.eq.0.and.nv.eq.0.and.lb.eq.0) iret=2
                !$omp parallel do private(n,xi,yi,i1,i2,wi1,wi2,j1,j2,wj1,wj2,cm11,cm21,cm12,cm22,sm11,sm21,sm12,sm22) &
                    !$omp schedule(static)
                do n=1,no
                    xi=xptb(n)
                    yi=yptb(n)
                    if(abs(xi-fill).gt.tinyreal.and.abs(yi-fill).gt.tinyreal) then
                        i1=int(xi)
                        i2=i1+1
                        wi2=xi-i1
                        wi1=1-wi2
                        j1=int(yi)
                        j2=j1+1
                        wj2=yi-j1
                        wj1=1-wj2
                        n11(n)=grid_in%field_pos(i1,j1)
                        n21(n)=grid_in%field_pos(i2,j1)
                        n12(n)=grid_in%field_pos(i1,j2)
                        n22(n)=grid_in%field_pos(i2,j2)
                        if(min(n11(n),n21(n),n12(n),n22(n)).gt.0) then
                            w11(n)=wi1*wj1
                            w21(n)=wi2*wj1
                            w12(n)=wi1*wj2
                            w22(n)=wi2*wj2
                            call movect(rlai(n11(n)),rloi(n11(n)),rlat(n),rlon(n),cm11,sm11)
                            call movect(rlai(n21(n)),rloi(n21(n)),rlat(n),rlon(n),cm21,sm21)
                            call movect(rlai(n12(n)),rloi(n12(n)),rlat(n),rlon(n),cm12,sm12)
                            call movect(rlai(n22(n)),rloi(n22(n)),rlat(n),rlon(n),cm22,sm22)
                            c11(n)=cm11*croi(n11(n))+sm11*sroi(n11(n))
                            s11(n)=sm11*croi(n11(n))-cm11*sroi(n11(n))
                            c21(n)=cm21*croi(n21(n))+sm21*sroi(n21(n))
                            s21(n)=sm21*croi(n21(n))-cm21*sroi(n21(n))
                            c12(n)=cm12*croi(n12(n))+sm12*sroi(n12(n))
                            s12(n)=sm12*croi(n12(n))-cm12*sroi(n12(n))
                            c22(n)=cm22*croi(n22(n))+sm22*sroi(n22(n))
                            s22(n)=sm22*croi(n22(n))-cm22*sroi(n22(n))
                        else
                            n11(n)=0
                            n21(n)=0
                            n12(n)=0
                            n22(n)=0
                        endif
                    else
                        n11(n)=0
                        n21(n)=0
                        n12(n)=0
                        n22(n)=0
                    endif
                enddo
                !$omp end parallel do
                ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                !  INTERPOLATE WITH OR WITHOUT BITMAPS
                !  KM IS OFTEN 1 .. DO NO PUT OMP PARALLEL DO HERE
                do k=1,km
                    !$omp parallel do private(n,u11,u12,u21,u22,ub,v11,v12,v21,v22,vb) schedule(static)
                    do n=1,no
                        if(n11(n).gt.0) then
                            if(ibi(k).eq.0) then
                                u11=c11(n)*ui(n11(n),k)-s11(n)*vi(n11(n),k)
                                v11=s11(n)*ui(n11(n),k)+c11(n)*vi(n11(n),k)
                                u21=c21(n)*ui(n21(n),k)-s21(n)*vi(n21(n),k)
                                v21=s21(n)*ui(n21(n),k)+c21(n)*vi(n21(n),k)
                                u12=c12(n)*ui(n12(n),k)-s12(n)*vi(n12(n),k)
                                v12=s12(n)*ui(n12(n),k)+c12(n)*vi(n12(n),k)
                                u22=c22(n)*ui(n22(n),k)-s22(n)*vi(n22(n),k)
                                v22=s22(n)*ui(n22(n),k)+c22(n)*vi(n22(n),k)
                                ub=w11(n)*u11+w21(n)*u21+w12(n)*u12+w22(n)*u22
                                vb=w11(n)*v11+w21(n)*v21+w12(n)*v12+w22(n)*v22
                                uo(n,k)=uo(n,k)+wb*ub
                                vo(n,k)=vo(n,k)+wb*vb
                                wo(n,k)=wo(n,k)+wb
                            else
                                if(li(n11(n),k)) then
                                    u11=c11(n)*ui(n11(n),k)-s11(n)*vi(n11(n),k)
                                    v11=s11(n)*ui(n11(n),k)+c11(n)*vi(n11(n),k)
                                    uo(n,k)=uo(n,k)+wb*w11(n)*u11
                                    vo(n,k)=vo(n,k)+wb*w11(n)*v11
                                    wo(n,k)=wo(n,k)+wb*w11(n)
                                endif
                                if(li(n21(n),k)) then
                                    u21=c21(n)*ui(n21(n),k)-s21(n)*vi(n21(n),k)
                                    v21=s21(n)*ui(n21(n),k)+c21(n)*vi(n21(n),k)
                                    uo(n,k)=uo(n,k)+wb*w21(n)*u21
                                    vo(n,k)=vo(n,k)+wb*w21(n)*v21
                                    wo(n,k)=wo(n,k)+wb*w21(n)
                                endif
                                if(li(n12(n),k)) then
                                    u12=c12(n)*ui(n12(n),k)-s12(n)*vi(n12(n),k)
                                    v12=s12(n)*ui(n12(n),k)+c12(n)*vi(n12(n),k)
                                    uo(n,k)=uo(n,k)+wb*w12(n)*u12
                                    vo(n,k)=vo(n,k)+wb*w12(n)*v12
                                    wo(n,k)=wo(n,k)+wb*w12(n)
                                endif
                                if(li(n22(n),k)) then
                                    u22=c22(n)*ui(n22(n),k)-s22(n)*vi(n22(n),k)
                                    v22=s22(n)*ui(n22(n),k)+c22(n)*vi(n22(n),k)
                                    uo(n,k)=uo(n,k)+wb*w22(n)*u22
                                    vo(n,k)=vo(n,k)+wb*w22(n)*v22
                                    wo(n,k)=wo(n,k)+wb*w22(n)
                                endif
                            endif
                        endif
                    enddo
                    !$omp end parallel do
                enddo
            endif
        enddo
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE OUTPUT BITMAPS AND FIELDS
        ! KM is often 1, do not put OMP PARALLEL here
        do k=1,km
            ibo(k)=ibi(k)
            !$omp parallel do private(n,urot,vrot) schedule(static)
            do n=1,no
                lo(n,k)=wo(n,k).ge.pmp*nb4
                if(lo(n,k)) then
                    uo(n,k)=uo(n,k)/wo(n,k)
                    vo(n,k)=vo(n,k)/wo(n,k)
                    urot=crot(n)*uo(n,k)-srot(n)*vo(n,k)
                    vrot=srot(n)*uo(n,k)+crot(n)*vo(n,k)
                    uo(n,k)=urot
                    vo(n,k)=vrot
                else
                    ibo(k)=1
                    uo(n,k)=0.
                    vo(n,k)=0.
                endif
            enddo
            !$omp end parallel do
        enddo

        select type(grid_out2)
        type is(ip_equid_cylind_grid)
            call polfixv(no,mo,km,rlat,rlon,ibo,lo,uo,vo)
        endselect
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    endsubroutine interpolate_budget_vector

endmodule budget_interp_mod
