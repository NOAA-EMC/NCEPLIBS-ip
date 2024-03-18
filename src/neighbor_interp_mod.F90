!> @file
!> @brief Interpolate scalar and vector fields with neighbor interpolation.
!> @author Mark Iredell @date 96-04-10

!> @brief Interpolate scalar fields (neighbor).
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | Iredell | initial
!> 1999-04-08 | Iredell | split ijkgds into two pieces
!> 2001-06-18 | Iredell | include spiral search option
!> 2006-01-04 | Gayno | minor bug fix
!> 2007-10-30 | Iredell | save weights and thread for performance
!> 2012-06-26 | Gayno | fix out-of-bounds error. see nceplibs ticket #9.
!> 2015-01-27 | Gayno | replace calls to gdswiz() with new merged version of gdswzd().
!> 2015-07-13 | Gayno | replace grib 1 kgds arrays with grib 2 grid definition template arrays.
!>
!> @author Mark Iredell @date 96-04-10
module neighbor_interp_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  implicit none

  private
  public :: interpolate_neighbor

  interface interpolate_neighbor
    module procedure interpolate_neighbor_scalar
    module procedure interpolate_neighbor_vector
  endinterface interpolate_neighbor

  ! Smallest positive real value (use for equality comparisons)
  real :: tinyreal=tiny(1.0)

contains

  !> Interpolate scalar fields (neighbor).
  !>
  !> This subprogram performs neighbor interpolation from any grid to
  !> any grid for scalar fields.
  !>
  !> Options allow choosing the width of the grid square (ipopt(1)) to
  !> search for valid data, which defaults to 1 (if ipopt(1)=-1). Odd
  !> width squares are centered on the nearest input grid point; even
  !> width squares are centered on the nearest four input grid points.
  !> Squares are searched for valid data in a spiral pattern starting
  !> from the center. No searching is done where the output grid is
  !> outside the input grid. Only horizontal interpolation is
  !> performed.
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
  !> latitudes and longitudes are also returned. On the other hand,
  !> the output can be a set of station points if igdtnumo<0, in which
  !> case the number of points and their latitudes and longitudes must
  !> be input.
  !>
  !> Input bitmaps will be interpolated to output bitmaps.
  !>
  !> Output bitmaps will also be created when the output grid extends
  !> outside of the domain of the input grid. The output field is set
  !> to 0 where the output bitmap is off.
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1) is width of
  !> square to examine in spiral search (defaults to 1 if ipopt(1)=-1)
  !> @param[in] grid_in The input grid.
  !> @param[in] grid_out The output grid.
  !> @param[in] mi skip number between input grid fields if km>1 or
  !> dimension of input grid fields if km=1.
  !> @param[in] mo skip number between output grid fields if km>1 or
  !> dimension of output grid fields if km=1
  !> @param[in] km number of fields to interpolate
  !> @param[in] ibi (km) input bitmap flags
  !> @param[in] li (mi,km) input bitmaps (if some ibi(k)=1)
  !> @param[in] gi (mi,km) input fields to interpolate
  !> @param[inout] no number of output points (only if igdtnumo<0)
  !> @param[inout] rlat (no) output latitudes in degrees (if igdtnumo<0)
  !> @param[inout] rlon (no) output longitudes in degrees (if igdtnumo<0)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] go (mo,km) output fields interpolated
  !> @param[out] iret return code
  !> - 0    successful interpolation
  !> - 2    unrecognized input grid or no grid overlap
  !> - 3    unrecognized output grid
  !>
  !> @author Mark Iredell @date 96-04-10
  !> @author Eric Engle @date 23-05-04
  subroutine interpolate_neighbor_scalar(ipopt,grid_in,grid_out, &
                                         mi,mo,km,ibi,li,gi, &
                                         no,rlat,rlon,ibo,lo,go,iret)
    class(ip_grid),intent(in) :: grid_in,grid_out
    integer,intent(in) :: ipopt(20)
    integer,intent(in) :: mi,mo,km
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
    real,parameter     :: fill=-9999.
    !
    integer                              :: i1,j1,ixs,jxs
    integer                              :: mspiral,n,k,nk
    integer                              :: nv
    integer                              :: mx,kxs,kxt,ix,jx,nx
    !
    logical                              :: same_gridi,same_grido
    !
    real                                 :: xpts(mo),ypts(mo)
    logical :: to_station_points

    integer,save :: nox=-1,iretx=-1
    integer,allocatable,save :: nxy(:)
    real,allocatable,save :: rlatx(:),rlonx(:),xptsx(:),yptsx(:)
    class(ip_grid),allocatable,save :: prev_grid_in,prev_grid_out
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    iret=0
    mspiral=max(ipopt(1),1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(.not.allocated(prev_grid_in).or..not.allocated(prev_grid_out)) then
      allocate(prev_grid_in,source=grid_in)
      allocate(prev_grid_out,source=grid_out)

      same_gridi=.false.
      same_grido=.false.
    else
      same_gridi=grid_in.eq.prev_grid_in
      same_grido=grid_out.eq.prev_grid_out

      if(.not.same_gridi.or..not.same_grido) then
        deallocate(prev_grid_in)
        deallocate(prev_grid_out)

        allocate(prev_grid_in,source=grid_in)
        allocate(prev_grid_out,source=grid_out)
      endif
    endif

    select type(grid_out)
    type is(ip_station_points_grid)
      to_station_points=.true.
    class default
      to_station_points=.false.
    endselect
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SAVE OR SKIP WEIGHT COMPUTATION
    if(iret.eq.0.and.(to_station_points.or..not.same_gridi.or..not.same_grido)) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
      call gdswzd(grid_out,0,mo,fill,xpts,ypts,rlon,rlat,no)
      if(no.eq.0) iret=3
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  LOCATE INPUT POINTS
      call gdswzd(grid_in,-1,no,fill,xpts,ypts,rlon,rlat,nv)
      if(iret.eq.0.and.nv.eq.0) iret=2
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  ALLOCATE AND SAVE GRID DATA
      if(nox.ne.no) then
        if(nox.ge.0) deallocate(rlatx,rlonx,xptsx,yptsx,nxy)
        allocate(rlatx(no),rlonx(no),xptsx(no),yptsx(no),nxy(no))
        nox=no
      endif
      iretx=iret
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  COMPUTE WEIGHTS
      if(iret.eq.0) then
        !$omp parallel do private(n) schedule(static)
        do n=1,no
          rlonx(n)=rlon(n)
          rlatx(n)=rlat(n)
          xptsx(n)=xpts(n)
          yptsx(n)=ypts(n)
          if(abs(xpts(n)-fill).gt.tinyreal.and.abs(ypts(n)-fill).gt.tinyreal) then
            nxy(n)=grid_in%field_pos(nint(xpts(n)),nint(ypts(n)))
          else
            nxy(n)=0
          endif
        enddo
      endif
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE OVER ALL FIELDS
    if(iret.eq.0.and.iretx.eq.0) then
      if(.not.to_station_points) then
        no=nox
        do n=1,no
          rlon(n)=rlonx(n)
          rlat(n)=rlatx(n)
        enddo
      endif
      do n=1,no
        xpts(n)=xptsx(n)
        ypts(n)=yptsx(n)
      enddo
      !$omp parallel do private(nk,k,n,i1,j1,ixs,jxs,mx,kxs,kxt,ix,jx,nx) schedule(static)
      do nk=1,no*km
        k=(nk-1)/no+1
        n=nk-no*(k-1)
        go(n,k)=0
        lo(n,k)=.false.
        if(nxy(n).gt.0) then
          if(ibi(k).eq.0.or.li(nxy(n),k)) then
            go(n,k)=gi(nxy(n),k)
            lo(n,k)=.true.
            ! SPIRAL AROUND UNTIL VALID DATA IS FOUND.
          elseif(mspiral.gt.1) then
            i1=nint(xpts(n))
            j1=nint(ypts(n))
            ixs=int(sign(1.,xpts(n)-i1))
            jxs=int(sign(1.,ypts(n)-j1))
            do mx=2,mspiral**2
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
              if(nx.gt.0) then
                if(li(nx,k)) then
                  go(n,k)=gi(nx,k)
                  lo(n,k)=.true.
                  exit
                endif
              endif
            enddo
          endif
        endif
      enddo

      do k=1,km
        ibo(k)=ibi(k)
        if(.not.all(lo(1:no,k))) ibo(k)=1
      enddo

      select type(grid_out)
      type is(ip_equid_cylind_grid)
        call polfixs(no,mo,km,rlat,ibo,lo,go)
      endselect
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else
      if(iret.eq.0) iret=iretx
      if(.not.to_station_points) no=0
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endsubroutine interpolate_neighbor_scalar

  !> Interpolate vector fields (neighbor).
  !>
  !> This subprogram performs neighbor interpolation
  !> from any grid to any grid for vector fields.
  !>
  !> Options allow choosing the width of the grid square (ipopt(1)) to
  !> search for valid data, which defaults to 1 (IF IPOPT(1)=-1). Odd
  !> width squares are centered on the nearest input grid point; even
  !> width squares are centered on the nearest four input grid points.
  !>
  !> Squares are searched for valid data in a spiral pattern starting
  !> from the center. no searching is done where the output grid is
  !> outside the input grid. Only horizontal interpolation is
  !> performed.
  !>
  !> The input and output grids are defined by their grib 2 grid
  !> definition template as decoded by the ncep g2 library.
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
  !> The input and output vectors are rotated so that they are either
  !> resolved relative to the defined grid in the direction of
  !> increasing x and y coordinates or resolved relative to easterly
  !> and northerly directions, as designated by their respective grid
  !> definition sections.
  !>
  !> As an added bonus the number of output grid points and their
  !> latitudes and longitudes are also returned along with their
  !> vector rotation parameters. On the other hand, the output can be
  !> a set of station points if igdtnumo<0, in which case the number
  !> of points and their latitudes and longitudes must be input along
  !> with their vector rotation parameters.
  !>
  !> Input bitmaps will be interpolated to output bitmaps. output
  !> bitmaps will also be created when the output grid extends outside
  !> of the domain of the input grid. The output field is set to 0
  !> where the output bitmap is off.
  !>
  !> @param[in] ipopt (20) interpolation options ipopt(1) is width of
  !> square to examine in spiral search (defaults to 1 if ipopt(1)=-1)
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
  !> @param[inout] no number of output points (only if igdtnumo>=0)
  !> @param[inout] rlat (mo) output latitudes in degrees (if igdtnumo<0)
  !> @param[inout] rlon (mo) output longitudes in degrees (if igdtnumo<0)
  !> @param[inout] crot (mo) vector rotation cosines (if igdtnumo<0)
  !> @param[inout] srot (mo) vector rotation sines (if igdtnumo<0)
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] ibo (km) output bitmap flags
  !> @param[out] lo (mo,km) output bitmaps (always output)
  !> @param[out] uo (mo,km) output u-component fields interpolated
  !> @param[out] vo (mo,km) output v-component fields interpolated
  !> @param[out] iret return code
  !> - 0    successful interpolation
  !> - 2    unrecognized input grid or no grid overlap
  !> - 3    unrecognized output grid
  !>
  !> @author Mark Iredell @date 96-04-10
  !> @author Eric Engle @date 23-05-04
  subroutine interpolate_neighbor_vector(ipopt,grid_in,grid_out, &
                                         mi,mo,km,ibi,li,ui,vi, &
                                         no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)

    class(ip_grid),intent(in) :: grid_in,grid_out
    integer,intent(in) :: ipopt(20)
    integer,intent(in) :: ibi(km),mi,mo,km
    integer,intent(inout) :: no
    integer,intent(out) :: iret,ibo(km)
    !
    logical*1,intent(in) :: li(mi,km)
    logical*1,intent(out) :: lo(mo,km)
    !
    real,intent(in) :: ui(mi,km),vi(mi,km)
    real,intent(inout) :: crot(mo),srot(mo)
    real,intent(inout) :: rlat(mo),rlon(mo)
    real,intent(out) :: uo(mo,km),vo(mo,km)
    !
    real,parameter     :: fill=-9999.
    !
    integer                           :: i1,j1,ixs,jxs,mx
    integer                           :: kxs,kxt,ix,jx,nx
    integer                           :: mspiral,n,k,nk,nv
    !
    logical                           :: same_gridi,same_grido
    !
    real                              :: cx,sx,cm,sm,urot,vrot
    real                              :: xpts(mo),ypts(mo)
    real                              :: croi(mi),sroi(mi)
    real                              :: xpti(mi),ypti(mi),rloi(mi),rlai(mi)

    logical :: to_station_points

    integer,save :: nox=-1,iretx=-1
    integer,allocatable,save :: nxy(:)

    real,allocatable,save :: rlatx(:),rlonx(:),xptsx(:),yptsx(:)
    real,allocatable,save :: crotx(:),srotx(:),cxy(:),sxy(:)
    class(ip_grid),allocatable,save :: prev_grid_in,prev_grid_out
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    iret=0
    mspiral=max(ipopt(1),1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(.not.allocated(prev_grid_in).or..not.allocated(prev_grid_out)) then
      allocate(prev_grid_in,source=grid_in)
      allocate(prev_grid_out,source=grid_out)

      same_gridi=.false.
      same_grido=.false.
    else
      same_gridi=grid_in.eq.prev_grid_in
      same_grido=grid_out.eq.prev_grid_out

      if(.not.same_gridi.or..not.same_grido) then
        deallocate(prev_grid_in)
        deallocate(prev_grid_out)

        allocate(prev_grid_in,source=grid_in)
        allocate(prev_grid_out,source=grid_out)
      endif
    endif

    select type(grid_out)
    type is(ip_station_points_grid)
      to_station_points=.true.
    class default
      to_station_points=.false.
    endselect

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SAVE OR SKIP WEIGHT COMPUTATION
    if(iret.eq.0.and.(to_station_points.or..not.same_gridi.or..not.same_grido)) then
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
      call gdswzd(grid_out,0,mo,fill,xpts,ypts,rlon,rlat,no,crot,srot)
      if(no.eq.0) iret=3
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  LOCATE INPUT POINTS
      call gdswzd(grid_in,-1,no,fill,xpts,ypts,rlon,rlat,nv)
      if(iret.eq.0.and.nv.eq.0) iret=2
      call gdswzd(grid_in,0,mi,fill,xpti,ypti,rloi,rlai,nv,croi,sroi)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  ALLOCATE AND SAVE GRID DATA
      if(nox.ne.no) then
        if(nox.ge.0) deallocate(rlatx,rlonx,xptsx,yptsx,crotx,srotx,nxy,cxy,sxy)
        allocate(rlatx(no),rlonx(no),xptsx(no),yptsx(no), &
                 crotx(no),srotx(no),nxy(no),cxy(no),sxy(no))
        nox=no
      endif
      iretx=iret
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  COMPUTE WEIGHTS
      if(iret.eq.0) then
        !$omp parallel do private(n,cm,sm) schedule(static)
        do n=1,no
          rlonx(n)=rlon(n)
          rlatx(n)=rlat(n)
          xptsx(n)=xpts(n)
          yptsx(n)=ypts(n)
          crotx(n)=crot(n)
          srotx(n)=srot(n)
          if(abs(xpts(n)-fill).gt.tinyreal.and.abs(ypts(n)-fill).gt.tinyreal) then
            nxy(n)=grid_in%field_pos(nint(xpts(n)),nint(ypts(n)))
            if(nxy(n).gt.0) then
              call movect(rlai(nxy(n)),rloi(nxy(n)),rlat(n),rlon(n),cm,sm)
              cxy(n)=cm*croi(nxy(n))+sm*sroi(nxy(n))
              sxy(n)=sm*croi(nxy(n))-cm*sroi(nxy(n))
            endif
          else
            nxy(n)=0
          endif
        enddo
      endif
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  INTERPOLATE OVER ALL FIELDS
    if(iret.eq.0.and.iretx.eq.0) then
      if(.not.to_station_points) then
        no=nox
        do n=1,no
          rlon(n)=rlonx(n)
          rlat(n)=rlatx(n)
          crot(n)=crotx(n)
          srot(n)=srotx(n)
        enddo
      endif
      do n=1,no
        xpts(n)=xptsx(n)
        ypts(n)=yptsx(n)
      enddo
      !$omp parallel do &
        !$omp private(nk,k,n,i1,j1,ixs,jxs,mx,kxs,kxt,ix,jx,nx) &
        !$omp private(cm,sm,cx,sx,urot,vrot) schedule(static)
      do nk=1,no*km
        k=(nk-1)/no+1
        n=nk-no*(k-1)
        uo(n,k)=0
        vo(n,k)=0
        lo(n,k)=.false.
        if(nxy(n).gt.0) then
          if(ibi(k).eq.0.or.li(nxy(n),k)) then
            urot=cxy(n)*ui(nxy(n),k)-sxy(n)*vi(nxy(n),k)
            vrot=sxy(n)*ui(nxy(n),k)+cxy(n)*vi(nxy(n),k)
            uo(n,k)=crot(n)*urot-srot(n)*vrot
            vo(n,k)=srot(n)*urot+crot(n)*vrot
            lo(n,k)=.true.
            ! SPIRAL AROUND UNTIL VALID DATA IS FOUND.
          elseif(mspiral.gt.1) then
            i1=nint(xpts(n))
            j1=nint(ypts(n))
            ixs=int(sign(1.,xpts(n)-i1))
            jxs=int(sign(1.,ypts(n)-j1))
            do mx=2,mspiral**2
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
              if(nx.gt.0) then
                if(li(nx,k)) then
                  call movect(rlai(nx),rloi(nx),rlat(n),rlon(n),cm,sm)
                  cx=cm*croi(nx)+sm*sroi(nx)
                  sx=sm*croi(nx)-cm*sroi(nx)
                  urot=cx*ui(nx,k)-sx*vi(nx,k)
                  vrot=sx*ui(nx,k)+cx*vi(nx,k)
                  uo(n,k)=crot(n)*urot-srot(n)*vrot
                  vo(n,k)=srot(n)*urot+crot(n)*vrot
                  lo(n,k)=.true.
                  exit
                endif
              endif
            enddo
          endif
        endif
      enddo
      do k=1,km
        ibo(k)=ibi(k)
        if(.not.all(lo(1:no,k))) ibo(k)=1
      enddo

      select type(grid_out)
      type is(ip_equid_cylind_grid)
        call polfixv(no,mo,km,rlat,rlon,ibo,lo,uo,vo)
      endselect

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    else
      if(iret.eq.0) iret=iretx
      if(.not.to_station_points) no=0
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endsubroutine interpolate_neighbor_vector

endmodule neighbor_interp_mod
