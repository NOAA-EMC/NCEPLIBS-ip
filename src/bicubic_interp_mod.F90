!> @file
!> @brief Bicubic interpolation routines for scalars and vectors
!> @author Mark Iredell, Kyle Gerheiser, Eric Engle

!> @brief Bicubic interpolation routines for scalars and vectors.
!>
!> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
module bicubic_interp_mod
  use gdswzd_mod
  use polfix_mod
  use ip_grids_mod
  implicit none

  private
  public :: interpolate_bicubic

  interface interpolate_bicubic
    module procedure interpolate_bicubic_scalar
    module procedure interpolate_bicubic_vector
  endinterface interpolate_bicubic

  ! Smallest positive real value (use for equality comparisons)
  real :: tinyreal=tiny(1.0)

contains

  !> This subprogram performs bicubic interpolation
  !> from any grid to any grid for scalar fields.
  !>
  !> @details Bitmaps are now allowed even when invalid points are within
  !> the bicubic template provided the minimum weight is reached.
  !>
  !> Options allow choices between straight bicubic (ipopt(1)=0)
  !> and constrained bicubic (ipopt(1)=1) where the value is
  !> confined within the range of the surrounding 16 points.
  !>
  !> Another option is the minimum percentage for mask,
  !> i.e. percent valid input data required to make output data,
  !> (ipopt(2)) which defaults to 50 (if ipopt(2)=-1).
  !>
  !> Bilinear used within one grid length of boundaries.
  !> Only horizontal interpolation is performed.
  !>
  !> The code recognizes the following projections, where
  !> for the input and output grids, respectively:
  !> as an added bonus the number of output grid points
  !> and their latitudes and longitudes are also returned.
  !> On the other hand, the output can be a set of station points
  !> if igdtnumo<0, in which case the number of points
  !> and their latitudes and longitudes must be input.
  !> output bitmaps will only be created when the output grid
  !> extends outside of the domain of the input grid.
  !> the output field is set to 0 where the output bitmap is off.
  !>
  !> @param[in] ipopt Interpolation options.
  !> - ipopt(1)=0 For straight bicubic;
  !> - ipopt(1)=1 For constrained bicubic where value is confined within the range of the surrounding 4 points.
  !> - ipopt(2) Is minimum percentage for mask (defaults to 50 if ipopt(2)=-1)
  !>
  !> @param[in] grid_in Input grid.
  !> @param[in] grid_out Output grid.
  !> @param[in]  mi Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !> @param[out] mo Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !> @param[in]  km Number of fields to interpolate.
  !> @param[in]  ibi Input bitmap flags.
  !> @param[in]  li Input bitmaps (if some ibi(k)=1).
  !> @param[in]  gi Input fields to interpolate.
  !> @param[in,out] no  Number of output points (only if igdtnumo<0).
  !> @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0).
  !> @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0).
  !> @param[out] ibo Output bitmap flags.
  !> @param[out] lo Output bitmaps (always output).
  !> @param[out] go Output fields interpolated.
  !> @param[out] iret Return code.
  !> - 0 successful interpolation,
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !>
  !> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
  subroutine interpolate_bicubic_scalar(ipopt,grid_in,grid_out, &
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
    integer                               :: ijx(4),ijy(4)
    integer                               :: mcon,mp,n,i,j,k
    integer                               :: nk,nv
    logical                               :: same_gridi,same_grido
    !
    real                                  :: pmp,xij,yij,xf,yf
    real                                  :: g,w,gmin,gmax
    real                                  :: wx(4),wy(4)
    real                                  :: xpts(mo),ypts(mo)
    logical :: to_station_points

    ! Save coeffecients between calls and only compute if grids have changed
    real,allocatable,save  :: rlatx(:),rlonx(:)
    real,allocatable,save  :: wxy(:,:,:)
    integer,save  :: nox=-1,iretx=-1
    integer,allocatable,save  :: nxy(:,:,:),nc(:)
    class(ip_grid),allocatable,save :: prev_grid_in,prev_grid_out

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    iret=0
    mcon=ipopt(1)
    mp=ipopt(2)
    if(mp.eq.-1.or.mp.eq.0) mp=50
    if(mp.lt.0.or.mp.gt.100) iret=32
    pmp=mp*0.01

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
        if(nox.ge.0) deallocate(rlatx,rlonx,nc,nxy,wxy)
        allocate(rlatx(no),rlonx(no),nc(no),nxy(4,4,no),wxy(4,4,no))
        nox=no
      endif
      iretx=iret
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  COMPUTE WEIGHTS
      if(iret.eq.0) then
        !$omp parallel do private(n,xij,yij,ijx,ijy,xf,yf,j,i,wx,wy) schedule(static)
        do n=1,no
          rlonx(n)=rlon(n)
          rlatx(n)=rlat(n)
          xij=xpts(n)
          yij=ypts(n)
          if(abs(xij-fill).gt.tinyreal.and.abs(yij-fill).gt.tinyreal) then
            ijx(1:4)=floor(xij-1)+(/0,1,2,3/)
            ijy(1:4)=floor(yij-1)+(/0,1,2,3/)
            xf=xij-ijx(2)
            yf=yij-ijy(2)
            do j=1,4
              do i=1,4
                nxy(i,j,n)=grid_in%field_pos(ijx(i),ijy(j))
              enddo
            enddo
            if(minval(nxy(1:4,1:4,n)).gt.0) then
              !  BICUBIC WHERE 16-POINT STENCIL IS AVAILABLE
              nc(n)=1
              wx(1)=xf*(1-xf)*(2-xf)/(-6.)
              wx(2)=(xf+1)*(1-xf)*(2-xf)/2.
              wx(3)=(xf+1)*xf*(2-xf)/2.
              wx(4)=(xf+1)*xf*(1-xf)/(-6.)
              wy(1)=yf*(1-yf)*(2-yf)/(-6.)
              wy(2)=(yf+1)*(1-yf)*(2-yf)/2.
              wy(3)=(yf+1)*yf*(2-yf)/2.
              wy(4)=(yf+1)*yf*(1-yf)/(-6.)
            else
              !  BILINEAR ELSEWHERE NEAR THE EDGE OF THE GRID
              nc(n)=2
              wx(1)=0
              wx(2)=(1-xf)
              wx(3)=xf
              wx(4)=0
              wy(1)=0
              wy(2)=(1-yf)
              wy(3)=yf
              wy(4)=0
            endif
            do j=1,4
              do i=1,4
                wxy(i,j,n)=wx(i)*wy(j)
              enddo
            enddo
          else
            nc(n)=0
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
      !$omp parallel do private(nk,k,n,g,w,gmin,gmax,j,i) schedule(static)
      do nk=1,no*km
        k=(nk-1)/no+1
        n=nk-no*(k-1)
        if(nc(n).gt.0) then
          g=0
          w=0
          if(mcon.gt.0) gmin=huge(gmin)
          if(mcon.gt.0) gmax=-huge(gmax)
          do j=nc(n),5-nc(n)
            do i=nc(n),5-nc(n)
              if(nxy(i,j,n).gt.0) then
                if(ibi(k).eq.0.or.li(nxy(i,j,n),k)) then
                  g=g+wxy(i,j,n)*gi(nxy(i,j,n),k)
                  w=w+wxy(i,j,n)
                  if(mcon.gt.0) gmin=min(gmin,gi(nxy(i,j,n),k))
                  if(mcon.gt.0) gmax=max(gmax,gi(nxy(i,j,n),k))
                endif
              endif
            enddo
          enddo
          lo(n,k)=w.ge.pmp
          if(lo(n,k)) then
            go(n,k)=g/w
            if(mcon.gt.0) go(n,k)=min(max(go(n,k),gmin),gmax)
          else
            go(n,k)=0.
          endif
        else
          lo(n,k)=.false.
          go(n,k)=0.
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
    else
      if(iret.eq.0) iret=iretx
      if(.not.to_station_points) no=0
    endif
  endsubroutine interpolate_bicubic_scalar

  !> This subprogram performs bicubic interpolation from any grid to
  !> any grid for vector fields.
  !>
  !> Bitmaps are now allowed even when invalid points are within the
  !> bicubic template provided the minimum weight is reached.
  !>
  !> Options allow choices between straight bicubic (ipopt(1)=0)
  !> and constrained bicubic (ipopt(1)=1) where the value is
  !> confined within the range of the surrounding 16 points.
  !>
  !> Another option is the minimum percentage for mask,
  !> i.e. percent valid input data required to make output data,
  !> (ipopt(2)) which defaults to 50 (if ipopt(2)=-1).
  !>
  !> Bilinear used within one grid length of boundaries.
  !> Only horizontal interpolation is performed.
  !>
  !> The code recognizes the following projections, where
  !> for the input and output grids, respectively:
  !> as an added bonus the number of output grid points
  !> and their latitudes and longitudes are also returned.
  !> On the other hand, the output can be a set of station points
  !> if igdtnumo<0, in which case the number of points
  !> and their latitudes and longitudes must be input.
  !>
  !> Output bitmaps will only be created when the output grid
  !> extends outside of the domain of the input grid.
  !> the output field is set to 0 where the output bitmap is off.
  !>
  !> @param[in] ipopt integer (20) interpolation options
  !> - ipopt(1)=0 for straight bicubic;
  !> - ipopt(1)=1 for constrained bicubic where value is confined within the range of the surrounding 4 points.
  !> - ipopt(2) is minimum percentage for mask (defaults to 50 if ipopt(2)=-1)
  !> @param[in] grid_in Input grid.
  !> @param[in] grid_out Output grid.
  !> @param[in] mi Skip number between input grid fields if km>1 or dimension of input grid fields if km=1.
  !> @param[out] mo Skip number between output grid fields if km>1 or dimension of output grid fields if km=1.
  !> @param[in] km Number of fields to interpolate.
  !> @param[in] ibi Input bitmap flags.
  !> @param[in] li Input bitmaps (if some ibi(k)=1).
  !> @param[in] ui Input u-component fields to interpolate.
  !> @param[in] vi Input v-component fields to interpolate.
  !> @param[in,out] no Number of output points (only if igdtnumo<0).
  !> @param[in,out] rlat Output latitudes in degrees (if igdtnumo<0).
  !> @param[in,out] rlon Output longitudes in degrees (if igdtnumo<0).
  !> @param[in,out] crot Vector rotation cosines (if igdtnumo<0) ugrid=crot*uearth-srot*vearth.
  !> @param[in,out] srot Vector rotation sines (if igdtnumo<0) vgrid=srot*uearth+crot*vearth).
  !> @param[out] ibo Output bitmap flags.
  !> @param[out] lo Output bitmaps (always output).
  !> @param[out] uo Output u-component fields interpolated.
  !> @param[out] vo Output v-component fields interpolated.
  !> @param[out] iret Return code.
  !> - 0 successful interpolation
  !> - 2 unrecognized input grid or no grid overlap
  !> - 3 unrecognized output grid
  !>
  !> @author George Gayno, Mark Iredell, Kyle Gerheiser, Eric Engle
  subroutine interpolate_bicubic_vector(ipopt,grid_in,grid_out, &
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
    real,intent(inout) :: rlat(mo),rlon(mo),crot(mo),srot(mo)
    real,intent(out) :: uo(mo,km),vo(mo,km)
    !
    real,parameter     :: fill=-9999.
    !
    integer                           :: ijx(4),ijy(4)
    integer                           :: mcon,mp,n,i,j,k,nk,nv
    !
    logical                           :: same_gridi,same_grido
    !
    real                              :: cm,sm,urot,vrot
    real                              :: pmp,xij,yij,xf,yf
    real                              :: u,v,w,umin,umax,vmin,vmax
    real                              :: xpts(mo),ypts(mo)
    real                              :: wx(4),wy(4)
    real                              :: xpti(mi),ypti(mi),rloi(mi),rlai(mi)
    real                              :: croi(mi),sroi(mi)

    logical :: to_station_points

    ! Save coeffecients between calls and only compute if grids have changed
    real,allocatable,save  :: rlatx(:),rlonx(:),crotx(:),srotx(:)
    real,allocatable,save  :: wxy(:,:,:),cxy(:,:,:),sxy(:,:,:)
    integer,save  :: nox=-1,iretx=-1
    integer,allocatable,save  :: nxy(:,:,:),nc(:)
    class(ip_grid),allocatable,save :: prev_grid_in,prev_grid_out
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  SET PARAMETERS
    iret=0
    mcon=ipopt(1)
    mp=ipopt(2)
    if(mp.eq.-1.or.mp.eq.0) mp=50
    if(mp.lt.0.or.mp.gt.100) iret=32
    pmp=mp*0.01

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
        if(nox.ge.0) deallocate(rlatx,rlonx,crotx,srotx,nc,nxy,wxy,cxy,sxy)
        allocate(rlatx(no),rlonx(no),crotx(no),srotx(no),nc(no), &
                 nxy(4,4,no),wxy(4,4,no),cxy(4,4,no),sxy(4,4,no))
        nox=no
      endif
      iretx=iret
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  COMPUTE WEIGHTS
      if(iret.eq.0) then
        !$omp parallel do private(n,xij,yij,ijx,ijy,xf,yf,j,i,wx,wy,cm,sm) schedule(static)
        do n=1,no
          rlonx(n)=rlon(n)
          rlatx(n)=rlat(n)
          crotx(n)=crot(n)
          srotx(n)=srot(n)
          xij=xpts(n)
          yij=ypts(n)
          if(abs(xij-fill).gt.tinyreal.and.abs(yij-fill).gt.tinyreal) then
            ijx(1:4)=floor(xij-1)+(/0,1,2,3/)
            ijy(1:4)=floor(yij-1)+(/0,1,2,3/)
            xf=xij-ijx(2)
            yf=yij-ijy(2)
            do j=1,4
              do i=1,4
                nxy(i,j,n)=grid_in%field_pos(ijx(i),ijy(j))
              enddo
            enddo
            if(minval(nxy(1:4,1:4,n)).gt.0) then
              !  BICUBIC WHERE 16-POINT STENCIL IS AVAILABLE
              nc(n)=1
              wx(1)=xf*(1-xf)*(2-xf)/(-6.)
              wx(2)=(xf+1)*(1-xf)*(2-xf)/2.
              wx(3)=(xf+1)*xf*(2-xf)/2.
              wx(4)=(xf+1)*xf*(1-xf)/(-6.)
              wy(1)=yf*(1-yf)*(2-yf)/(-6.)
              wy(2)=(yf+1)*(1-yf)*(2-yf)/2.
              wy(3)=(yf+1)*yf*(2-yf)/2.
              wy(4)=(yf+1)*yf*(1-yf)/(-6.)
            else
              !  BILINEAR ELSEWHERE NEAR THE EDGE OF THE GRID
              nc(n)=2
              wx(1)=0
              wx(2)=(1-xf)
              wx(3)=xf
              wx(4)=0
              wy(1)=0
              wy(2)=(1-yf)
              wy(3)=yf
              wy(4)=0
            endif
            do j=1,4
              do i=1,4
                wxy(i,j,n)=wx(i)*wy(j)
                if(nxy(i,j,n).gt.0) then
                  call movect(rlai(nxy(i,j,n)),rloi(nxy(i,j,n)), &
                              rlat(n),rlon(n),cm,sm)
                  cxy(i,j,n)=cm*croi(nxy(i,j,n))+sm*sroi(nxy(i,j,n))
                  sxy(i,j,n)=sm*croi(nxy(i,j,n))-cm*sroi(nxy(i,j,n))
                endif
              enddo
            enddo
          else
            nc(n)=0
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
      !$omp parallel do private(nk,k,n,u,v,w,umin,umax,vmin,vmax,urot,vrot,j,i) schedule(static)
      do nk=1,no*km
        k=(nk-1)/no+1
        n=nk-no*(k-1)
        if(nc(n).gt.0) then
          u=0
          v=0
          w=0
          if(mcon.gt.0) umin=huge(umin)
          if(mcon.gt.0) umax=-huge(umax)
          if(mcon.gt.0) vmin=huge(vmin)
          if(mcon.gt.0) vmax=-huge(vmax)
          do j=nc(n),5-nc(n)
            do i=nc(n),5-nc(n)
              if(nxy(i,j,n).gt.0) then
                if(ibi(k).eq.0.or.li(nxy(i,j,n),k)) then
                  urot=cxy(i,j,n)*ui(nxy(i,j,n),k)-sxy(i,j,n)*vi(nxy(i,j,n),k)
                  vrot=sxy(i,j,n)*ui(nxy(i,j,n),k)+cxy(i,j,n)*vi(nxy(i,j,n),k)
                  u=u+wxy(i,j,n)*urot
                  v=v+wxy(i,j,n)*vrot
                  w=w+wxy(i,j,n)
                  if(mcon.gt.0) umin=min(umin,urot)
                  if(mcon.gt.0) umax=max(umax,urot)
                  if(mcon.gt.0) vmin=min(vmin,vrot)
                  if(mcon.gt.0) vmax=max(vmax,vrot)
                endif
              endif
            enddo
          enddo
          lo(n,k)=w.ge.pmp
          if(lo(n,k)) then
            urot=crot(n)*u-srot(n)*v
            vrot=srot(n)*u+crot(n)*v
            uo(n,k)=urot/w
            vo(n,k)=vrot/w
            if(mcon.gt.0) uo(n,k)=min(max(uo(n,k),umin),umax)
            if(mcon.gt.0) vo(n,k)=min(max(vo(n,k),vmin),vmax)
          else
            uo(n,k)=0.
            vo(n,k)=0.
          endif
        else
          lo(n,k)=.false.
          uo(n,k)=0.
          vo(n,k)=0.
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
  endsubroutine interpolate_bicubic_vector

endmodule bicubic_interp_mod
