!> @file
!> @brief Make multiple pole scalar values consistent.
!> @author Kyle Gerheiser @date 2021-07-21

!> @brief Make multiple pole scalar values consistent.
!> @author Iredell, Kyle Gerheiser
module polfix_mod
  implicit none

  private
  public :: polfixs,polfixv

contains

  !> Make multiple pole scalar values consistent.
  !>
  !> This subprogram averages multiple pole scalar values
  !> on a latitude/longitude grid. Bitmaps may be averaged too.
  !>
  !> @param[in] nm integer number of grid points
  !> @param[in] nx integer leading dimension of fields
  !> @param[in] km integer number of fields
  !> @param[in] rlat real (no) latitudes in degrees
  !> @param[in] ib integer (km) bitmap flags
  !> @param[out] lo logical*1 (nx,km) bitmaps (if some ib(k)=1)
  !> @param[out] go real (nx,km) fields
  !>
  !> @author Iredell @date 96-04-10
  subroutine polfixs(nm,nx,km,rlat,ib,lo,go)
    implicit none
    !
    integer,intent(in) :: nm,nx,km
    integer,intent(in) :: ib(km)
    !
    logical*1,intent(inout) :: lo(nx,km)
    !
    real,intent(in) :: rlat(nm)
    real,intent(inout) :: go(nx,km)
    !
    real,parameter     :: rlatnp=89.9995
    real,parameter     :: rlatsp=-rlatnp
    !
    integer                   :: k,n
    !
    real                      :: wnp,gnp,tnp,wsp,gsp,tsp
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,km
      wnp=0.
      gnp=0.
      tnp=0.
      wsp=0.
      gsp=0.
      tsp=0.
      !  AVERAGE MULTIPLE POLE VALUES
      !$omp parallel do private(n) reduction(+:wnp,gnp,tnp,wsp,gsp,tsp) schedule(static)
      do n=1,nm
        if(rlat(n).ge.rlatnp) then
          wnp=wnp+1
          if(ib(k).eq.0.or.lo(n,k)) then
            gnp=gnp+go(n,k)
            tnp=tnp+1
          endif
        elseif(rlat(n).le.rlatsp) then
          wsp=wsp+1
          if(ib(k).eq.0.or.lo(n,k)) then
            gsp=gsp+go(n,k)
            tsp=tsp+1
          endif
        endif
      enddo
      !$omp end parallel do
      !  DISTRIBUTE AVERAGE VALUES BACK TO MULTIPLE POLES
      if(wnp.gt.1) then
        if(tnp.ge.wnp/2) then
          gnp=gnp/tnp
        else
          gnp=0.
        endif
        !$omp parallel do private(n) schedule(static)
        do n=1,nm
          if(rlat(n).ge.rlatnp) then
            if(ib(k).ne.0) lo(n,k)=tnp.ge.wnp/2
            go(n,k)=gnp
          endif
        enddo
        !$omp end parallel do
      endif
      if(wsp.gt.1) then
        if(tsp.ge.wsp/2) then
          gsp=gsp/tsp
        else
          gsp=0.
        endif
        !$omp parallel do private(n) schedule(static)
        do n=1,nm
          if(rlat(n).le.rlatsp) then
            if(ib(k).ne.0) lo(n,k)=tsp.ge.wsp/2
            go(n,k)=gsp
          endif
        enddo
        !$omp end parallel do
      endif
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endsubroutine polfixs

  !> Make multiple pole vector values consistent,
  !>
  !> This subprogram averages multiple pole vector values on a
  !> latitude/longitude grid. Bitmaps may be averaged too. Vectors
  !> are rotated with respect to their longitude.
  !>
  !> @param[in] nm integer number of grid points
  !> @param[in] nx integer leading dimension of fields
  !> @param[in] km integer number of fields
  !> @param[in] rlat real (nm) latitudes in degrees
  !> @param[in] rlon real (nm) longitudes in degrees
  !> @param[in] ib integer (km) bitmap flags
  !> @param[inout] lo logical*1 (nx,km) bitmaps (if some ib(k)=1)
  !> @param[inout] uo real (nx,km) u-winds
  !> @param[inout] vo real (nx,km) v-winds
  !>
  !> @author Iredell @date 96-04-10
  subroutine polfixv(nm,nx,km,rlat,rlon,ib,lo,uo,vo)
    implicit none
    !
    integer,intent(in) :: ib(km),nm,nx,km
    !
    logical*1,intent(inout) :: lo(nx,km)
    !
    real,intent(in) :: rlat(nm),rlon(nm)
    real,intent(inout) :: uo(nx,km),vo(nx,km)
    !
    real,parameter     :: rlatnp=89.9995
    real,parameter     :: rlatsp=-rlatnp
    real,parameter     :: pi=3.14159265358979
    real,parameter     :: dpr=180./pi
    !
    integer                     :: k,n
    !
    real                        :: clon(nm),slon(nm)
    real                        :: tnp,unp,vnp,wnp
    real                        :: tsp,usp,vsp,wsp
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !$omp parallel do private(n) schedule(static)
    do n=1,nm
      clon(n)=cos(rlon(n)/dpr)
      slon(n)=sin(rlon(n)/dpr)
    enddo
    !$omp end parallel do
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,km
      wnp=0.
      unp=0.
      vnp=0.
      tnp=0.
      wsp=0.
      usp=0.
      vsp=0.
      tsp=0.
      !  AVERAGE MULTIPLE POLE VALUES
      !$omp parallel do private(n) reduction(+:wnp,unp,vnp,tnp,wsp,usp,vsp,tsp) schedule(static)
      do n=1,nm
        if(rlat(n).ge.rlatnp) then
          wnp=wnp+1
          if(ib(k).eq.0.or.lo(n,k)) then
            unp=unp+(clon(n)*uo(n,k)-slon(n)*vo(n,k))
            vnp=vnp+(slon(n)*uo(n,k)+clon(n)*vo(n,k))
            tnp=tnp+1
          endif
        elseif(rlat(n).le.rlatsp) then
          wsp=wsp+1
          if(ib(k).eq.0.or.lo(n,k)) then
            usp=usp+(clon(n)*uo(n,k)+slon(n)*vo(n,k))
            vsp=vsp+(-slon(n)*uo(n,k)+clon(n)*vo(n,k))
            tsp=tsp+1
          endif
        endif
      enddo
      !$omp end parallel do
      !  DISTRIBUTE AVERAGE VALUES BACK TO MULTIPLE POLES
      if(wnp.gt.1) then
        if(tnp.ge.wnp/2) then
          unp=unp/tnp
          vnp=vnp/tnp
        else
          unp=0.
          vnp=0.
        endif
        !$omp parallel do private(n) schedule(static)
        do n=1,nm
          if(rlat(n).ge.rlatnp) then
            if(ib(k).ne.0) lo(n,k)=tnp.ge.wnp/2
            uo(n,k)=clon(n)*unp+slon(n)*vnp
            vo(n,k)=-slon(n)*unp+clon(n)*vnp
          endif
        enddo
        !$omp end parallel do
      endif
      if(wsp.gt.1) then
        if(tsp.ge.wsp/2) then
          usp=usp/wsp
          vsp=vsp/wsp
        else
          usp=0.
          vsp=0.
        endif
        !$omp parallel do private(n) schedule(static)
        do n=1,nm
          if(rlat(n).le.rlatsp) then
            if(ib(k).ne.0) lo(n,k)=tsp.ge.wsp/2
            uo(n,k)=clon(n)*usp-slon(n)*vsp
            vo(n,k)=slon(n)*usp+clon(n)*vsp
          endif
        enddo
        !$omp end parallel do
      endif
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endsubroutine polfixv
endmodule polfix_mod
