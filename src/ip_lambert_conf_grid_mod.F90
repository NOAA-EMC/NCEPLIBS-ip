!> @file
!> @brief GDS wizard for lambert conformal conical.
!>
!> @author Iredell @date 96-04-10

!> @brief Lambert conformal grib decoder and grid coordinate
!! transformations.
!!
!! Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.30
!! Lambert
!! conformal](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-30.shtml).
!!
!> @author Iredell @date 96-04-10
module ip_lambert_conf_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  use ip_constants_mod
  implicit none

  private
  public :: ip_lambert_conf_grid

  type,extends(ip_grid) :: ip_lambert_conf_grid
    real :: rlat1 !< La1― latitude of first grid point. GRIB2, Section 3.30, octet 39-42.
    real :: rlon1 !< Lo1― longitude of first grid point. GRIB2, Section 3.30, octet 43-46.
    real :: rlati1 !< First latitude from the pole at which the secant cone cuts the sphere. GRIB2, Section 3, octets 66-69.
    real :: rlati2 !< Second latitude from the pole at which the scant cone cuts the sphere. GRIB2, Section 3, octets 70-73.
    real :: orient !<  Longitude of meridian parallel to y-axis along which latitude increases at the latitude increases. GRIB2, Section 3, octets 52-55.
    real :: dxs !< x-direction grid length adjusted for scan mode. GRIB2, Section 3, octets 56-59.
    real :: dys !< y-direction grid length adjusted for scan model. GRIB2, Section 3, octets 60-63.
    real :: h !< Hemisphere flag. 1-NH, minus 1-SH.
    integer :: irot !< vector rotation flag. When "1", vectors are grid relative. When "0", vectors are earth relative. GRIB2, Section 3, octet 55.
  contains
    !> Initializes a gaussian grid given a grib1_descriptor object. @return N/A
    procedure :: init_grib1
    !> Initializes a gaussian grid given a grib2_descriptor object. @return N/A
    procedure :: init_grib2
    !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
    !> for Gaussian grids. @return N/A
    procedure :: gdswzd=>gdswzd_lambert_conf
  endtype ip_lambert_conf_grid

  integer :: irot !< vector rotation flag. When "1", vectors are grid relative. When "0", vectors are earth relative. GRIB2, Section 3, octet 55.
  real :: an !< Cone factor
  real :: dxs !< x-direction grid length adjusted for scan mode. GRIB2, Section 3, octets 56-59.
  real :: dys !< y-direction grid length adjusted for scan model. GRIB2, Section 3, octets 60-63.
  real :: h !<  Hemisphere flag. 1-NH, minus 1-SH.
  real :: rerth !< Radius of the earth. GRIB2, Section 3, octets 15-30.
  real :: tinyreal=tiny(1.0) !< Smallest positive real value (use for equality comparisons)

contains

  !> Initializes a Lambert Conformal grid given a grib1_descriptor object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Iredell @date 96-04-10
  subroutine init_grib1(self,g1_desc)
    class(ip_lambert_conf_grid),intent(inout) :: self
    type(grib1_descriptor),intent(in) :: g1_desc

    real :: dx,dy,hi,hj
    integer :: iproj,iscan,jscan

    associate(kgds=>g1_desc%gds)
      self%rerth=6.3712e6
      self%eccen_squared=0.0

      self%im=kgds(2)
      self%jm=kgds(3)

      self%rlat1=kgds(4)*1.e-3
      self%rlon1=kgds(5)*1.e-3

      self%irot=mod(kgds(6)/8,2)
      self%orient=kgds(7)*1.e-3

      dx=kgds(8)
      dy=kgds(9)

      iproj=mod(kgds(10)/128,2)
      iscan=mod(kgds(11)/128,2)
      jscan=mod(kgds(11)/64,2)

      self%rlati1=kgds(12)*1.e-3
      self%rlati2=kgds(13)*1.e-3
      self%h=(-1.)**iproj

      hi=(-1.)**iscan
      hj=(-1.)**(1-jscan)
      self%dxs=dx*hi
      self%dys=dy*hj

      self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0
      self%nscan=mod(kgds(11)/32,2)
      self%nscan_field_pos=self%nscan
      self%kscan=0
    endassociate

  endsubroutine init_grib1

  !> Initializes a Lambert Conformal grid given a grib2_descriptor object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g2_desc A grib2_descriptor
  !!
  !! @author Iredell @date 96-04-10
  subroutine init_grib2(self,g2_desc)
    class(ip_lambert_conf_grid),intent(inout) :: self
    type(grib2_descriptor),intent(in) :: g2_desc

    real :: dx,dy,hi,hj
    integer :: iproj,iscan,jscan

    associate(igdtmpl=>g2_desc%gdt_tmpl,igdtlen=>g2_desc%gdt_len)
      call earth_radius(igdtmpl,igdtlen,self%rerth,self%eccen_squared)

      self%im=igdtmpl(8)
      self%jm=igdtmpl(9)

      self%rlat1=float(igdtmpl(10))*1.0e-6
      self%rlon1=float(igdtmpl(11))*1.0e-6

      self%irot=mod(igdtmpl(12)/8,2)
      self%orient=float(igdtmpl(14))*1.0e-6

      dx=float(igdtmpl(15))*1.0e-3
      dy=float(igdtmpl(16))*1.0e-3

      iproj=mod(igdtmpl(17)/128,2)
      iscan=mod(igdtmpl(18)/128,2)
      jscan=mod(igdtmpl(18)/64,2)

      self%rlati1=float(igdtmpl(19))*1.0e-6
      self%rlati2=float(igdtmpl(20))*1.0e-6

      self%h=(-1.)**iproj
      hi=(-1.)**iscan
      hj=(-1.)**(1-jscan)
      self%dxs=dx*hi
      self%dys=dy*hj

      self%nscan=mod(igdtmpl(18)/32,2)
      self%nscan_field_pos=self%nscan
      self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0
      self%kscan=0
    endassociate
  endsubroutine init_grib2

  !> GDS wizard for lambert conformal conical.
  !>
  !> This subprogram decodes the grib 2 grid definition template
  !> (passed in integer form as decoded by the ncep g2 library) and
  !> returns one of the following:
  !> - (iopt=+1) earth coordinates of selected grid coordinates
  !> - (iopt=-1) grid coordinates of selected earth coordinates
  !>
  !> Works for lambert conformal conical projections.
  !>
  !> If the selected coordinates are more than one gridpoint beyond
  !> the the edges of the grid domain, then the relevant output
  !> elements are set to fill values.
  !>
  !> The actual number of valid points computed is returned too.
  !>
  !> Optionally, the vector rotations, map jacobians and grid box
  !> areas for this grid may be returned as well.
  !>
  !> To compute the vector rotations, the optional arguments 'srot'
  !> and 'crot' must be present. To compute the map jacobians, the
  !> optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be
  !> present. To compute the grid box areas the optional argument
  !> 'area' must be present.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 96-04-10 | iredell | Initial.
  !> 96-10-01 | iredell | protected against unresolvable points
  !> 97-10-20 | iredell | include map options
  !> 1999-04-27 | gilbert | corrected minor error calculating variable an for the secant projection case (rlati1.ne.rlati2).
  !> 2012-08-14 | gayno | fix problem with sh grids. Ensure grid box area always positive.
  !> 2015-01-21 | gayno | merger of gdswiz03() and gdswzd03(). Make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. Make part of a module. Move vector rotation, map jacobian and grid box area computations to separate subroutines.
  !> 2015-07-13 | gayno | Convert to grib 2. Replace grib 1 kgds array with grib 2 grid definition template array. Rename routine.
  !> 2018-07-20 | wesley | add threads.
  !>
  !> @param[in] self ip_lambert_conf_grid object.
  !> @param[in] iopt option flag
  !> - 1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param[in] npts maximum number of coordinates
  !> @param[in] fill fill value to set invalid output data (must be
  !> impossible value; suggested value: -9999.)
  !> @param[inout] xpts (npts) grid x point coordinates if iopt>0
  !> @param[inout] ypts (npts) grid y point coordinates if iopt>0
  !> @param[inout] rlon (npts) earth longitudes in degrees e if iopt<0
  !> (acceptable range: -360. to 360.)
  !> @param[inout] rlat (npts) earth latitudes in degrees n if iopt<0
  !> (acceptable range: -90. to 90.)
  !> @param[out] nret number of valid points computed
  !> @param[out] crot optional (npts) clockwise vector rotation cosines
  !> @param[out] srot optional (npts) clockwise vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !> @param[out] xlon optional (npts) dx/dlon in 1/degrees
  !> @param[out] xlat optional (npts) dx/dlat in 1/degrees
  !> @param[out] ylon optional (npts) dy/dlon in 1/degrees
  !> @param[out] ylat optional (npts) dy/dlat in 1/degrees
  !> @param[out] area optional (npts) area weights in m**2
  !> (proportional to the square of the map factor)
  !>
  !> @author Iredell @date 96-04-10
  subroutine gdswzd_lambert_conf(self,iopt,npts,fill, &
                                 xpts,ypts,rlon,rlat,nret, &
                                 crot,srot,xlon,xlat,ylon,ylat,area)
    implicit none
    !
    class(ip_lambert_conf_grid),intent(in) :: self
    integer,intent(in) :: iopt,npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon(npts),rlat(npts)
    real,intent(inout) :: xpts(npts),ypts(npts)
    real,optional,intent(out) :: crot(npts),srot(npts)
    real,optional,intent(out) :: xlon(npts),xlat(npts)
    real,optional,intent(out) :: ylon(npts),ylat(npts),area(npts)
    !
    integer                       :: im,jm,n
    !
    logical                       :: lrot,lmap,larea
    !
    real                          :: antr,di,dj
    real                          :: dlon1
    real                          :: de,de2,dr2
    real                          :: orient,rlat1,rlon1
    real                          :: rlati1,rlati2
    real                          :: xmax,xmin,ymax,ymin,xp,yp
    real                          :: dlon,dr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(present(crot)) crot=fill
    if(present(srot)) srot=fill
    if(present(xlon)) xlon=fill
    if(present(xlat)) xlat=fill
    if(present(ylon)) ylon=fill
    if(present(ylat)) ylat=fill
    if(present(area)) area=fill

    im=self%im
    jm=self%jm

    rlat1=self%rlat1
    rlon1=self%rlon1

    irot=self%irot
    orient=self%orient

    rlati1=self%rlati1
    rlati2=self%rlati2

    h=self%h
    dxs=self%dxs
    dys=self%dys

    rerth=self%rerth

    if(abs(rlati1-rlati2).lt.tinyreal) then
      an=sin(rlati1/dpr)
    else
      an=log(cos(rlati1/dpr)/cos(rlati2/dpr))/ &
          log(tan((90-rlati1)/2/dpr)/tan((90-rlati2)/2/dpr))
    endif
    de=rerth*cos(rlati1/dpr)*tan((rlati1+90)/2/dpr)**an/an
    if(abs(h*rlat1-90).lt.tinyreal) then
      xp=1
      yp=1
    else
      dr=de/tan((rlat1+90)/2/dpr)**an
      dlon1=mod(rlon1-orient+180+3600,360.)-180
      xp=1-sin(an*dlon1/dpr)*dr/dxs
      yp=1+cos(an*dlon1/dpr)*dr/dys
    endif
    antr=1/(2*an)
    de2=de**2
    xmin=0
    xmax=im+1
    ymin=0
    ymax=jm+1
    nret=0
    if(present(crot).and.present(srot)) then
      lrot=.true.
    else
      lrot=.false.
    endif
    if(present(xlon).and.present(xlat).and.present(ylon).and.present(ylat)) then
      lmap=.true.
    else
      lmap=.false.
    endif
    if(present(area)) then
      larea=.true.
    else
      larea=.false.
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    if(iopt.eq.0.or.iopt.eq.1) then
      !$omp parallel do private(n,di,dj,dr2,dr,dlon) reduction(+:nret) schedule(static)
      do n=1,npts
        if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
           ypts(n).ge.ymin.and.ypts(n).le.ymax) then
          di=h*(xpts(n)-xp)*dxs
          dj=h*(ypts(n)-yp)*dys
          dr2=di**2+dj**2
          dr=sqrt(dr2)
          if(dr2.lt.de2*1.e-6) then
            rlon(n)=0.
            rlat(n)=h*90.
          else
            rlon(n)=mod(orient+1./an*dpr*atan2(di,-dj)+3600,360.)
            rlat(n)=(2*dpr*atan((de2/dr2)**antr)-90)
          endif
          nret=nret+1
          dlon=mod(rlon(n)-orient+180+3600,360.)-180
          if(lrot) call lambert_conf_vect_rot(dlon,crot(n),srot(n))
          if(lmap) call lambert_conf_map_jacob(rlat(n),fill,dlon,dr, &
                                               xlon(n),xlat(n),ylon(n),ylat(n))
          if(larea) call lambert_conf_grid_area(rlat(n),fill,dr,area(n))
        else
          rlon(n)=fill
          rlat(n)=fill
        endif
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    elseif(iopt.eq.-1) then
      !$omp parallel do private(n,dr,dlon) reduction(+:nret) schedule(static)
      do n=1,npts
        if(abs(rlon(n)).lt.(360.+tinyreal).and.abs(rlat(n)).lt.(90.+tinyreal).and. &
           abs(h*rlat(n)+90).gt.tinyreal) then
          dr=h*de*tan((90-rlat(n))/2/dpr)**an
          dlon=mod(rlon(n)-orient+180+3600,360.)-180
          xpts(n)=xp+h*sin(an*dlon/dpr)*dr/dxs
          ypts(n)=yp-h*cos(an*dlon/dpr)*dr/dys
          if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
             ypts(n).ge.ymin.and.ypts(n).le.ymax) then
            nret=nret+1
            if(lrot) call lambert_conf_vect_rot(dlon,crot(n),srot(n))
            if(lmap) call lambert_conf_map_jacob(rlat(n),fill,dlon,dr, &
                                                 xlon(n),xlat(n),ylon(n),ylat(n))
            if(larea) call lambert_conf_grid_area(rlat(n),fill,dr,area(n))
          else
            xpts(n)=fill
            ypts(n)=fill
          endif
        else
          xpts(n)=fill
          ypts(n)=fill
        endif
      enddo
      !$omp end parallel do
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  endsubroutine gdswzd_lambert_conf

  !> Vector rotation fields for lambert conformal conical.
  !>
  !> This subprogram computes the vector rotation sines and
  !> cosines for a lambert conformal conical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "lambert_conf_vect_rot"
  !> 2018-07-20 | wesley | pass in dlon for threading.
  !>
  !> @param[in] dlon from orientation longitude (real)
  !> @param[out] crot vector rotation cosines (real)
  !> @param[out] srot vector rotation sines (real)
  !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
  !>
  !> @author Gayno @date 2015-01-21
  subroutine lambert_conf_vect_rot(dlon,crot,srot)
    implicit none
    real,intent(in) :: dlon
    real,intent(out) :: crot,srot

    if(irot.eq.1) then
      crot=cos(an*dlon/dpr)
      srot=sin(an*dlon/dpr)
    else
      crot=1.
      srot=0.
    endif

  endsubroutine lambert_conf_vect_rot

  !> Map jacobians for lambert conformal conical.
  !>
  !> This subprogram computes the map jacobians for a lambert
  !> conformal conical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | Gayno | initial version
  !> 2015-09-17 | Gayno | rename as "lambert_conf_map_jacob"
  !> 2018-07-20 | Wesley | pass dlon and dr for threading.
  !>
  !> @param[in] rlat grid point latitude in degrees (real)
  !> @param[in] fill fill value for undefined points (real)
  !> @param[in] dlon distance from orientation longitude (real)
  !> @param[in] dr distance from pole point (real)
  !> @param[out] xlon dx/dlon in 1/degrees (real)
  !> @param[out] xlat dx/dlat in 1/degrees (real)
  !> @param[out] ylon dy/dlon in 1/degrees (real)
  !> @param[out] ylat dy/dlat in 1/degrees (real)
  !>
  !> @author Gayno @date 2015-01-21
  subroutine lambert_conf_map_jacob(rlat,fill,dlon,dr,xlon,xlat,ylon,ylat)
    implicit none

    real,intent(in) :: rlat,fill,dlon,dr
    real,intent(out) :: xlon,xlat,ylon,ylat

    real                          :: clat

    clat=cos(rlat/dpr)
    if(clat.le.0.or.dr.le.0) then
      xlon=fill
      xlat=fill
      ylon=fill
      ylat=fill
    else
      xlon=h*cos(an*dlon/dpr)*an/dpr*dr/dxs
      xlat=-h*sin(an*dlon/dpr)*an/dpr*dr/dxs/clat
      ylon=h*sin(an*dlon/dpr)*an/dpr*dr/dys
      ylat=h*cos(an*dlon/dpr)*an/dpr*dr/dys/clat
    endif

  endsubroutine lambert_conf_map_jacob

  !> Grid box area for lambert conformal conical.
  !>
  !> This subprogram computes the grid box area for a lambert
  !> conformal conical grid.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | Gayno | initial version
  !> 2015-09-17 | Gayno | rename as "lambert_conf_grid_area"
  !> 2018-07-20 | Wesley | pass in dr for threading.
  !>
  !> @param[in] rlat latitude of grid point in degrees (real)
  !> @param[in] fill fill value for undefined points (real)
  !> @param[in] dr distance from pole point (real)
  !> @param[out] area area weights in m**2 (real)
  !>
  !> @author Gayno @date 2015-01-21
  subroutine lambert_conf_grid_area(rlat,fill,dr,area)
    implicit none

    real,intent(in) :: rlat
    real,intent(in) :: fill
    real,intent(in) :: dr
    real,intent(out) :: area

    real                          :: clat

    clat=cos(rlat/dpr)
    if(clat.le.0.or.dr.le.0) then
      area=fill
    else
      area=rerth**2*clat**2*abs(dxs)*abs(dys)/(an*dr)**2
    endif

  endsubroutine lambert_conf_grid_area

endmodule ip_lambert_conf_grid_mod

