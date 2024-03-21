!> @file
!! @brief Equidistant cylindrical grib decoder and grid coordinate
!! transformations.
!!
!! @author Mark Iredell, George Gayno, Kyle Gerheiser
!! @date July 2021

!> @brief Equidistant cylindrical grib decoder and grid coordinate
!! transformations.
!!
!! Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.0
!! Latitude/Longitude or equidistant cylindrical, or Plate
!! Carree](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-0.shtml).
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
!! @date July 2021
module ip_equid_cylind_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  implicit none

  private
  public :: ip_equid_cylind_grid

  type,extends(ip_grid) :: ip_equid_cylind_grid
    real :: hi !< Scan mode in the 'i' direction. GRIB2, Section 3, octet 72.
    real :: rlat1 !< Latitude of first grid point. GRIB2, Section 3, octets 47-50.
    real :: rlon1 !< Longitude of first grid point. GRIB2, Section 3, octets 51-54.
    real :: rlat2 !< Latitude of last grid point. GRIB2, Section 3, octets 56-59.
    real :: rlon2 !< Longitude of last grid point. GRIB2, Section 3, octets 60-63.
    real :: dlat !< Di — i direction increment. GRIB2, Section 3, octets 64-67.
    real :: dlon !< Dj — j direction increment. GRIB2, Section 3, octets 68-71.
  contains
    procedure :: init_grib1 !< Init GRIB1. @return N/A
    procedure :: init_grib2 !< Init GRIB2. @return N/A
    procedure :: gdswzd=>gdswzd_equid_cylind !< See gdswzd_equid_cylind(). @return N/A
  endtype ip_equid_cylind_grid

  real :: dlat !< Grid resolution in degrees n/s direction.
  real :: dlon !< Grid resolution in degrees e/w direction.
  real :: rerth !< Radius of the Earth.

contains

  !> Initializes an equidistant cylindrical grid given a grib1_descriptor object.
  !!
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  subroutine init_grib1(self,g1_desc)
    class(ip_equid_cylind_grid),intent(inout) :: self
    type(grib1_descriptor),intent(in) :: g1_desc

    integer :: iscan

    associate(kgds=>g1_desc%gds)
      self%im=kgds(2)
      self%jm=kgds(3)
      self%rlat1=kgds(4)*1.e-3
      self%rlon1=kgds(5)*1.e-3
      self%rlat2=kgds(7)*1.e-3
      self%rlon2=kgds(8)*1.e-3
      iscan=mod(kgds(11)/128,2)
      self%hi=(-1.)**iscan
      self%dlon=self%hi*(mod(self%hi*(self%rlon2-self%rlon1)-1+3600,360.)+1)/(self%im-1)
      self%dlat=(self%rlat2-self%rlat1)/(self%jm-1)

      ! defaults
      self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0
      self%nscan=mod(kgds(11)/32,2)
      self%nscan_field_pos=self%nscan
      self%kscan=0

      self%iwrap=nint(360/abs(self%dlon))

      if(self%im.lt.self%iwrap) self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0
      if(self%iwrap.gt.0.and.mod(self%iwrap,2).eq.0) then
        if(abs(self%rlat1).gt.90-0.25*self%dlat) then
          self%jwrap1=2
        elseif(abs(self%rlat1).gt.90-0.75*self%dlat) then
          self%jwrap1=1
        endif
        if(abs(self%rlat2).gt.90-0.25*self%dlat) then
          self%jwrap2=2*self%jm
        elseif(abs(self%rlat2).gt.90-0.75*self%dlat) then
          self%jwrap2=2*self%jm+1
        endif
      endif

      self%rerth=6.3712e6
      self%eccen_squared=0.0
    endassociate

  endsubroutine init_grib1

  !> Initializes an equidistant cylindrical grid given a grib2_descriptor object.
  !! @param[inout] self The grid to initialize
  !! @param[in] g2_desc A grib2_descriptor
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  subroutine init_grib2(self,g2_desc)
    class(ip_equid_cylind_grid),intent(inout) :: self
    type(grib2_descriptor),intent(in) :: g2_desc

    integer :: iscale,iscan

    associate(igdtmpl=>g2_desc%gdt_tmpl,igdtlen=>g2_desc%gdt_len)
      self%im=igdtmpl(8)
      self%jm=igdtmpl(9)
      iscale=igdtmpl(10)*igdtmpl(11)
      if(iscale.eq.0) iscale=10**6
      self%rlat1=float(igdtmpl(12))/float(iscale)
      self%rlon1=float(igdtmpl(13))/float(iscale)
      self%rlat2=float(igdtmpl(15))/float(iscale)
      self%rlon2=float(igdtmpl(16))/float(iscale)
      iscan=mod(igdtmpl(19)/128,2)
      self%hi=(-1.)**iscan
      self%dlon=self%hi*(mod(self%hi*(self%rlon2-self%rlon1)-1+3600,360.)+1)/(self%im-1)
      self%dlat=(self%rlat2-self%rlat1)/(self%jm-1)

      self%nscan=mod(igdtmpl(19)/32,2)
      self%nscan_field_pos=self%nscan
      self%kscan=0
      self%iwrap=nint(360/abs(self%dlon))

      if(self%im.lt.self%iwrap) self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0

      if(self%im.lt.self%iwrap) self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0
      if(self%iwrap.gt.0.and.mod(self%iwrap,2).eq.0) then
        if(abs(self%rlat1).gt.90-0.25*self%dlat) then
          self%jwrap1=2
        elseif(abs(self%rlat1).gt.90-0.75*self%dlat) then
          self%jwrap1=1
        endif
        if(abs(self%rlat2).gt.90-0.25*self%dlat) then
          self%jwrap2=2*self%jm
        elseif(abs(self%rlat2).gt.90-0.75*self%dlat) then
          self%jwrap2=2*self%jm+1
        endif
      endif

      call earth_radius(igdtmpl,igdtlen,self%rerth,self%eccen_squared)

    endassociate
  endsubroutine init_grib2

  !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
  !! for equidistant cylindrical grids.
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !! Optionally, the vector rotations, the map jacobians and
  !! the grid box areas may be returned as well.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present.
  !!
  !! To compute the map jacobians, the optional arguments
  !! 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! @param[in] self The grid object gdswzd was called on.
  !! @param[in] iopt option flag
  !!            - +1 to compute earth coords of selected grid coords.
  !!            - -1 o compute grid coords of selected earth coords.
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!            Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates if iopt>0.
  !! @param[inout] ypts Grid y point coordinates if iopt>0.
  !! @param[inout] rlon Earth longitudes in degrees e if iopt<0
  !!                   (Acceptable range: -360. to 360.)
  !! @param[inout] rlat Earth latitudes in degrees n if iopt<0
  !!                (Acceptable range: -90. to 90.)
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Optional clockwise vector rotation cosines.
  !! @param[out] srot Optional clockwise vector rotation sines.
  !! @param[out] xlon Optional dx/dlon in 1/degrees.
  !! @param[out] xlat Optional dx/dlat in 1/degrees.
  !! @param[out] ylon Optional dy/dlon in 1/degrees.
  !! @param[out] ylat Optional dy/dlat in 1/degrees.
  !! @param[out] area Optional area weights in m**2.
  !!
  !! @author Mark Iredell, George Gayno, Kyle Gerheiser
  !! @date July 2021
  subroutine gdswzd_equid_cylind(self,iopt,npts,fill, &
                                 xpts,ypts,rlon,rlat,nret, &
                                 crot,srot,xlon,xlat,ylon,ylat,area)
    implicit none
    !
    class(ip_equid_cylind_grid),intent(in) :: self
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
    integer                            :: im,jm,n
    !
    logical                            :: lrot,lmap,larea
    !
    real                               :: hi,rlat1,rlon1,rlat2,rlon2
    real                               :: xmax,xmin,ymax,ymin
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
    rlat2=self%rlat2
    rlon2=self%rlon2

    hi=self%hi

    rerth=self%rerth
    dlat=self%dlat
    dlon=self%dlon

    xmin=0
    xmax=im+1
    if(im.eq.nint(360/abs(dlon))) xmax=im+2
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
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    if(iopt.eq.0.or.iopt.eq.1) then
      !$omp parallel do private(n) reduction(+:nret) schedule(static)
      do n=1,npts
        if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
           ypts(n).ge.ymin.and.ypts(n).le.ymax) then
          rlon(n)=mod(rlon1+dlon*(xpts(n)-1)+3600,360.)
          rlat(n)=min(max(rlat1+dlat*(ypts(n)-1),-90.),90.)
          nret=nret+1
          if(lrot) call equid_cylind_vect_rot(crot(n),srot(n))
          if(lmap) call equid_cylind_map_jacob(xlon(n),xlat(n),ylon(n),ylat(n))
          if(larea) call equid_cylind_grid_area(rlat(n),area(n))
        else
          rlon(n)=fill
          rlat(n)=fill
        endif
      enddo
      !$omp end parallel do
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    elseif(iopt.eq.-1) then
      !$omp parallel do private(n) reduction(+:nret) schedule(static)
      do n=1,npts
        if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90) then
          xpts(n)=1+hi*mod(hi*(rlon(n)-rlon1)+3600,360.)/dlon
          ypts(n)=1+(rlat(n)-rlat1)/dlat
          if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
             ypts(n).ge.ymin.and.ypts(n).le.ymax) then
            nret=nret+1
            if(lrot) call equid_cylind_vect_rot(crot(n),srot(n))
            if(lmap) call equid_cylind_map_jacob(xlon(n),xlat(n),ylon(n),ylat(n))
            if(larea) call equid_cylind_grid_area(rlat(n),area(n))
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
  endsubroutine gdswzd_equid_cylind

  !> Computes the vector rotation sines and
  !! cosines for a equidistant cylindrical grid.
  !!
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !!
  !! @note
  !! - ugrid=crot*uearth-srot*vearth;
  !! - vgrid=srot*uearth+crot*vearth
  !!
  !! @author George Gayno
  !! @date July 2021
  subroutine equid_cylind_vect_rot(crot,srot)
    implicit none

    real,intent(out) :: crot,srot

    crot=1.0
    srot=0.0

  endsubroutine equid_cylind_vect_rot

  !> Computes the map jacobians for a equidistant cylindrical grid.
  !!
  !! @param[out] xlon dx/dlon in 1/degrees.
  !! @param[out] xlat dx/dlat in 1/degrees.
  !! @param[out] ylon dy/dlon in 1/degrees.
  !! @param[out] ylat dy/dlat in 1/degrees.
  !!
  !! @author George Gayno
  !! @date July 2021
  subroutine equid_cylind_map_jacob(xlon,xlat,ylon,ylat)
    real,intent(out) :: xlon,xlat,ylon,ylat

    xlon=1.0/dlon
    xlat=0.
    ylon=0.
    ylat=1.0/dlat

  endsubroutine equid_cylind_map_jacob

  !> Computes the grid box area for a equidistant cylindrical grid.
  !!
  !! @param[in] rlat Latitude of grid point in degrees.
  !! @param[out] area Area weights in m^2.
  !!
  !! @author Mark Iredell, George Gayno
  !! @date July 2021
  subroutine equid_cylind_grid_area(rlat,area)
    implicit none

    real,intent(in) :: rlat
    real,intent(out) :: area

    real,parameter     :: pi=3.14159265358979
    real,parameter     :: dpr=180./pi

    real                               :: dslat,rlatu,rlatd

    rlatu=min(max(rlat+dlat/2,-90.),90.)
    rlatd=min(max(rlat-dlat/2,-90.),90.)
    dslat=sin(rlatu/dpr)-sin(rlatd/dpr)
    area=rerth**2*abs(dslat*dlon)/dpr

  endsubroutine equid_cylind_grid_area

endmodule ip_equid_cylind_grid_mod

