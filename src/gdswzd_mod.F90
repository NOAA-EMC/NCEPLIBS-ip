!> @file
!! @brief Driver module for gdswzd routines.
!!
!! @date Jan 2015
!! @author George Gayno, Kyle Gerheiser

!> @brief Driver module for gdswzd routines.
!!
!! These routines do the following for several map projections:
!! - Convert from earth to grid coordinates or vice versa.
!! - Compute vector rotation sines and cosines.
!! - Compute map jacobians.
!! - Compute grid box area.
!!
!! Map projections include:
!! - Equidistant Cyclindrical
!! - Mercator Cylindrical
!! - Gaussian Cylindrical
!! - Polar stereographic
!! - Lambert Conformal Conic
!! - Rotated Equidistant Cyclindrical ("E" and non-"E" staggers)
!!
!! @author Mark Iredell, George Gayno, Kyle Gerheiser
!! @date Jan 2015
module gdswzd_mod
  use ip_grid_descriptor_mod
  use ip_grids_mod
  use ip_grid_mod
  use ip_grid_factory_mod

  implicit none

  private

  public :: gdswzd_2d_array_grib1,gdswzd_grib1,gdswzd

  interface gdswzd
    module procedure gdswzd_1d_array
    module procedure gdswzd_2d_array
    module procedure gdswzd_scalar
    module procedure gdswzd_grib1
    module procedure gdswzd_2d_array_grib1
    module procedure gdswzd_grid
  endinterface gdswzd

contains

  !> Returns one of the following for a grid object:
  !! - iopt=0 Grid and earth coordinates of all grid points.
  !! - iopt=+1 Earth coordinates of selected grid coordinates.
  !! - iopt=-1 Grid coordinates of selected earth coordinates.
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values.  Also if iopt=0,
  !! if the number of grid points exceeds the number allotted,
  !! then all the output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !!
  !! Optionally, the vector rotations, map jacobians and
  !! grid box areas may be returned.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present. To compute the map jacobians, the
  !! optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! @param[in] grid Grid to call gdswzd on.
  !!
  !! @param[in] iopt Option flag.
  !! - 0 Earth coords of all the grid points.
  !! - 1 Earth coords of selected grid coords.
  !! - -1 Grid coords of selected earth coords
  !!
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!      Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates.
  !! @param[inout] ypts Grid y point coordinates.
  !! @param[inout] rlon Earth longitudes in degrees E.
  !! (Acceptable range: -360. to 360.)
  !!
  !! @param[inout] rlat Earth latitudes in degrees N.
  !! (Acceptable range: -90. to 90.)
  !!
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @param[out] xlon dx/dlat in 1/degrees
  !! @param[out] xlat dy/dlon in 1/degrees
  !! @param[out] ylon dy/dlon in 1/degrees
  !! @param[out] ylat dy/dlat in 1/degrees
  !! @param[out] area Area weights in m^2.
  !! Proportional to the square of the map factor in the case of
  !! conformal projections
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  subroutine gdswzd_grid(grid,iopt,npts,fill, &
                         xpts,ypts,rlon,rlat,nret, &
                         crot,srot,xlon,xlat,ylon,ylat,area)

    class(ip_grid),intent(in) :: grid
    integer,intent(in) :: iopt,npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon(npts),rlat(npts)
    real,intent(inout) :: xpts(npts),ypts(npts)
    real,optional,intent(out) :: crot(npts),srot(npts)
    real,optional,intent(out) :: xlon(npts),xlat(npts)
    real,optional,intent(out) :: ylon(npts),ylat(npts),area(npts)

    integer                       :: is1,im,jm,nm,kscan,nscan,n
    integer                       :: iopf,nn,i,j

    !  COMPUTE GRID COORDINATES FOR ALL GRID POINTS
    if(iopt.eq.0) then
      iopf=1

      if(grid%descriptor%grid_num.eq.-1) then
        nm=npts
      else
        im=grid%im
        jm=grid%jm
        nm=im*jm
      endif
      nscan=grid%nscan
      kscan=grid%kscan

      if(nm.gt.npts) then
        rlat=fill
        rlon=fill
        xpts=fill
        ypts=fill
        return
      endif

      select type(grid)
      type is(ip_rot_equid_cylind_egrid)
        if(kscan.eq.0) then
          is1=(jm+1)/2
        else
          is1=jm/2
        endif

        do n=1,nm
          if(nscan.eq.0) then
            j=(n-1)/im+1
            i=(n-im*(j-1))*2-mod(j+kscan,2)
          else
            nn=(n*2)-1+kscan
            i=(nn-1)/jm+1
            j=mod(nn-1,jm)+1
            if(mod(jm,2).eq.0.and.mod(i,2).eq.0.and.kscan.eq.0) j=j+1
            if(mod(jm,2).eq.0.and.mod(i,2).eq.0.and.kscan.eq.1) j=j-1
          endif
          xpts(n)=is1+(i-(j-kscan))/2
          ypts(n)=(i+(j-kscan))/2
        enddo
      type is(ip_station_points_grid)
        do n=1,nm
          xpts(n)=fill
          ypts(n)=fill
        enddo
      class default
        do n=1,nm
          if(nscan.eq.0) then
            j=(n-1)/im+1
            i=n-im*(j-1)
          else
            i=(n-1)/jm+1
            j=n-jm*(i-1)
          endif
          xpts(n)=i
          ypts(n)=j
        enddo
      endselect

      do n=nm+1,npts
        xpts(n)=fill
        ypts(n)=fill
      enddo

    else  ! IOPT /= 0
      iopf=iopt
    endif ! IOPT CHECK

    call grid%gdswzd(iopf,npts,fill, &
                     xpts,ypts,rlon,rlat,nret, &
                     crot,srot,xlon,xlat,ylon,ylat,area)

  endsubroutine gdswzd_grid

  !> Decodes the grib 2 grid definition template and returns
  !! one of the following (for scalars):
  !! - iopt=0 Grid and earth coordinates of all grid points.
  !! - iopt=+1 Earth coordinates of selected grid coordinates.
  !! - iopt=-1 Grid coordinates of selected earth coordinates.
  !!
  !! The current code recognizes the following projections,
  !! where "igdtnum" is the grid definition template number:
  !! - igdtnum=00 Equidistant Cylindrical
  !! - igdtnum=01 Rotated Equidistant Cylindrical. "E" and non-"E" staggered
  !! - igdtnum=10 Mercator Cyclindrical
  !! - igdtnum=20 Polar Stereographic Azimuthal
  !! - igdtnum=30 Lambert Conformal Conical
  !! - igdtnum=40 Gaussian Equidistant Cyclindrical
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values. Also if iopt=0,
  !! if the number of grid points exceeds the number allotted,
  !! then all the output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !!
  !! Optionally, the vector rotations, map jacobians and
  !! grid box areas may be returned.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present. To compute the map jacobians, the
  !! optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! @param[in] igdtnum Grid definition template number.
  !!
  !! @param[in] igdtmpl Grid definition template array.
  !! Corresponds to the gfld%igdtmpl component of the
  !! NCEPLIBS-g2 gridmod data structure
  !! See igdtmpl definition in gdswzd_1d_array() for full details.
  !!
  !! @param[in] igdtlen Number of elements of the grid definition
  !! template array. Corresponds to the gfld%igdtlen
  !! component of the ncep g2 library gridmod data structure.
  !!
  !! @param[in] iopt Option flag.
  !! - 0 Earth coords of all the grid points.
  !! - 1 Earth coords of selected grid coords.
  !! - -1 Grid coords of selected earth coords
  !!
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!      Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates.
  !! @param[inout] ypts Grid y point coordinates.
  !! @param[inout] rlon Earth longitudes in degrees E.
  !! (Acceptable range: -360. to 360.)
  !!
  !! @param[inout] rlat Earth latitudes in degrees N.
  !! (Acceptable range: -90. to 90.)
  !!
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @param[out] xlon dx/dlat in 1/degrees
  !! @param[out] xlat dy/dlon in 1/degrees
  !! @param[out] ylon dy/dlon in 1/degrees
  !! @param[out] ylat dy/dlat in 1/degrees
  !! @param[out] area Area weights in m^2.
  !! Proportional to the square of the map factor in the case of
  !! conformal projections
  !!
  !! @author George Gayno, Mark Iredell
  !! @date Jan 2015
  subroutine gdswzd_scalar(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                           xpts,ypts,rlon,rlat,nret, &
                           crot,srot,xlon,xlat,ylon,ylat,area)

    implicit none
    !
    integer,intent(in) :: igdtnum,igdtlen
    integer,intent(in) :: igdtmpl(igdtlen)
    integer,intent(in) :: iopt,npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon,rlat
    real,intent(inout) :: xpts,ypts
    real,optional,intent(out) :: crot,srot
    real,optional,intent(out) :: xlon,xlat
    real,optional,intent(out) :: ylon,ylat,area

    real                          :: rlona(1),rlata(1)
    real                          :: xptsa(1),yptsa(1)
    real                          :: crota(1),srota(1)
    real                          :: xlona(1),xlata(1)
    real                          :: ylona(1),ylata(1),areaa(1)

    rlona(1)=rlon
    rlata(1)=rlat
    xptsa(1)=xpts
    yptsa(1)=ypts

    nret=0

    ! CALL WITHOUT EXTRA FIELDS.

    if(.not.present(crot).and. &
       .not.present(srot).and. &
       .not.present(xlon).and. &
       .not.present(xlat).and. &
       .not.present(ylon).and. &
       .not.present(ylat).and. &
       .not.present(area)) then

      call gdswzd_1d_array(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                           xptsa,yptsa,rlona,rlata,nret)

      rlon=rlona(1)
      rlat=rlata(1)
      xpts=xptsa(1)
      ypts=yptsa(1)

    endif

    ! MIMIC CALL TO OLD 'GDSWIZ' ROUTINES.

    if(present(crot).and. &
       present(srot).and. &
       .not.present(xlon).and. &
       .not.present(xlat).and. &
       .not.present(ylon).and. &
       .not.present(ylat).and. &
       .not.present(area)) then

      call gdswzd_1d_array(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                           xptsa,yptsa,rlona,rlata,nret,crota,srota)

      rlon=rlona(1)
      rlat=rlata(1)
      xpts=xptsa(1)
      ypts=yptsa(1)
      crot=crota(1)
      srot=srota(1)

    endif

    ! MIMIC CALL TO OLD 'GDSWZD' ROUTINES.

    if(present(crot).and. &
       present(srot).and. &
       present(xlon).and. &
       present(xlat).and. &
       present(ylon).and. &
       present(ylat).and. &
       present(area)) then

      call gdswzd_1d_array(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                           xptsa,yptsa,rlona,rlata,nret, &
                           crota,srota,xlona,xlata,ylona,ylata,areaa)

      rlon=rlona(1)
      rlat=rlata(1)
      xpts=xptsa(1)
      ypts=yptsa(1)
      crot=crota(1)
      srot=srota(1)
      xlon=xlona(1)
      xlat=xlata(1)
      ylon=ylona(1)
      ylat=ylata(1)
      area=areaa(1)

    endif

    return

  endsubroutine gdswzd_scalar

  !> Decodes the grib 2 grid definition template and returns
  !! one of the following (for 2d-arrays):
  !! - iopt=0 Grid and earth coordinates of all grid points.
  !! - iopt=+1 Earth coordinates of selected grid coordinates.
  !! - iopt=-1 Grid coordinates of selected earth coordinates.
  !!
  !! The current code recognizes the following projections,
  !! where "igdtnum" is the grid definition template number:
  !! - igdtnum=00 Equidistant Cylindrical
  !! - igdtnum=01 Rotated Equidistant Cylindrical. "E" and non-"E" staggered
  !! - igdtnum=10 Mercator Cyclindrical
  !! - igdtnum=20 Polar Stereographic Azimuthal
  !! - igdtnum=30 Lambert Conformal Conical
  !! - igdtnum=40 Gaussian Equidistant Cyclindrical
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values. Also if iopt=0,
  !! if the number of grid points exceeds the number allotted,
  !! then all the output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !!
  !! Optionally, the vector rotations, map jacobians and
  !! grid box areas may be returned.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present. To compute the map jacobians, the
  !! optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! @param[in] igdtnum Grid definition template number.
  !!
  !! @param[in] igdtmpl Grid definition template array.
  !! Corresponds to the gfld%igdtmpl component of the
  !! NCEPLIBS-g2 gridmod data structure.
  !! See igdtmpl definition in gdswzd_1d_array() for full details.
  !!
  !! @param[in] igdtlen Number of elements of the grid definition
  !! template array. Corresponds to the gfld%igdtlen
  !! component of the ncep g2 library gridmod data structure.
  !!
  !! @param[in] iopt Option flag.
  !! - 0 Earth coords of all the grid points.
  !! - 1 Earth coords of selected grid coords.
  !! - -1 Grid coords of selected earth coords
  !!
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!      Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates.
  !! @param[inout] ypts Grid y point coordinates.
  !! @param[inout] rlon Earth longitudes in degrees E.
  !! (Acceptable range: -360. to 360.)
  !!
  !! @param[inout] rlat Earth latitudes in degrees N.
  !! (Acceptable range: -90. to 90.)
  !!
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @param[out] xlon dx/dlat in 1/degrees
  !! @param[out] xlat dy/dlon in 1/degrees
  !! @param[out] ylon dy/dlon in 1/degrees
  !! @param[out] ylat dy/dlat in 1/degrees
  !! @param[out] area Area weights in m^2.
  !! Proportional to the square of the map factor in the case of
  !! conformal projections
  !!
  !! @author George Gayno, Mark Iredell
  !! @date Jan 2015
  subroutine gdswzd_2d_array(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                             xpts,ypts,rlon,rlat,nret, &
                             crot,srot,xlon,xlat,ylon,ylat,area)

    implicit none
    !
    integer,intent(in) :: igdtnum,igdtlen
    integer,intent(in) :: igdtmpl(igdtlen)
    integer,intent(in) :: iopt,npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon(:,:),rlat(:,:)
    real,intent(inout) :: xpts(:,:),ypts(:,:)
    real,optional,intent(out) :: crot(:,:),srot(:,:)
    real,optional,intent(out) :: xlon(:,:),xlat(:,:)
    real,optional,intent(out) :: ylon(:,:),ylat(:,:),area(:,:)

    call gdswzd_1d_array(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                         xpts,ypts,rlon,rlat,nret, &
                         crot,srot,xlon,xlat,ylon,ylat,area)

  endsubroutine gdswzd_2d_array

  !> Decodes the grib 2 grid definition template and returns one of the following:
  !! - iopt=0 Grid and earth coordinates of all grid points.
  !! - iopt=+1 Earth coordinates of selected grid coordinates.
  !! - iopt=-1 Grid coordinates of selected earth coordinates.
  !!
  !! The current code recognizes the following projections,
  !! where "igdtnum" is the grid definition template number:
  !! - igdtnum=00 Equidistant Cylindrical
  !! - igdtnum=01 Rotated Equidistant Cylindrical. "E" and non-"E" staggered
  !! - igdtnum=10 Mercator Cyclindrical
  !! - igdtnum=20 Polar Stereographic Azimuthal
  !! - igdtnum=30 Lambert Conformal Conical
  !! - igdtnum=40 Gaussian Equidistant Cyclindrical
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values. Also if iopt=0,
  !! if the number of grid points exceeds the number allotted,
  !! then all the output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !!
  !! Optionally, the vector rotations, map jacobians and
  !! grid box areas may be returned.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present. To compute the map jacobians, the
  !! optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! @param[in] igdtnum Grid definition template number.
  !! Corresponds to the gfld%igdtnum component of the ncep g2 library
  !! gridmod data structure:
  !! - 00 - Equidistant Cylindrical
  !! - 01 - Rotated Equidistant Cylindrical. "E" and non-"E" staggered
  !! - 10 - Mercator Cyclindrical
  !! - 20 - Polar Stereographic Azimuthal
  !! - 30 - Lambert Conformal Conical
  !! - 40 - Gaussian Equidistant Cyclindrical
  !!
  !! @param[in] igdtmpl Grid definition template array.
  !! Corresponds to the gfld%igdtmpl component of the
  !! NCEPLIBS-g2 gridmod data structure
  !!
  !! Section 3 Info:
  !!
  !! All Map Projections:
  !! - 1: Shape of earth, octet 15.
  !! - 2: Scale factor of spherical earth radius, octet 16.
  !! - 3: Scaled value of radius of spherical earth, octets 17-20.
  !! - 4: Scale factor of major axis of elliptical earth, octet 21.
  !! - 5: Scaled value of major axis of elliptical earth, octets 22-25.
  !! - 6: Scale factor of minor axis of elliptical earth, octet 26.
  !! - 7: Scaled value of minor axis of elliptical earth, octets 27-30.
  !!
  !! Equidistant Cyclindrical:
  !! - 8:  Number of points along a parallel, octs 31-34.
  !! - 9:  Number of points along a meridian, octs 35-38.
  !! - 10: Basic angle of initial production domain, octets 39-42.
  !! - 11: Subdivisions of basic angle, octets 43-46.
  !! - 12: Latitude of first grid point, octets 47-50.
  !! - 13: Longitude of first grid point, octets 51-54.
  !! - 14: Resolution and component flags, octet 55.
  !! - 15: Latitude of last grid point, octets 56-59.
  !! - 16: Longitude of last grid point, octets 60-63.
  !! - 17: i-direction increment, octets 64-67.
  !! - 18: j-direction increment, octets 68-71.
  !! - 19: Scanning mode, octet 72.
  !!
  !! Mercator Cyclindrical:
  !! - 8:  Number of points along a parallel, octs 31-34.
  !! - 9:  Number of points along a meridian, octs 35-38.
  !! - 10: Latitude of first point, octets 39-42.
  !! - 11: Longitude of first point, octets 43-46.
  !! - 12: Resolution and component flags, octet 47.
  !! - 13: Tangent latitude, octets 48-51.
  !! - 14: Latitude of last point, octets 52-55.
  !! - 15: Longitude of last point, octets 56-59.
  !! - 16: Scanning mode flags, octet 60.
  !! - 17: Orientation of grid, octets 61-64.
  !! - 18: Longitudinal grid length, octets 65-68.
  !! - 19: Latitudinal grid length, octets 69-72.
  !!
  !! Lambert Conformal Conical:
  !! - 8:  Number of points along x-axis, octs 31-34.
  !! - 9:  Number of points along y-axis, octs 35-38.
  !! - 10: Latitude of first point, octets 39-42.
  !! - 11: Longitude of first point, octets 43-46.
  !! - 12: Resolution of component flag, octet 47.
  !! - 13: Latitude where grid lengths specified,octets 48-51.
  !! - 14: Longitude of meridian that is parallel to y-axis, octets 52-55.
  !! - 15: x-direction grid length, octets 56-59.
  !! - 16: y-direction grid length, octets 60-63.
  !! - 17: Projection center flag, octet 64.
  !! - 18: Scanning mode, octet 65.
  !! - 19: First tangent latitude from pole, octets 66-69.
  !! - 20: Second tangent latitude from pole, octets 70-73.
  !! - 21: Latitude of south pole of projection, octets 74-77.
  !! - 22: Longitude of south pole of projection, octets 78-81.
  !!
  !! Gaussian Cylindrical:
  !! - 8:  Number of points along a parallel, octs 31-34.
  !! - 9:  Number of points along a meridian, octs 35-38.
  !! - 10: Basic angle of initial production domain, octets 39-42.
  !! - 11: Subdivisions of basic angle, octets 43-46.
  !! - 12: Latitude of first grid point, octets 47-50.
  !! - 13: Longitude of first grid point, octets 51-54.
  !! - 14: Resolution and component flags, octet 55.
  !! - 15: Latitude of last grid point, octets 56-59.
  !! - 16: Longitude of last grid point, octets 60-63.
  !! - 17: i-direction increment, octets 64-67.
  !! - 18: Number of parallels between pole and equator, octets 68-71.
  !! - 19: Scanning mode, octet 72.
  !!
  !! Polar Stereographic Azimuthal:
  !! - 8:  Number of points along x-axis, octets 31-34.
  !! - 9:  Number of points along y-axis, octets 35-38.
  !! - 10: Latitude of first grid point, octets 39-42.
  !! - 11: Longitude of first grid point, octets 43-46.
  !! - 12: Resolution and component flags, octet 47.
  !! - 13: True latitude, octets 48-51.
  !! - 14: Orientation longitude, octets 52-55.
  !! - 15: x-direction grid length, octets 56-59.
  !! - 16: y-direction grid length, octets 60-63.
  !! - 17: Projection center flag, octet 64.
  !! - 18: Scanning mode flags, octet 65.
  !!
  !! Rotated Equidistant Cyclindrical:
  !! - 8:  Number of points along a parallel, octs 31-34.
  !! - 9:  Number of points along a meridian, octs 35-38.
  !! - 10: Basic angle of initial production domain, octets 39-42.
  !! - 11: Subdivisions of basic angle, octets 43-46.
  !! - 12: Latitude of first grid point, octets 47-50.
  !! - 13: Longitude of first grid point, octets 51-54.
  !! - 14: Resolution and component flags, octet 55.
  !! - 15: Latitude of last grid point, octets 56-59.
  !! - 16: Longitude of last grid point, octets 60-63.
  !! - 17: i-direction increment, octets 64-67.
  !! - 18: j-direction increment, octets 68-71.
  !! - 19: Scanning mode, octet 72.
  !! - 20: Latitude of southern pole of projection, octets 73-76.
  !! - 21: Longitude of southern pole of projection, octets 77-80.
  !! - 22: Angle of rotation of projection, octs 81-84.
  !!
  !! @param[in] igdtlen Number of elements of the grid definition
  !! template array. Corresponds to the gfld%igdtlen
  !! component of the ncep g2 library gridmod data structure.
  !!
  !! @param[in] iopt Option flag.
  !! - 0 Earth coords of all the grid points.
  !! - 1 Earth coords of selected grid coords.
  !! - -1 Grid coords of selected earth coords
  !!
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!      Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates.
  !! @param[inout] ypts Grid y point coordinates.
  !! @param[inout] rlon Earth longitudes in degrees E.
  !! (Acceptable range: -360. to 360.)
  !!
  !! @param[inout] rlat Earth latitudes in degrees N.
  !! (Acceptable range: -90. to 90.)
  !!
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @param[out] xlon dx/dlat in 1/degrees
  !! @param[out] xlat dy/dlon in 1/degrees
  !! @param[out] ylon dy/dlon in 1/degrees
  !! @param[out] ylat dy/dlat in 1/degrees
  !! @param[out] area Area weights in m^2.
  !! Proportional to the square of the map factor in the case of
  !! conformal projections
  !!
  !! @author George Gayno, Mark Iredell
  !! @date Jan 2015
  subroutine gdswzd_1d_array(igdtnum,igdtmpl,igdtlen,iopt,npts,fill, &
                             xpts,ypts,rlon,rlat,nret, &
                             crot,srot,xlon,xlat,ylon,ylat,area)
    integer,intent(in) :: igdtnum,igdtlen
    integer,intent(in) :: igdtmpl(igdtlen)
    integer,intent(in) :: iopt,npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon(npts),rlat(npts)
    real,intent(inout) :: xpts(npts),ypts(npts)
    real,optional,intent(out) :: crot(npts),srot(npts)
    real,optional,intent(out) :: xlon(npts),xlat(npts)
    real,optional,intent(out) :: ylon(npts),ylat(npts),area(npts)

    type(grib2_descriptor) :: desc
    class(ip_grid),allocatable :: grid

    desc=init_descriptor(igdtnum,igdtlen,igdtmpl)
    call init_grid(grid,desc)

    call gdswzd_grid(grid,iopt,npts,fill, &
                     xpts,ypts,rlon,rlat,nret, &
                     crot,srot,xlon,xlat,ylon,ylat,area)

  endsubroutine gdswzd_1d_array

  !> Decodes the grib grid description section and
  !! returns one of the following (for 1-d arrays):
  !! - iopt=0 Grid and earth coordinates of all grid points.
  !! - iopt=+1 Earth coordinates of selected grid coordinates.
  !! - iopt=-1 Grid coordinates of selected earth coordinates.
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values.  Also if iopt=0,
  !! if the number of grid points exceeds the number allotted,
  !! then all the output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !!
  !! Optionally, the vector rotations, map jacobians and
  !! grid box areas may be returned.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present. To compute the map jacobians, the
  !! optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! The current code recognizes the following projections:
  !! - kgds(1)=000 Equidistant Cylindrical
  !! - kgds(1)=001 Mercator Cylindrical
  !! - kgds(1)=003 lambert Conformal Conical
  !! - kgds(1)=004 Gaussian Cylindrical
  !! - kgds(1)=005 Polar Stereographic azimuthal
  !! - kgds(1)=203 E-staggered Rotated Equidistant Cylindrical
  !! - kgds(1)=205 B-staggered Rotated Equidistant Cylindrical
  !!
  !! @param[in] kgds GDS parameters as decoded by w3fi63.
  !! @param[in] iopt Option flag.
  !! - 0 Earth coords of all the grid points.
  !! - 1 Earth coords of selected grid coords.
  !! - -1 Grid coords of selected earth coords
  !!
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!      Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates.
  !! @param[inout] ypts Grid y point coordinates.
  !! @param[inout] rlon Earth longitudes in degrees E.
  !! (Acceptable range: -360. to 360.)
  !!
  !! @param[inout] rlat Earth latitudes in degrees N.
  !! (Acceptable range: -90. to 90.)
  !!
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @param[out] xlon dx/dlat in 1/degrees
  !! @param[out] xlat dy/dlon in 1/degrees
  !! @param[out] ylon dy/dlon in 1/degrees
  !! @param[out] ylat dy/dlat in 1/degrees
  !! @param[out] area Area weights in m^2.
  !! Proportional to the square of the map factor in the case of
  !! conformal projections
  !!
  !! @author George Gayno, Mark Iredell
  !! @date April 1996
  subroutine gdswzd_grib1(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret, &
                          crot,srot,xlon,xlat,ylon,ylat,area)
    integer,intent(in) :: iopt,kgds(200),npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon(npts),rlat(npts)
    real,intent(inout) :: xpts(npts),ypts(npts)
    real,optional,intent(out) :: crot(npts),srot(npts)
    real,optional,intent(out) :: xlon(npts),xlat(npts)
    real,optional,intent(out) :: ylon(npts),ylat(npts),area(npts)

    type(grib1_descriptor) :: desc
    class(ip_grid),allocatable :: grid

    desc=init_descriptor(kgds)
    call init_grid(grid,desc)

    call gdswzd_grid(grid,iopt,npts,fill, &
                     xpts,ypts,rlon,rlat,nret, &
                     crot,srot,xlon,xlat,ylon,ylat,area)

  endsubroutine gdswzd_grib1

  !> Decodes the grib grid description section and returns
  !! one of the following (for 2-d arrays):
  !! - iopt=0 Grid and earth coordinates of all grid points.
  !! - iopt=+1 Earth coordinates of selected grid coordinates.
  !! - iopt=-1 Grid coordinates of selected earth coordinates.
  !!
  !! If the selected coordinates are more than one gridpoint
  !! beyond the the edges of the grid domain, then the relevant
  !! output elements are set to fill values.  Also if iopt=0,
  !! if the number of grid points exceeds the number allotted,
  !! then all the output elements are set to fill values.
  !!
  !! The actual number of valid points computed is returned too.
  !!
  !! Optionally, the vector rotations, map jacobians and
  !! grid box areas may be returned.
  !!
  !! To compute the vector rotations, the optional arguments 'srot' and 'crot'
  !! must be present. To compute the map jacobians, the
  !! optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be present.
  !!
  !! To compute the grid box areas, the optional argument
  !! 'area' must be present.
  !!
  !! The current code recognizes the following projections:
  !! - kgds(1)=000 Equidistant Cylindrical
  !! - kgds(1)=001 Mercator Cylindrical
  !! - kgds(1)=003 lambert Conformal Conical
  !! - kgds(1)=004 Gaussian Cylindrical
  !! - kgds(1)=005 Polar Stereographic azimuthal
  !! - kgds(1)=203 E-staggered Rotated Equidistant Cylindrical
  !! - kgds(1)=205 B-staggered Rotated Equidistant Cylindrical
  !!
  !! @param[in] kgds GDS parameters as decoded by w3fi63.
  !! @param[in] iopt Option flag.
  !! - 0 Earth coords of all the grid points.
  !! - 1 Earth coords of selected grid coords.
  !! - -1 Grid coords of selected earth coords
  !!
  !! @param[in] npts Maximum number of coordinates.
  !! @param[in] fill Fill value to set invalid output data.
  !!      Must be impossible value; suggested value: -9999.
  !! @param[inout] xpts Grid x point coordinates.
  !! @param[inout] ypts Grid y point coordinates.
  !! @param[inout] rlon Earth longitudes in degrees E.
  !! (Acceptable range: -360. to 360.)
  !!
  !! @param[inout] rlat Earth latitudes in degrees N.
  !! (Acceptable range: -90. to 90.)
  !!
  !! @param[out] nret Number of valid points computed.
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @param[out] xlon dx/dlat in 1/degrees
  !! @param[out] xlat dy/dlon in 1/degrees
  !! @param[out] ylon dy/dlon in 1/degrees
  !! @param[out] ylat dy/dlat in 1/degrees
  !! @param[out] area Area weights in m^2.
  !! Proportional to the square of the map factor in the case of
  !! conformal projections
  !!
  !! @author George Gayno, Mark Iredell
  !! @date April 1996
  subroutine gdswzd_2d_array_grib1(kgds,iopt,npts,fill,xpts,ypts,rlon,rlat,nret, &
                                   crot,srot,xlon,xlat,ylon,ylat,area)

    !$$$
    integer,intent(in) :: iopt,kgds(200),npts
    integer,intent(out) :: nret
    !
    real,intent(in) :: fill
    real,intent(inout) :: rlon(:,:),rlat(:,:)
    real,intent(inout) :: xpts(:,:),ypts(:,:)
    real,optional,intent(out) :: crot(:,:),srot(:,:)
    real,optional,intent(out) :: xlon(:,:),xlat(:,:)
    real,optional,intent(out) :: ylon(:,:),ylat(:,:),area(:,:)

    type(grib1_descriptor) :: desc
    class(ip_grid),allocatable :: grid

    desc=init_descriptor(kgds)
    call init_grid(grid,desc)

    call gdswzd_grid(grid,iopt,npts,fill, &
                     xpts,ypts,rlon,rlat,nret, &
                     crot,srot,xlon,xlat,ylon,ylat,area)

  endsubroutine gdswzd_2d_array_grib1

endmodule gdswzd_mod
