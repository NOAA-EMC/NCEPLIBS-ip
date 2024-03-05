!> @file
!> @brief GDS wizard for mercator cylindrical.
!>
!> @author Iredell @date 96-04-10

!> @brief GDS wizard for mercator cylindrical.
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.10 -
!> Mercator](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-10.shtml).
!>
!> @author Iredell @date 96-04-10
module ip_mercator_grid_mod
    use ip_grid_descriptor_mod
    use ip_grid_mod
    use ip_constants_mod,only:dpr,pi
    use earth_radius_mod
    implicit none

    private
    public :: ip_mercator_grid

    type,extends(ip_grid) :: ip_mercator_grid
        real :: rlat1 !< Latitude of first grid point. Section 3, octets 39-42.
        real :: rlon1 !< Longitude of first grid point. Section 3, octets 43-46.
        real :: rlon2 !< Longitude of last grid point. Section 3, octets 56-59.
        real :: rlati !< Latitude at which the Mercator projection intersects the Earth. Section 3, octets 48-51.
        real :: hi !< Scan mode in the 'i' direction. Section 3, octet 60.
        real :: dlon !< Longitudinal direction grid length. Section 3, octets 65-68.
        real :: dphi !< Latitudinal direction grid length. Section 3, octets 69-72.
    contains
        !> Initializes a gaussian grid given a grib1_descriptor object. @return N/A
        procedure :: init_grib1
        !> Initializes a gaussian grid given a grib2_descriptor object. @return N/A
        procedure :: init_grib2
        !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
        !> for Gaussian grids. @return N/A
        procedure :: gdswzd=>gdswzd_mercator !< gdswzd() @return N/A
    endtype ip_mercator_grid

    real :: dlon !< Longitudinal direction grid length.
    real :: dphi !< Latitudinal direction grid length.
    real :: rerth !< Radius of the Earth.

contains

    !> Initializes a mercator grid given a grib1_descriptor object.
    !>
    !> @param[inout] self ip_mercator_grid object.
    !> @param[in] g1_desc GRIB1 descriptor.
    !>
    !> @author Iredell @date 96-04-10
    subroutine init_grib1(self,g1_desc)
        class(ip_mercator_grid),intent(inout) :: self
        type(grib1_descriptor),intent(in) :: g1_desc

        integer :: iscan,jscan
        real :: dy,hj

        associate(kgds=>g1_desc%gds)
            self%rerth=6.3712e6
            self%eccen_squared=0.0

            self%im=kgds(2)
            self%jm=kgds(3)

            self%rlat1=kgds(4)*1.e-3
            self%rlon1=kgds(5)*1.e-3
            self%rlon2=kgds(8)*1.e-3
            self%rlati=kgds(9)*1.e-3

            iscan=mod(kgds(11)/128,2)
            jscan=mod(kgds(11)/64,2)

            dy=kgds(13)
            self%hi=(-1.)**iscan
            hj=(-1.)**(1-jscan)
            self%dlon=self%hi*(mod(self%hi*(self%rlon2-self%rlon1)-1+3600,360.)+1)/(self%im-1)
            self%dphi=hj*dy/(self%rerth*cos(self%rlati/dpr))

            ! defaults
            self%iwrap=0
            self%jwrap1=0
            self%jwrap2=0
            self%nscan=mod(kgds(11)/32,2)
            self%nscan_field_pos=self%nscan
            self%kscan=0

            self%iwrap=nint(360/abs(self%dlon))
            if(self%im.lt.self%iwrap) self%iwrap=0
        endassociate

    endsubroutine init_grib1

    !> Init GRIB2.
    !>
    !> @param[inout] self ip_mercator_grid object.
    !> @param[in] g2_desc GRIB2 descriptor.
    !>
    !> @author Iredell @date 96-04-10
    subroutine init_grib2(self,g2_desc)
        class(ip_mercator_grid),intent(inout) :: self
        type(grib2_descriptor),intent(in) :: g2_desc

        integer :: iscan,jscan
        real :: hj,dy

        associate(igdtmpl=>g2_desc%gdt_tmpl,igdtlen=>g2_desc%gdt_len)

            call earth_radius(igdtmpl,igdtlen,self%rerth,self%eccen_squared)

            self%im=igdtmpl(8)
            self%jm=igdtmpl(9)

            self%rlat1=float(igdtmpl(10))*1.0e-6
            self%rlon1=float(igdtmpl(11))*1.0e-6
            self%rlon2=float(igdtmpl(15))*1.0e-6
            self%rlati=float(igdtmpl(13))*1.0e-6

            iscan=mod(igdtmpl(16)/128,2)
            jscan=mod(igdtmpl(16)/64,2)

            dy=float(igdtmpl(19))*1.0e-3
            self%hi=(-1.)**iscan
            hj=(-1.)**(1-jscan)
            self%dlon=self%hi*(mod(self%hi*(self%rlon2-self%rlon1)-1+3600,360.)+1)/(self%im-1)
            self%dphi=hj*dy/(self%rerth*cos(self%rlati/dpr))

            self%jwrap1=0
            self%jwrap2=0
            self%kscan=0
            self%nscan=mod(igdtmpl(16)/32,2)
            self%nscan_field_pos=self%nscan

            self%iwrap=nint(360/abs(self%dlon))
            if(self%im.lt.self%iwrap) self%iwrap=0

        endassociate
    endsubroutine init_grib2

    !> GDS wizard for mercator cylindrical.
    !>
    !> This routine decodes the grib 2 grid definition template (passed
    !> in integer form as decoded by the ncep g2 library) and returns
    !> one of the following:
    !> - (iopt=+1) earth coordinates of selected grid coordinates
    !> - (iopt=-1) grid coordinates of selected earth coordinates
    !>
    !> Works for mercator cylindrical projections.
    !>
    !> If the selected coordinates are more than one gridpoint beyond
    !> the the edges of the grid domain, then the relevant output
    !> elements are set to fill values.
    !>
    !> The actual number of valid points computed is returned too.
    !>
    !> Optionally, the vector rotations, map jacobians and the grid box
    !> areas may be returned. To compute the vector rotations, the
    !> optional arguments 'srot' and 'crot' must be present. To compute
    !> the map jacobians, the optional arguments 'xlon', 'xlat', 'ylon',
    !> 'ylat' must be present. to compute the grid box areas, the
    !> optional argument 'area' must be present.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | iredell | Initial
    !> 96-10-01 | iredell | protected against unresolvable points
    !> 97-10-20 | iredell | include map options
    !> 2015-01-21 | gayno | merger of gdswiz01() and gdswzd01(). Make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. Make part of a module. move vector rotation, map jacobian and grid box area computations to separate subroutines.
    !> 2015-07-13 | gayno | convert to grib 2. replace grib 1 kgds array with grib 2 grid definition template array. Rename.
    !> 2018-07-20 | wesley | add threads.
    !>
    !> @param[in] self grid descriptor.
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
    subroutine gdswzd_mercator(self,iopt,npts,fill, &
                               xpts,ypts,rlon,rlat,nret, &
                               crot,srot,xlon,xlat,ylon,ylat,area)
        implicit none
        !
        class(ip_mercator_grid),intent(in) :: self
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
        integer                          :: im,jm,n
        !
        logical                          :: lrot,lmap,larea
        !
        real                             :: hi
        real                             :: rlat1,rlon1,rlon2,rlati
        real                             :: xmax,xmin,ymax,ymin
        real                             :: ye
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if(present(crot)) crot=fill
        if(present(srot)) srot=fill
        if(present(xlon)) xlon=fill
        if(present(xlat)) xlat=fill
        if(present(ylon)) ylon=fill
        if(present(ylat)) ylat=fill
        if(present(area)) area=fill
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        im=self%im
        jm=self%jm

        rlat1=self%rlat1
        rlon1=self%rlon1
        rlon2=self%rlon2
        rlati=self%rlati

        hi=self%hi

        dlon=self%dlon
        dphi=self%dphi
        rerth=self%rerth

        ye=1-log(tan((rlat1+90)/2/dpr))/dphi
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
                    rlat(n)=2*atan(exp(dphi*(ypts(n)-ye)))*dpr-90
                    nret=nret+1
                    if(lrot) call mercator_vect_rot(crot(n),srot(n))
                    if(lmap) call mercator_map_jacob(rlat(n),xlon(n),xlat(n),ylon(n),ylat(n))
                    if(larea) call mercator_grid_area(rlat(n),area(n))
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
                if(abs(rlon(n)).le.360.and.abs(rlat(n)).lt.90) then
                    xpts(n)=1+hi*mod(hi*(rlon(n)-rlon1)+3600,360.)/dlon
                    ypts(n)=ye+log(tan((rlat(n)+90)/2/dpr))/dphi
                    if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
                       ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                        nret=nret+1
                        if(lrot) call mercator_vect_rot(crot(n),srot(n))
                        if(lmap) call mercator_map_jacob(rlat(n),xlon(n),xlat(n),ylon(n),ylat(n))
                        if(larea) call mercator_grid_area(rlat(n),area(n))
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
    endsubroutine gdswzd_mercator

    !> Vector rotation fields for mercator cylindrical grids.
    !>
    !> This subprogram computes the vector rotation sines and cosines
    !> for a mercator cylindrical grid.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "mercator_vect_rot".
    !>
    !> @param[in] crot clockwise vector rotation cosines (real)
    !> @param[in] srot clockwise vector rotation sines (real)
    !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine mercator_vect_rot(crot,srot)
        implicit none

        real,intent(out) :: crot,srot

        crot=1.0
        srot=0.0

    endsubroutine mercator_vect_rot

    !> Map jacobians for mercator cylindrical grids.
    !>
    !> This subprogram computes the map jacobians for a mercator
    !> cylindrical grid.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "mercator_map_jacob"
    !>
    !> @param[in] rlat latitude in degrees (real)
    !> @param[out] xlon dx/dlon in 1/degrees (real)
    !> @param[out] xlat dx/dlat in 1/degrees (real)
    !> @param[out] ylon dy/dlon in 1/degrees (real)
    !> @param[out] ylat dy/dlat in 1/degrees (real)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine mercator_map_jacob(rlat,xlon,xlat,ylon,ylat)
        implicit none

        real,intent(in) :: rlat
        real,intent(out) :: xlon,xlat,ylon,ylat

        xlon=1./dlon
        xlat=0.
        ylon=0.
        ylat=1./dphi/cos(rlat/dpr)/dpr

    endsubroutine mercator_map_jacob

    !> Grid box area for mercator cylindrical grids.
    !>
    !> This subprogram computes the grid box area for a mercator
    !> cylindrical grid.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-01-21 | gayno | initial version
    !> 2015-09-17 | gayno | rename as "mercator_grid_area"
    !>
    !> @param[in] rlat latitude of grid point in degrees (real)
    !> @param[out] area area weights in m**2 (real)
    !>
    !> @author Gayno @date 2015-01-21
    subroutine mercator_grid_area(rlat,area)
        implicit none

        real,intent(in) :: rlat
        real,intent(out) :: area

        area=rerth**2*cos(rlat/dpr)**2*dphi*dlon/dpr

    endsubroutine mercator_grid_area

endmodule ip_mercator_grid_mod

