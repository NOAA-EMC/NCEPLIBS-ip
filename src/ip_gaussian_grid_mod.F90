!> @file
!! @brief Gaussian grid coordinate transformations.
!! @author Mark Iredell, George Gayno, Kyle Gerheiser
!! @date July 2021

!> @brief Gaussian grid coordinate transformations.
!!
!! Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.40
!! Gaussian
!! Latitude/Longitude](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-40.shtml).
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
!! @date July 2021
module ip_gaussian_grid_mod
    use ip_grid_descriptor_mod
    use ip_grid_mod
    use earth_radius_mod
    use ip_constants_mod
    use sp_mod
    implicit none

    private
    public :: ip_gaussian_grid

    type,extends(ip_grid) :: ip_gaussian_grid
        integer :: jh !< Scan mode flag in 'j' direction. When '1' points scan from N to S. When "-1" points scan from S to N.
        real :: dlon !< "i"-direction increment. GRIB2 Section 3, octets 64-67.
        real :: rlat1 !<  Latitude of first grid point. GRIB2 Section 3, octets 47-50.
        real :: rlon1 !< Longitude of first grid point. GRIB2 Section 3, octets 51-54.
        real :: rlon2 !< Longitude of last grid point. GRIB2 Section 3, octets 60-63.
        real :: hi !< Scan mode flag in 'i' direction. When '1' points scan from W to E. When "-1" points scan from E to W.
        integer :: jg !< Number of parallels between the equator and pole times 2. GRIB2 Section 3, octets 68-71.
        integer :: jscan !< Scanning mode in the 'j' direction. GRIB2 Section 3, octet 72.
    contains
        !> Initializes a gaussian grid given a grib1_descriptor object. @return N/A
        procedure :: init_grib1
        !> Initializes a gaussian grid given a grib2_descriptor object. @return N/A
        procedure :: init_grib2
        !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
        !> for Gaussian grids. @return N/A
        procedure :: gdswzd=>gdswzd_gaussian
    endtype ip_gaussian_grid

    integer :: j1 !< 'j' index of first grid point within the global array of latitudes.
    integer :: jh !< Scan mode flag in 'j' direction. When '1' points scan from N to S. When "-1" points scan from S to N.
    real,allocatable :: blat(:) !< Gaussian latitude for each parallel.
    real :: dlon !< "i"-direction increment. GRIB2 Section 3, octets 64-67.
    real :: rerth !< Radius of the earth. GRIB2 Section 3, octets 15-30.
    real,allocatable :: ylat_row(:) !< dy/dlat for each row in 1/degrees.

contains

    !> Initializes a gaussian grid given a grib1_descriptor object.
    !>
    !> @param[inout] self The grid to initialize
    !> @param[in] g1_desc A grib1_descriptor
    !>
    !> @author Kyle Gerheiser
    !> @date July 2021
    subroutine init_grib1(self,g1_desc)
        class(ip_gaussian_grid),intent(inout) :: self
        type(grib1_descriptor),intent(in) :: g1_desc

        integer :: iscan,jg

        associate(kgds=>g1_desc%gds)
            self%rerth=6.3712e6
            self%eccen_squared=0.0

            self%im=kgds(2)
            self%jm=kgds(3)
            self%rlat1=kgds(4)*1.e-3
            self%rlon1=kgds(5)*1.e-3
            self%rlon2=kgds(8)*1.e-3
            self%jg=kgds(10)*2
            iscan=mod(kgds(11)/128,2)
            self%jscan=mod(kgds(11)/64,2)
            self%hi=(-1.)**iscan
            self%jh=(-1)**self%jscan
            self%dlon=self%hi*(mod(self%hi*(self%rlon2-self%rlon1)-1+3600,360.)+1)/(self%im-1)

            self%iwrap=0
            self%jwrap1=0
            self%jwrap2=0
            self%nscan=mod(kgds(11)/32,2)
            self%nscan_field_pos=self%nscan
            self%kscan=0

            self%iwrap=nint(360/abs(self%dlon))
            if(self%im.lt.self%iwrap) self%iwrap=0

            if(self%iwrap.gt.0.and.mod(self%iwrap,2).eq.0) then
                jg=kgds(10)*2
                if(self%jm.eq.self%jg) then
                    self%jwrap1=1
                    self%jwrap2=2*self%jm+1
                endif
            endif

        endassociate
    endsubroutine init_grib1

    !> Initializes a gaussian grid given a grib2_descriptor object.
    !> @param[inout] self The grid to initialize
    !> @param[in] g2_desc A grib2_descriptor
    !>
    !> @author Kyle Gerheiser
    !> @date July 2021
    subroutine init_grib2(self,g2_desc)
        class(ip_gaussian_grid),intent(inout) :: self
        type(grib2_descriptor),intent(in) :: g2_desc

        integer :: iscale,iscan,jg

        associate(igdtmpl=>g2_desc%gdt_tmpl,igdtlen=>g2_desc%gdt_len)
            call earth_radius(igdtmpl,igdtlen,self%rerth,self%eccen_squared)

            self%im=igdtmpl(8)
            self%jm=igdtmpl(9)
            iscale=igdtmpl(10)*igdtmpl(11)
            if(iscale.eq.0) iscale=10**6
            self%rlat1=float(igdtmpl(12))/float(iscale)
            self%rlon1=float(igdtmpl(13))/float(iscale)
            self%rlon2=float(igdtmpl(16))/float(iscale)
            self%jg=igdtmpl(18)*2
            iscan=mod(igdtmpl(19)/128,2)
            self%jscan=mod(igdtmpl(19)/64,2)
            self%hi=(-1.)**iscan
            self%jh=(-1)**self%jscan
            self%dlon=self%hi*(mod(self%hi*(self%rlon2-self%rlon1)-1+3600,360.)+1)/(self%im-1)

            self%iwrap=nint(360/abs(self%dlon))
            if(self%im.lt.self%iwrap) self%iwrap=0
            self%jwrap1=0
            self%jwrap2=0
            if(self%iwrap.gt.0.and.mod(self%iwrap,2).eq.0) then
                jg=igdtmpl(18)*2
                if(self%jm.eq.jg) then
                    self%jwrap1=1
                    self%jwrap2=2*self%jm+1
                endif
            endif
            self%nscan=mod(igdtmpl(19)/32,2)
            self%nscan_field_pos=self%nscan
            self%kscan=0
        endassociate

    endsubroutine init_grib2

    !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
    !> for Gaussian grids.
    !>
    !> If the selected coordinates are more than one gridpoint
    !> beyond the the edges of the grid domain, then the relevant
    !> output elements are set to fill values.
    !>
    !> The actual number of valid points computed is returned too.
    !> Optionally, the vector rotations, the map jacobians and
    !> the grid box areas may be returned as well.
    !>
    !> To compute the vector rotations, the optional arguments 'srot' and 'crot'
    !> must be present.
    !>
    !> To compute the map jacobians, the optional arguments
    !> 'xlon', 'xlat', 'ylon', 'ylat' must be present.
    !>
    !> To compute the grid box areas, the optional argument
    !> 'area' must be present.
    !>
    !> @param[in] self The grid object gdswzd was called on.
    !> @param[in] iopt option flag
    !>            - +1 to compute earth coords of selected grid coords.
    !>            - -1 o compute grid coords of selected earth coords.
    !> @param[in] npts Maximum number of coordinates.
    !> @param[in] fill Fill value to set invalid output data.
    !>            Must be impossible value; suggested value: -9999.
    !> @param[inout] xpts Grid x point coordinates if iopt>0.
    !> @param[inout] ypts Grid y point coordinates if iopt>0.
    !> @param[inout] rlon Earth longitudes in degrees e if iopt<0
    !>                   (Acceptable range: -360. to 360.)
    !> @param[inout] rlat Earth latitudes in degrees n if iopt<0
    !>                (Acceptable range: -90. to 90.)
    !> @param[out] nret Number of valid points computed.
    !> @param[out] crot Optional clockwise vector rotation cosines.
    !> @param[out] srot Optional clockwise vector rotation sines.
    !> @param[out] xlon Optional dx/dlon in 1/degrees.
    !> @param[out] xlat Optional dx/dlat in 1/degrees.
    !> @param[out] ylon Optional dy/dlon in 1/degrees.
    !> @param[out] ylat Optional dy/dlat in 1/degrees.
    !> @param[out] area Optional area weights in m**2.
    !>
    !> @author Mark Iredell, George Gayno, Kyle Gerheiser
    !> @date July 2021
    subroutine gdswzd_gaussian(self,iopt,npts,fill, &
                               xpts,ypts,rlon,rlat,nret, &
                               crot,srot,xlon,xlat,ylon,ylat,area)
        implicit none
        !
        class(ip_gaussian_grid),intent(in) :: self
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
        integer                        :: jscan,im,jm
        integer                        :: j,ja,jg
        integer                        :: n
        !
        logical                        :: lrot,lmap,larea
        !
        real,allocatable   :: alat(:),alat_jscan(:)
        real,allocatable   :: alat_temp(:),blat_temp(:)
        real                           :: hi,rlata,rlatb,rlat1,rlon1,rlon2
        real                           :: xmax,xmin,ymax,ymin,yptsa,yptsb
        real                           :: wb

        if(present(crot)) crot=fill
        if(present(srot)) srot=fill
        if(present(xlon)) xlon=fill
        if(present(xlat)) xlat=fill
        if(present(ylon)) ylon=fill
        if(present(ylat)) ylat=fill
        if(present(area)) area=fill

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

        im=self%im
        jm=self%jm

        rlat1=self%rlat1
        rlon1=self%rlon1
        rlon2=self%rlon2

        jg=self%jg
        jscan=self%jscan
        hi=self%hi

        jh=self%jh
        dlon=self%dlon
        rerth=self%rerth

        allocate(alat_temp(jg))
        allocate(blat_temp(jg))
        call splat(4,jg,alat_temp,blat_temp)
        allocate(alat(0:jg+1))
        allocate(blat(0:jg+1))
        !$omp parallel do private(ja) schedule(static)
        do ja=1,jg
            alat(ja)=real(dpr*asin(alat_temp(ja)))
            blat(ja)=blat_temp(ja)
        enddo
        !$omp end parallel do
        deallocate(alat_temp,blat_temp)
        alat(0)=180.-alat(1)
        alat(jg+1)=-alat(0)
        blat(0)=-blat(1)
        blat(jg+1)=blat(0)
        j1=1
        do while(j1.lt.jg.and.rlat1.lt.(alat(j1)+alat(j1+1))/2)
            j1=j1+1
        enddo
        if(lmap) then
            allocate(alat_jscan(jg))
            do ja=1,jg
                alat_jscan(j1+jh*(ja-1))=alat(ja)
            enddo
            allocate(ylat_row(0:jg+1))
            do ja=2,(jg-1)
                ylat_row(ja)=2.0/(alat_jscan(ja+1)-alat_jscan(ja-1))
            enddo
            ylat_row(1)=1.0/(alat_jscan(2)-alat_jscan(1))
            ylat_row(0)=ylat_row(1)
            ylat_row(jg)=1.0/(alat_jscan(jg)-alat_jscan(jg-1))
            ylat_row(jg+1)=ylat_row(jg)
            deallocate(alat_jscan)
        endif
        xmin=0
        xmax=im+1
        if(im.eq.nint(360/abs(dlon))) xmax=im+2
        ymin=0.5
        ymax=jm+0.5
        nret=0

        !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
        if(iopt.eq.0.or.iopt.eq.1) then
            !$omp parallel do private(n,j,wb,rlata,rlatb) reduction(+:nret) schedule(static)
            do n=1,npts
                if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
                   ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                    rlon(n)=mod(rlon1+dlon*(xpts(n)-1)+3600,360.)
                    j=int(ypts(n))
                    wb=ypts(n)-j
                    rlata=alat(j1+jh*(j-1))
                    rlatb=alat(j1+jh*j)
                    rlat(n)=rlata+wb*(rlatb-rlata)
                    nret=nret+1
                    if(lrot) call gaussian_vect_rot(crot(n),srot(n))
                    if(lmap) call gaussian_map_jacob(ypts(n), &
                                                     xlon(n),xlat(n),ylon(n),ylat(n))
                    if(larea) call gaussian_grid_area(ypts(n),area(n))
                else
                    rlon(n)=fill
                    rlat(n)=fill
                endif
            enddo
            !$omp end parallel do

            !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
        elseif(iopt.eq.-1) then
            !$omp parallel do private(n,ja,yptsa,yptsb,wb) reduction(+:nret) schedule(static)
            do n=1,npts
                xpts(n)=fill
                ypts(n)=fill
                if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90) then
                    xpts(n)=1+hi*mod(hi*(rlon(n)-rlon1)+3600,360.)/dlon
                    ja=min(int((jg+1)/180.*(90-rlat(n))),jg)
                    if(rlat(n).gt.alat(ja)) ja=max(ja-2,0)
                    if(rlat(n).lt.alat(ja+1)) ja=min(ja+2,jg)
                    if(rlat(n).gt.alat(ja)) ja=ja-1
                    if(rlat(n).lt.alat(ja+1)) ja=ja+1
                    yptsa=1+jh*(ja-j1)
                    yptsb=1+jh*(ja+1-j1)
                    wb=(alat(ja)-rlat(n))/(alat(ja)-alat(ja+1))
                    ypts(n)=yptsa+wb*(yptsb-yptsa)
                    if(xpts(n).ge.xmin.and.xpts(n).le.xmax.and. &
                       ypts(n).ge.ymin.and.ypts(n).le.ymax) then
                        nret=nret+1
                        if(lrot) call gaussian_vect_rot(crot(n),srot(n))
                        if(lmap) call gaussian_map_jacob(ypts(n), &
                                                         xlon(n),xlat(n),ylon(n),ylat(n))
                        if(larea) call gaussian_grid_area(ypts(n),area(n))
                    else
                        xpts(n)=fill
                        ypts(n)=fill
                    endif
                endif
            enddo
            !$omp end parallel do
        endif
        deallocate(alat,blat)
        if(allocated(ylat_row)) deallocate(ylat_row)

    endsubroutine gdswzd_gaussian

    !> Computes the vector rotation sines and cosines for a gaussian
    !> cylindrical grid.
    !>
    !> @param[out] crot Clockwise vector rotation cosines.
    !> @param[out] srot Clockwise vector rotation sines.
    !>
    !> @note
    !> ugrid=crot*uearth-srot*vearth;
    !> vgrid=srot*uearth+crot*vearth)
    !>
    !> @author George Gayno
    !> @date July 2021
    subroutine gaussian_vect_rot(crot,srot)
        implicit none

        real,intent(out) :: crot,srot

        crot=1.0
        srot=0.0

    endsubroutine gaussian_vect_rot

    !> Computes the map jacobians for a gaussian cylindrical grid.
    !>
    !> @param[in] ypts y-index of grid point.
    !> @param[out] xlon dx/dlon in 1/degrees.
    !> @param[out] xlat dx/dlat in 1/degrees.
    !> @param[out] ylon dy/dlon in 1/degrees.
    !> @param[out] ylat dy/dlat in 1/degrees.
    !>
    !> @author George Gayno
    !> @date July 2021
    subroutine gaussian_map_jacob(ypts,xlon,xlat,ylon,ylat)
        implicit none

        real,intent(in) :: ypts
        real,intent(out) :: xlon,xlat,ylon,ylat

        xlon=1/dlon
        xlat=0.
        ylon=0.
        ylat=ylat_row(nint(ypts))

    endsubroutine gaussian_map_jacob

    !> Computes the grid box area for a gaussian cylindrical grid.
    !>
    !> @param[in] ypts y-index of grid point.
    !> @param[out] area Area weights in m^2
    !>
    !> @author Mark Iredell, George Gayno
    !> @date July 2021
    subroutine gaussian_grid_area(ypts,area)
        implicit none

        real,intent(in) :: ypts
        real,intent(out) :: area

        integer                        :: j

        real                           :: wb,wlat,wlata,wlatb

        j=int(ypts)
        wb=ypts-j
        wlata=blat(j1+jh*(j-1))
        wlatb=blat(j1+jh*j)
        wlat=wlata+wb*(wlatb-wlata)
        area=real(rerth**2*wlat*dlon/dpr)

    endsubroutine gaussian_grid_area
endmodule ip_gaussian_grid_mod

