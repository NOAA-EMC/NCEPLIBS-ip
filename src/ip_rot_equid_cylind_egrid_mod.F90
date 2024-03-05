!> @file
!> @brief Rotated equidistant cylindrical GRIB decoder and grid
!> coordinate transformations for Arakawa grid E.
!>
!> @author Mark Iredell, George Gayno, Kyle Gerheiser
!> @date July 2021

!> Rotated equidistant cylindrical GRIB decoder and grid coordinate
!> transformations for Arakawa grid E. (To handle the A through D
!> grids, see ip_rot_equid_cylind_grid_mod).
!>
!> The E stagger is a bit odd because the 'wind' points shift by
!> half a grid box in each row. That makes the logic tricky. So the
!> routine does its computations by rotating the grid by 45 degrees.
!>
!> See more info about [Awakawa
!> grids](https://en.wikipedia.org/wiki/Arakawa_grids).
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.1 Rotate
!> Latitude/Longitude](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-1.shtml).
!>
!> @author George Gayno, Mark Iredell, Kyle Gerheiser
!> @date July 2021
module ip_rot_equid_cylind_egrid_mod
    use iso_fortran_env,only:real64
    use ip_grid_descriptor_mod
    use ip_grid_mod
    use ip_constants_mod,only:dpr,pi
    use earth_radius_mod
    implicit none

    private
    public :: ip_rot_equid_cylind_egrid

    integer,parameter :: kd=real64 !< Kind of reals.

    type,extends(ip_grid) :: ip_rot_equid_cylind_egrid
        real(kd) :: rlon0 !< Longitude of southern pole of projection.
        real(kd) :: rlon1 !< Longitude of first grid point.
        real(kd) :: rlat1 !< Latitude of first grid point.
        real(kd) :: clat0 !< Cosine of the latitude of the southern pole of projection.
        real(kd) :: slat0 !< Sine of the latitude of the southern pole of projection.
        real(kd) :: dlats !< 'J'-direction grid increment.
        real(kd) :: dlons !< 'I'-direction grid increment.
        real(kd) :: hi !< Scan mode in the 'i' direction.
        !> Rotation flag. When '0' the u/v vector components are relative
        !> to north/east. When '1' the u/v vector components are grid
        !> relative.
        integer :: irot
    contains
        !> Initializes a rotated equidistant cylindrical grid given a
        !> grib1_descriptor object. @return N/A
        procedure :: init_grib1
        !> Initializes a rotated equidistant cylindrical grid given a
        !> grib2_descriptor object. @return N/A
        procedure :: init_grib2
        !> Calculates Earth coordinates (iopt = 1) or grid coorindates
        !> (iopt = -1). @return N/A
        procedure :: gdswzd=>gdswzd_rot_equid_cylind_egrid
    endtype ip_rot_equid_cylind_egrid

    integer :: irot !< Local copy of irot.

    real(KIND=kd) :: clat !< Cosine of the latitude.
    real(KIND=kd) :: clat0 !< Local copy of clat0.
    real(KIND=kd) :: clatr !< Cosine of the rotated latitude.
    real(KIND=kd) :: clon !< Cosine of the difference between rlon and rlon0.
    real(KIND=kd) :: dlats !< Local copy of dlats.
    real(KIND=kd) :: dlons !< Local copy of dlons.
    real(KIND=kd) :: rerth !< Radius of the Earth.
    real(KIND=kd) :: rlon0 !< Local copy of rlon0.
    real(KIND=kd) :: slat !< Sine of the latitude.
    real(KIND=kd) :: slat0 !< Local copy of slat0.
    real(KIND=kd) :: slatr !< Sine of the rotated latitude.

contains

    !> Initializes a rotated equidistant cylindrical grid given a
    !> grib1_descriptor object.
    !>
    !> @param[inout] self The grid to initialize
    !> @param[in] g1_desc A grib1_descriptor
    !>
    !> @author Kyle Gerheiser
    !> @date July 2021
    subroutine init_grib1(self,g1_desc)
        class(ip_rot_equid_cylind_egrid),intent(inout) :: self
        type(grib1_descriptor),intent(in) :: g1_desc

        integer :: iscan
        real(kd) :: rlat0

        real(kd) :: rlat1,rlon1,rlon0,slat1,clat1,slat0,clat0,clon1
        real(kd) :: slatr,clatr,clonr,rlatr,rlonr,dlats,dlons,hs,hi
        integer :: im,jm

        integer :: is1,kscan,irot

        associate(kgds=>g1_desc%gds)
            self%rerth=6.3712e6_kd
            self%eccen_squared=0.0

            im=kgds(2)
            jm=kgds(3)

            self%nscan_field_pos=3
            self%nscan=mod(kgds(11)/32,2)

            rlat1=kgds(4)*1.e-3_kd
            rlon1=kgds(5)*1.e-3_kd
            rlat0=kgds(7)*1.e-3_kd
            rlon0=kgds(8)*1.e-3_kd

            irot=mod(kgds(6)/8,2)
            kscan=mod(kgds(11)/256,2)
            iscan=mod(kgds(11)/128,2)
            hi=(-1.)**iscan
            slat1=sin(rlat1/dpr)
            clat1=cos(rlat1/dpr)
            slat0=sin(rlat0/dpr)
            clat0=cos(rlat0/dpr)
            hs=sign(1._kd,mod(rlon1-rlon0+180+3600,360._kd)-180)
            clon1=cos((rlon1-rlon0)/dpr)
            slatr=clat0*slat1-slat0*clat1*clon1
            clatr=sqrt(1-slatr**2)
            clonr=(clat0*clat1*clon1+slat0*slat1)/clatr
            rlatr=dpr*asin(slatr)
            rlonr=hs*dpr*acos(clonr)
            dlats=rlatr/(-(jm-1)/2)
            dlons=rlonr/(-((im*2-1)-1)/2)

            if(kscan.eq.0) then
                is1=(jm+1)/2
            else
                is1=jm/2
            endif

            self%im=im
            self%jm=jm
            self%rlon0=rlon0
            self%rlon1=rlon1
            self%rlat1=rlat1
            self%clat0=clat0
            self%slat0=slat0
            self%dlats=dlats
            self%dlons=dlons
            self%hi=hi
            self%irot=irot
            self%kscan=kscan

        endassociate

    endsubroutine init_grib1

    !> Initializes a rotated equidistant cylindrical grid given a grib2_descriptor object.
    !> @param[inout] self The grid to initialize
    !> @param[in] g2_desc A grib2_descriptor
    !>
    !> @author Kyle Gerheiser
    !> @date July 2021
    subroutine init_grib2(self,g2_desc)
        class(ip_rot_equid_cylind_egrid),intent(inout) :: self
        type(grib2_descriptor),intent(in) :: g2_desc

        integer :: iscale,iscan
        real(kd) :: rlat0
        integer :: i_offset_odd!, i_offset_even

        associate(igdtmpl=>g2_desc%gdt_tmpl,igdtlen=>g2_desc%gdt_len)
            call earth_radius(igdtmpl,igdtlen,self%rerth,self%eccen_squared)

            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! ROUTINE ONLY WORKS FOR "E"-STAGGER GRIDS.
            !   "V" GRID WHEN BIT 5 IS '1' AND BIT 6 IS '0'.
            !   "H" GRID WHEN BIT 5 IS '0' AND BIT 6 IS '1'.
            ! I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
            ! I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
            ! IF(I_OFFSET_ODD==I_OFFSET_EVEN) THEN
            !    CALL ROT_EQUID_CYLIND_EGRID_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
            !    RETURN
            ! ENDIF
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            self%im=igdtmpl(8)
            self%jm=igdtmpl(9)

            self%nscan=mod(igdtmpl(16)/32,2)
            self%nscan_field_pos=3

            iscale=igdtmpl(10)*igdtmpl(11)
            if(iscale.eq.0) iscale=10**6

            self%rlon0=float(igdtmpl(21))/float(iscale)
            self%dlats=float(igdtmpl(18))/float(iscale)
            ! THE GRIB2 CONVENTION FOR "I" RESOLUTION IS TWICE WHAT THIS ROUTINE ASSUMES.
            self%dlons=float(igdtmpl(17))/float(iscale)*0.5_kd

            self%irot=mod(igdtmpl(14)/8,2)

            i_offset_odd=mod(igdtmpl(19)/8,2)
            self%kscan=i_offset_odd
            iscan=mod(igdtmpl(19)/128,2)

            self%hi=(-1.)**iscan

            rlat0=float(igdtmpl(20))/float(iscale)
            rlat0=rlat0+90.0_kd

            self%slat0=sin(rlat0/dpr)
            self%clat0=cos(rlat0/dpr)

            self%rlat1=float(igdtmpl(12))/float(iscale)
            self%rlon1=float(igdtmpl(13))/float(iscale)
        endassociate
    endsubroutine init_grib2

    !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
    !> for rotated equidistant cylindrical grids.
    !>
    !> Works for e-staggered rotated equidistant cylindrical projections.
    !> The scan mode determines whether this is an "h" or "v" grid.
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
    !> @date Jan 2015
    subroutine gdswzd_rot_equid_cylind_egrid(self,iopt,npts, &
                                             fill,xpts,ypts,rlon,rlat,nret, &
                                             crot,srot,xlon,xlat,ylon,ylat,area)
        implicit none
        !
        class(ip_rot_equid_cylind_egrid),intent(in) :: self

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
        integer                        :: im,jm,is1,n
        integer                        :: kscan
        !    INTEGER                        :: I_OFFSET_ODD, I_OFFSET_EVEN
        !
        logical                        :: lrot,lmap,larea
        !
        real(KIND=kd)                  :: rlat1,rlon1
        real(KIND=kd)                  :: clonr
        real(KIND=kd)                  :: rlatr,rlonr,sbd,wbd,hs,hi
        real                           :: xmax,xmin,ymax,ymin,xptf,yptf
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if(present(crot)) crot=fill
        if(present(srot)) srot=fill
        if(present(xlon)) xlon=fill
        if(present(xlat)) xlat=fill
        if(present(ylon)) ylon=fill
        if(present(ylat)) ylat=fill
        if(present(area)) area=fill

        rlon0=self%rlon0
        irot=self%irot
        im=self%im*2-1
        jm=self%jm
        dlats=self%dlats
        dlons=self%dlons
        kscan=self%kscan
        hi=self%hi
        slat0=self%slat0
        clat0=self%clat0
        rlat1=self%rlat1
        rlon1=self%rlon1

        rerth=self%rerth

        ! IS THE EARTH RADIUS DEFINED?
        if(rerth.lt.0.) then
            call rot_equid_cylind_egrid_error(iopt,fill,rlat,rlon,xpts,ypts,npts)
            return
        endif

        sbd=rlat1
        wbd=rlon1

        if(wbd.gt.180.0) wbd=wbd-360.0
        if(kscan.eq.0) then
            is1=(jm+1)/2
        else
            is1=jm/2
        endif

        xmin=0
        xmax=im+2
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
            do n=1,npts
                xptf=ypts(n)+(xpts(n)-is1)
                yptf=ypts(n)-(xpts(n)-is1)+kscan
                if(xptf.ge.xmin.and.xptf.le.xmax.and. &
                   yptf.ge.ymin.and.yptf.le.ymax) then
                    hs=hi*sign(1.,xptf-(im+1)/2)
                    select type(desc=>self%descriptor)
                    type is(grib1_descriptor)
                        rlonr=(xptf-(im+1)/2)*dlons
                        rlatr=(yptf-(jm+1)/2)*dlats
                    type is(grib2_descriptor)
                        rlonr=(xptf-1.0_kd)*dlons+wbd
                        rlatr=(yptf-1.0_kd)*dlats+sbd
                    endselect
                    clonr=cos(rlonr/dpr)
                    slatr=sin(rlatr/dpr)
                    clatr=cos(rlatr/dpr)
                    slat=clat0*slatr+slat0*clatr*clonr
                    if(slat.le.-1) then
                        clat=0.
                        clon=cos(rlon0/dpr)
                        rlon(n)=0
                        rlat(n)=-90
                    elseif(slat.ge.1) then
                        clat=0.
                        clon=cos(rlon0/dpr)
                        rlon(n)=0
                        rlat(n)=90
                    else
                        clat=sqrt(1-slat**2)
                        clon=(clat0*clatr*clonr-slat0*slatr)/clat
                        clon=min(max(clon,-1._kd),1._kd)
                        rlon(n)=real(mod(rlon0+hs*dpr*acos(clon)+3600,360._kd))
                        rlat(n)=real(dpr*asin(slat))
                    endif
                    nret=nret+1
                    if(lrot) call rot_equid_cylind_egrid_vect_rot(rlon(n),crot(n),srot(n))
                    if(lmap) call rot_equid_cylind_egrid_map_jacob(fill,rlon(n), &
                                                                   xlon(n),xlat(n),ylon(n),ylat(n))
                    if(larea) call rot_equid_cylind_egrid_grid_area(fill,area(n))
                else
                    rlon(n)=fill
                    rlat(n)=fill
                endif
            enddo
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
        elseif(iopt.eq.-1) then
            do n=1,npts
                if(abs(rlon(n)).le.360.and.abs(rlat(n)).le.90) then
                    hs=sign(1._kd,mod(rlon(n)-rlon0+180+3600,360._kd)-180)
                    clon=cos((rlon(n)-rlon0)/dpr)
                    slat=sin(rlat(n)/dpr)
                    clat=cos(rlat(n)/dpr)
                    slatr=clat0*slat-slat0*clat*clon
                    if(slatr.le.-1) then
                        clatr=0.
                        rlonr=0
                        rlatr=-90
                    elseif(slatr.ge.1) then
                        clatr=0.
                        rlonr=0
                        rlatr=90
                    else
                        clatr=sqrt(1-slatr**2)
                        clonr=(clat0*clat*clon+slat0*slat)/clatr
                        clonr=min(max(clonr,-1._kd),1._kd)
                        rlonr=hs*dpr*acos(clonr)
                        rlatr=dpr*asin(slatr)
                    endif
                    select type(desc=>self%descriptor)
                    type is(grib1_descriptor)
                        xptf=real(((rlonr-wbd)/dlons)+1.0_kd)
                        yptf=real(((rlatr-sbd)/dlats)+1.0_kd)
                    type is(grib2_descriptor)
                        xptf=real((im+1)/2+rlonr/dlons)
                        yptf=real((jm+1)/2+rlatr/dlats)
                    endselect

                    if(xptf.ge.xmin.and.xptf.le.xmax.and. &
                       yptf.ge.ymin.and.yptf.le.ymax) then
                        xpts(n)=is1+(xptf-(yptf-kscan))/2
                        ypts(n)=(xptf+(yptf-kscan))/2
                        nret=nret+1
                        if(lrot) call rot_equid_cylind_egrid_vect_rot(rlon(n),crot(n),srot(n))
                        if(lmap) call rot_equid_cylind_egrid_map_jacob(fill,rlon(n), &
                                                                       xlon(n),xlat(n),ylon(n),ylat(n))
                        if(larea) call rot_equid_cylind_egrid_grid_area(fill,area(n))
                    else
                        xpts(n)=fill
                        ypts(n)=fill
                    endif
                else
                    xpts(n)=fill
                    ypts(n)=fill
                endif
            enddo
        endif
    endsubroutine gdswzd_rot_equid_cylind_egrid

    !> Error handler.
    !>
    !> UPON AN ERROR, THIS SUBPROGRAM ASSIGNS A "FILL" VALUE TO THE
    !> OUTPUT FIELDS.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 2015-07-13 | GAYNO | Initial version
    !> 2015-09-17 | GAYNO | Rename as "rot_equid_cylind_egrid_error"
    !>
    !> @param[in] iopt option flag
    !> - 1 to compute earth coords of selected grid coords
    !> - -1 to compute grid coords of selected earth coords
    !> @param[in] fill fill value to set invalid output data (must be
    !> impossible value; suggested value: -9999.)
    !> @param[out] rlat (npts) earth latitudes in degrees n if iopt<0
    !> @param[out] rlon (npts) earth longitudes in degrees e if iopt<0
    !> @param[out] xpts (npts) grid x point coordinates if iopt>0
    !> @param[out] ypts (npts) grid y point coordinates if iopt>0
    !> @param[in] npts maximum number of coordinates
    !>
    !> @author GAYNO @date 2015-07-13
    subroutine rot_equid_cylind_egrid_error(iopt,fill,rlat,rlon,xpts,ypts,npts)
        implicit none
        !
        integer,intent(in) :: iopt,npts
        !
        real,intent(in) :: fill
        real,intent(out) :: rlat(npts),rlon(npts)
        real,intent(out) :: xpts(npts),ypts(npts)

        if(iopt.ge.0) then
            rlon=fill
            rlat=fill
        endif
        if(iopt.le.0) then
            xpts=fill
            ypts=fill
        endif
    endsubroutine rot_equid_cylind_egrid_error

    !> Computes the vector rotation sines and
    !> cosines for a rotated equidistant cylindrical grid.
    !>
    !> @param[in] rlon Longitude in degrees.
    !> @param[out] crot Clockwise vector rotation cosines.
    !> @param[out] srot Clockwise vector rotation sines.
    !>
    !> @note
    !> ugrid=crot*uearth-srot*vearth;
    !> vgrid=srot*uearth+crot*vearth)
    !>
    !> @author George Gayno
    !> @date Jan 2015
    subroutine rot_equid_cylind_egrid_vect_rot(rlon,crot,srot)
        implicit none

        real,intent(in) :: rlon
        real,intent(out) :: crot,srot

        real(KIND=kd)                   :: slon

        if(irot.eq.1) then
            if(clatr.le.0) then
                crot=real(-sign(1._kd,slatr*slat0))
                srot=0.
            else
                slon=sin((rlon-rlon0)/dpr)
                crot=real((clat0*clat+slat0*slat*clon)/clatr)
                srot=real(slat0*slon/clatr)
            endif
        else
            crot=1.
            srot=0.
        endif

    endsubroutine rot_equid_cylind_egrid_vect_rot

    !> Computes the map jacobians for a rotated equidistant cylindrical grid.
    !>
    !> @param[in] fill Fill value for undefined points.
    !> @param[in] rlon Longitude in degrees.
    !> @param[out] xlon dx/dlon in 1/degrees.
    !> @param[out] xlat dx/dlat in 1/degrees.
    !> @param[out] ylon dy/dlon in 1/degrees.
    !> @param[out] ylat dy/dlat in 1/degrees.
    !>
    !> @author George Gayno
    !> @date Jan 2015
    subroutine rot_equid_cylind_egrid_map_jacob(fill,rlon, &
                                                xlon,xlat,ylon,ylat)
        implicit none

        real,intent(in) :: fill,rlon
        real,intent(out) :: xlon,xlat,ylon,ylat

        real(KIND=kd)                   :: slon,term1,term2
        real(KIND=kd)                   :: xlatf,xlonf,ylatf,ylonf

        if(clatr.le.0._kd) then
            xlon=fill
            xlat=fill
            ylon=fill
            ylat=fill
        else
            slon=sin((rlon-rlon0)/dpr)
            term1=(clat0*clat+slat0*slat*clon)/clatr
            term2=slat0*slon/clatr
            xlonf=term1*clat/(dlons*clatr)
            xlatf=-term2/(dlons*clatr)
            ylonf=term2*clat/dlats
            ylatf=term1/dlats
            xlon=real(xlonf-ylonf)
            xlat=real(xlatf-ylatf)
            ylon=real(xlonf+ylonf)
            ylat=real(xlatf+ylatf)
        endif

    endsubroutine rot_equid_cylind_egrid_map_jacob

    !> Computes the grid box area for a rotated equidistant cylindrical grid.
    !>
    !> @param[in] fill Fill value for undefined points.
    !> @param[out] area Area weights in m^2.
    !>
    !> @author George Gayno
    !> @date Jan 2015
    subroutine rot_equid_cylind_egrid_grid_area(fill,area)
        implicit none

        real,intent(in) :: fill
        real,intent(out) :: area

        if(clatr.le.0._kd) then
            area=fill
        else
            area=real(rerth**2*clatr*dlats*dlons)*2/dpr**2
        endif

    endsubroutine rot_equid_cylind_egrid_grid_area

endmodule ip_rot_equid_cylind_egrid_mod

