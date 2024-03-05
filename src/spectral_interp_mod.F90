!> @file
!! @brief Interpolate spectral.
!! @author Mark Iredell @date 96-04-10

!> @brief Interpolate spectral.
!!
!! @author Mark Iredell @date 96-04-10
module spectral_interp_mod
    use gdswzd_mod
    use ip_grid_mod
    use ip_grid_descriptor_mod
    use ip_grid_factory_mod
    use earth_radius_mod
    use sp_mod
    implicit none

    private
    public :: interpolate_spectral

    interface interpolate_spectral
        module procedure interpolate_spectral_scalar
        module procedure interpolate_spectral_vector
    endinterface interpolate_spectral

    interface polates4
        module procedure polates4_grib1
        module procedure polates4_grib2
    endinterface polates4

    interface polatev4
        module procedure polatev4_grib1
        module procedure polatev4_grib2
    endinterface polatev4

contains

    !> Interpolate spectral scalar.
    !>
    !> @param[in] ipopt interpolation options; ipopt(1)=0 for triangular;
    !> ipopt(1)=1 for rhomboidal; ipopt(2) is truncation number
    !> (defaults to a sensible truncation if ipopt(2)=-1).
    !> @param[in] grid_in input grid descriptor.
    !> @param[in] grid_out output grid descriptor.
    !> @param[in] MI skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1.
    !> @param[in] MO skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1.
    !> @param[in] KM number of fields to interpolate.
    !> @param[in] IBI input bitmap flags (Must be all 0. Routine does
    !> not do bitmapped interpolation.)
    !> @param[in] GI input fields to interpolate.
    !> @param[out] NO number of output points.
    !> @param[inout] RLAT output latitudes in degrees.
    !> @param[inout] RLON output longitudes in degrees.
    !> @param[out] IBO output bitmap flags.
    !> @param[out] LO output bitmaps.
    !> @param[out] GO output fields interpolated.
    !> @param[out] IRET return code. 0/non-0 - successful/not successful.
    !>
  !! @author Mark Iredell @date 96-04-10
    subroutine interpolate_spectral_scalar(ipopt,grid_in,grid_out, &
                                           mi,mo,km,ibi,gi, &
                                           no,rlat,rlon,ibo,lo,go,iret)
        integer,intent(in) :: ipopt(20)
        class(ip_grid),intent(in) :: grid_in,grid_out
        integer,intent(in) :: mi,mo
        integer,intent(in) :: ibi(km),km
        integer,intent(out) :: ibo(km),iret,no
        !
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: gi(mi,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: go(mo,km)

        select type(desc_in=>grid_in%descriptor)
        type is(grib1_descriptor)
            select type(desc_out=>grid_out%descriptor)
            type is(grib1_descriptor)
                call polates4(ipopt,desc_in%gds,desc_out%gds,mi,mo,km,ibi,gi,no,rlat,rlon,ibo,lo,go,iret)
            endselect

        type is(grib2_descriptor)
            select type(desc_out=>grid_out%descriptor)
            type is(grib2_descriptor)
                call polates4(ipopt,desc_in%gdt_num,desc_in%gdt_tmpl,desc_in%gdt_len, &
                              desc_out%gdt_num,desc_out%gdt_tmpl,desc_out%gdt_len, &
                              mi,mo,km,ibi,gi,no,rlat,rlon,ibo,lo,go,iret)
            endselect
        endselect
    endsubroutine interpolate_spectral_scalar

    !> Interpolate spectral vector.
    !>
    !> @param ipopt interpolation options; ipopt(1)=0 for triangular;
    !> ipopt(1)=1 for rhomboidal; ipopt(2) is truncation number
    !> (defaults to a sensible truncation if ipopt(2)=-1).
    !> @param grid_in input grid descriptor.
    !> @param grid_out output grid descriptor.
    !> @param MI skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1.
    !> @param MO skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1.
    !> @param KM number of fields to interpolate.
    !> @param IBI input bitmap flags (Must be all 0. Routine does not do
    !> bitmapped interpolation.)
    !> @param UI input u-component fields to interpolate.
    !> @param VI input v-component fields to interpolate.
    !> @param NO number of output points.
    !> @param RLAT output latitudes in degrees.
    !> @param RLON output longitudes in degrees.
    !> @param CROT vector rotation cosines.
    !> @param SROT vector rotation sines.
    !> @param IBO output bitmap flags.
    !> @param LO output bitmaps.
    !> @param UO output u-component fields interpolated.
    !> @param VO output v-component fields interpolated.
    !> @param IRET return code. 0/non-0 - successful/not successful.
    !>
  !! @author Mark Iredell @date 96-04-10
    subroutine interpolate_spectral_vector(ipopt,grid_in,grid_out, &
                                           mi,mo,km,ibi,ui,vi, &
                                           no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
        class(ip_grid),intent(in) :: grid_in,grid_out
        integer,intent(in) :: ipopt(20),ibi(km)
        integer,intent(in) :: km,mi,mo
        integer,intent(out) :: iret,ibo(km),no
        !
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: ui(mi,km),vi(mi,km)
        real,intent(out) :: uo(mo,km),vo(mo,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: crot(mo),srot(mo)

        select type(desc_in=>grid_in%descriptor)
        type is(grib1_descriptor)
            select type(desc_out=>grid_out%descriptor)
            type is(grib1_descriptor)
                call polatev4_grib1(ipopt,desc_in%gds,desc_out%gds,mi,mo,km,ibi,ui,vi, &
                                    no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
            endselect

        type is(grib2_descriptor)
            select type(desc_out=>grid_out%descriptor)
            type is(grib2_descriptor)
                call polatev4(ipopt,desc_in%gdt_num,desc_in%gdt_tmpl,desc_in%gdt_len, &
                              desc_out%gdt_num,desc_out%gdt_tmpl,desc_out%gdt_len, &
                              mi,mo,km,ibi,ui,vi, &
                              no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
            endselect
        endselect

    endsubroutine interpolate_spectral_vector

    !> Interpolate scalar fields (spectral).
    !>
    !> This subprogram performs spectral interpolation from any grid to
    !> any grid for scalar fields. It requires that the input fields be
    !> uniformly global. Options allow choices between triangular shape
    !> (ipopt(1)=0) and rhomboidal shape (ipopt(1)=1) which has no
    !> default; a second option is the truncation (ipopt(2)) which
    !> defaults to a sensible truncation for the input grid (if
    !> opt(2)=-1).
    !>
    !> @note If the output grid is not found in a special list, then the
    !> transform back to grid is not very fast. This special list
    !> contains global cylindrical grids, polar stereographic grids
    !> centered at the pole and mercator grids.
    !>
    !> Only horizontal interpolation is performed.
    !>
    !> The code recognizes the following projections, where "igdtnumi/o"
    !> is the GRIB2 grid defintion template number for the input and
    !> onutput grids, respectively:
    !> - igdtnumi/o = 00 equidistant cylindrical
    !> - igdtnumo = 01 rotated equidistant cylindrical. "e" and non-"e" staggered
    !> - igdtnumo = 10 mercator cylindrical
    !> - igdtnumo = 20 polar stereographic azimuthal
    !> - igdtnumo = 30 lambert conformal conical
    !> - igdtnumi/o = 40 gaussian cylindrical
    !>
    !> As an added bonus the number of output grid points and their
    !> latitudes and longitudes are also returned. On the other hand,
    !> the output can be a set of station points if igdtnumo < 0, in which
    !> case the number of points and their latitudes and longitudes must
    !> be input. Output bitmaps will not be created.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !>   96-04-10 | Iredell | initial
    !> 2001-06-18 | Iredell | improve detection of special fast transform
    !> 2015-01-27 | Gayno | replace calls to gdswiz with new merged version of gdswzd.
    !> 2015-07-13 | Gayno | convert to grib 2. replace grib 1 kgds arrays with grib 2 grid definition template arrays.
    !>
    !> @param[in] ipopt (20) interpolation options; ipopt(1)=0 for
    !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
    !> number (defaults to sensible if ipopt(2)=-1).
    !> @param[in] igdtnumi grid definition template number - input
    !> grid. Corresponds to the gfld%igdtnum component of the
    !> [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) library
    !> gridmod data structure.
    !> - 00 - equidistant cylindrical
    !> - 01 - rotated equidistant cylindrical.  "e" and non-"e" staggered
    !> - 10 - mercator cyclindrical
    !> - 20 - polar stereographic azimuthal
    !> - 30 - lambert conformal conical
    !> - 40 - gaussian equidistant cyclindrical
    !> @param[in] igdtmpli (igdtleni) grid definition template array -
    !> input grid. corresponds to the gfld%igdtmpl component of the ncep
    !> g2 library gridmod data structure: (section 3 info). See comments
    !> in routine ipolates() for complete definition.
    !> @param[in] igdtleni number of elements of the grid definition
    !> template array - input grid.  corresponds to the gfld%igdtlen
    !> component of the ncep g2 library gridmod data structure.
    !> @param[in] igdtnumo grid definition template number - output
    !> grid. Corresponds to the gfld%igdtnum component of the
    !> [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) library
    !> gridmod data structure. igdtnumo<0 means interpolate to random
    !> station points. Otherwise, same definition as igdtnumi.
    !> @param[in] igdtmplo (igdtleno) grid definition template array -
    !> output grid. Corresponds to the gfld%igdtmpl component of the
    !> ncep g2 library gridmod data structure (section 3 info). See
    !> comments in routine ipolates() for complete definition.
    !> @param[in] igdtleno number of elements of the grid definition
    !> template array - output grid. Corresponds to the gfld%igdtlen
    !> component of the
    !> [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2) library
    !> gridmod data structure.
    !> @param[in] mi skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1
    !> @param[in] mo skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1
    !> @param[out] km number of fields to interpolate
    !> @param[out] ibi (km) input bitmap flags (must be all 0)
    !> @param[out] gi (mi,km) input fields to interpolate
    !> @param[out] no number of output points (only if igdtnumo>=0)
    !> @param[out] rlat (mo) output latitudes in degrees (if igdtnumo<0)
    !> @param[out] rlon (mo) output longitudes in degrees (if igdtnumo<0)
    !> @param[out] ibo (km) output bitmap flags
    !> @param[out] lo (mo,km) output bitmaps (always output)
    !> @param[out] go (mo,km) output fields interpolated
    !> @param[out] iret return code
    !> - 0 successful interpolation
    !> - 2 unrecognized input grid or no grid overlap
    !> - 3 unrecognized output grid
    !> - 41 invalid nonglobal input grid
    !> - 42 invalid spectral method parameters
    !>
    !>! @author Mark Iredell @date 96-04-10
    subroutine polates4_grib2(ipopt,igdtnumi,igdtmpli,igdtleni, &
                              igdtnumo,igdtmplo,igdtleno, &
                              mi,mo,km,ibi,gi, &
                              no,rlat,rlon,ibo,lo,go,iret)
        integer,intent(in) :: igdtnumi,igdtleni
        integer,intent(in) :: igdtmpli(igdtleni)
        integer,intent(in) :: igdtnumo,igdtleno
        integer,intent(in) :: igdtmplo(igdtleno)
        integer,intent(in) :: ipopt(20)
        integer,intent(in) :: mi,mo
        integer,intent(in) :: ibi(km),km
        integer,intent(out) :: ibo(km),iret,no
        !
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: gi(mi,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: go(mo,km)
        !
        real,parameter     :: fill=-9999.
        real,parameter     :: pi=3.14159265358979
        real,parameter     :: dpr=180./pi
        !
        integer                         :: idrti,idrto,ig,jg,im,jm
        integer                         :: igo,jgo,imo,jmo
        integer                         :: iscan,jscan,nscan
        integer                         :: iscano,jscano,nscano
        integer                         :: iskipi,jskipi,iscale
        integer                         :: imaxi,jmaxi,ispec
        integer                         :: ip,iprime,iproj,iromb,k
        integer                         :: maxwv,n,ni,nj,nps
        !
        real                            :: de,dr,dy
        real                            :: dlat,dlon,dlato,dlono
        real                            :: go2(mo,km),h,hi,hj
        real                            :: orient,slat,rerth,e2
        real                            :: rlat1,rlon1,rlat2,rlon2,rlati
        real                            :: xmesh,xp,yp
        real                            :: xpts(mo),ypts(mo)

        type(grib2_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out

        desc_in=init_descriptor(igdtnumi,igdtleni,igdtmpli)
        desc_out=init_descriptor(igdtnumo,igdtleno,igdtmplo)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        iret=0
        if(igdtnumo.ge.0) then
            !CALL GDSWZD(IGDTNUMO,IGDTMPLO,IGDTLENO, 0,MO,FILL,XPTS,YPTS,RLON,RLAT,NO)
            call gdswzd(grid_out,0,mo,fill,xpts,ypts,rlon,rlat,no)
            if(no.eq.0) iret=3
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  AFFIRM APPROPRIATE INPUT GRID
        !    LAT/LON OR GAUSSIAN
        !    NO BITMAPS
        !    FULL ZONAL COVERAGE
        !    FULL MERIDIONAL COVERAGE
        idrti=igdtnumi
        if(idrti.eq.40) idrti=4
        if(idrti.eq.0.or.idrti.eq.4) then
            im=igdtmpli(8)
            jm=igdtmpli(9)
            iscale=igdtmpli(10)*igdtmpli(11)
            if(iscale.eq.0) iscale=10**6
            rlon1=float(igdtmpli(13))/float(iscale)
            rlon2=float(igdtmpli(16))/float(iscale)
            iscan=mod(igdtmpli(19)/128,2)
            jscan=mod(igdtmpli(19)/64,2)
            nscan=mod(igdtmpli(19)/32,2)
        else
            iret=41
        endif
        do k=1,km
            if(ibi(k).ne.0) iret=41
        enddo
        if(iret.eq.0) then
            if(iscan.eq.0) then
                dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
            else
                dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
            endif
            ig=nint(360/abs(dlon))
            iprime=1+mod(-nint(rlon1/dlon)+ig,ig)
            imaxi=ig
            jmaxi=jm
            if(mod(ig,2).ne.0.or.im.lt.ig) iret=41
        endif
        if(iret.eq.0.and.idrti.eq.0) then
            iscale=igdtmpli(10)*igdtmpli(11)
            if(iscale.eq.0) iscale=10**6
            rlat1=float(igdtmpli(12))/float(iscale)
            rlat2=float(igdtmpli(15))/float(iscale)
            dlat=(rlat2-rlat1)/(jm-1)
            jg=nint(180/abs(dlat))
            if(jm.eq.jg) idrti=256
            if(jm.ne.jg.and.jm.ne.jg+1) iret=41
        elseif(iret.eq.0.and.idrti.eq.4) then
            jg=igdtmpli(18)*2
            if(jm.ne.jg) iret=41
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        if(iret.eq.0) then
            iromb=ipopt(1)
            maxwv=ipopt(2)
            if(maxwv.eq.-1) then
                if(iromb.eq.0.and.idrti.eq.4) maxwv=(jmaxi-1)
                if(iromb.eq.1.and.idrti.eq.4) maxwv=(jmaxi-1)/2
                if(iromb.eq.0.and.idrti.eq.0) maxwv=(jmaxi-3)/2
                if(iromb.eq.1.and.idrti.eq.0) maxwv=(jmaxi-3)/4
                if(iromb.eq.0.and.idrti.eq.256) maxwv=(jmaxi-1)/2
                if(iromb.eq.1.and.idrti.eq.256) maxwv=(jmaxi-1)/4
            endif
            if((iromb.ne.0.and.iromb.ne.1).or.maxwv.lt.0) iret=42
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  INTERPOLATE
        if(iret.eq.0) then
            if(nscan.eq.0) then
                iskipi=1
                jskipi=im
            else
                iskipi=jm
                jskipi=1
            endif
            if(iscan.eq.1) iskipi=-iskipi
            if(jscan.eq.0) jskipi=-jskipi
            ispec=0
            !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
            if((igdtnumo.eq.0.or.igdtnumo.eq.40).and. &
               mod(igdtmplo(8),2).eq.0.and.igdtmplo(13).eq.0.and.igdtmplo(19).eq.0) then
                idrto=igdtnumo
                if(idrto.eq.40) idrto=4
                imo=igdtmplo(8)
                jmo=igdtmplo(9)
                iscale=igdtmplo(10)*igdtmplo(11)
                if(iscale.eq.0) iscale=10**6
                rlon2=float(igdtmplo(16))/float(iscale)
                dlono=(mod(rlon2-1+3600,360.)+1)/(imo-1)
                igo=nint(360/abs(dlono))
                if(imo.eq.igo.and.idrto.eq.0) then
                    rlat1=float(igdtmplo(12))/float(iscale)
                    rlat2=float(igdtmplo(15))/float(iscale)
                    dlat=(rlat2-rlat1)/(jmo-1)
                    jgo=nint(180/abs(dlat))
                    if(jmo.eq.jgo) idrto=256
                    if(jmo.eq.jgo.or.jmo.eq.jgo+1) ispec=1
                elseif(imo.eq.igo.and.idrto.eq.4) then
                    jgo=igdtmplo(18)*2
                    if(jmo.eq.jgo) ispec=1
                endif
                if(ispec.eq.1) then
                    call sptrun(iromb,maxwv,idrti,imaxi,jmaxi,idrto,imo,jmo, &
                                km,iprime,iskipi,jskipi,mi,0,0,mo,0,gi,go)
                endif
                !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
            elseif(igdtnumo.eq.20.and. &
                   igdtmplo(8).eq.igdtmplo(9).and.mod(igdtmplo(8),2).eq.1.and. &
                   igdtmplo(15).eq.igdtmplo(16).and.igdtmplo(18).eq.64) then
                nps=igdtmplo(8)
                rlat1=float(igdtmplo(10))*1.e-6
                rlon1=float(igdtmplo(11))*1.e-6
                orient=float(igdtmplo(14))*1.e-6
                xmesh=float(igdtmplo(15))*1.e-3
                iproj=mod(igdtmplo(17)/128,2)
                ip=(nps+1)/2
                h=(-1.)**iproj
                slat=float(abs(igdtmplo(13)))*1.e-6
                call earth_radius(igdtmplo,igdtleno,rerth,e2)
                de=(1.+sin(slat/dpr))*rerth
                dr=de*cos(rlat1/dpr)/(1+h*sin(rlat1/dpr))
                xp=1-h*sin((rlon1-orient)/dpr)*dr/xmesh
                yp=1+cos((rlon1-orient)/dpr)*dr/xmesh
                if(nint(xp).eq.ip.and.nint(yp).eq.ip) then
                    if(iproj.eq.0) then
                        call sptruns(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                     iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                     slat,xmesh,orient,gi,go,go2)
                    else
                        call sptruns(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                     iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                     slat,xmesh,orient,gi,go2,go)
                    endif
                    ispec=1
                endif
                !  SPECIAL CASE OF MERCATOR GRID
            elseif(igdtnumo.eq.10) then
                ni=igdtmplo(8)
                nj=igdtmplo(9)
                rlat1=float(igdtmplo(10))*1.0e-6
                rlon1=float(igdtmplo(11))*1.0e-6
                rlon2=float(igdtmplo(15))*1.0e-6
                rlati=float(igdtmplo(13))*1.0e-6
                iscano=mod(igdtmplo(16)/128,2)
                jscano=mod(igdtmplo(16)/64,2)
                nscano=mod(igdtmplo(16)/32,2)
                dy=float(igdtmplo(19))*1.0e-3
                hi=(-1.)**iscano
                hj=(-1.)**(1-jscano)
                call earth_radius(igdtmplo,igdtleno,rerth,e2)
                dlono=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(ni-1)
                dlato=hj*dy/(rerth*cos(rlati/dpr))*dpr
                if(nscano.eq.0) then
                    call sptrunm(iromb,maxwv,idrti,imaxi,jmaxi,km,ni,nj, &
                                 iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                 rlat1,rlon1,dlato,dlono,gi,go)
                    ispec=1
                endif
            endif
            !  GENERAL SLOW CASE
            if(ispec.eq.0) then
                call sptrung(iromb,maxwv,idrti,imaxi,jmaxi,km,no, &
                             iprime,iskipi,jskipi,mi,mo,0,0,0,rlat,rlon,gi,go)
            endif
            do k=1,km
                ibo(k)=0
                do n=1,no
                    lo(n,k)=.true.
                enddo
            enddo
        else
            do k=1,km
                ibo(k)=1
                do n=1,no
                    lo(n,k)=.false.
                    go(n,k)=0.
                enddo
            enddo
        endif
    endsubroutine polates4_grib2

    !> Interpolate scalar fields (spectral).
    !>
    !> This subprogram performs spectral interpolation from any grid to
    !> any grid for scalar fields.  It requires that the input fields be
    !> uniformly global.
    !>
    !> Options allow choices between triangular shape (ipopt(1)=0) and
    !> rhomboidal shape (ipopt(1)=1) which has no default; a second
    !> option is the truncation (ipopt(2)) which defaults to a sensible
    !> truncation for the input grid (if opt(2)=-1).
    !>
    !> @note If the output grid is not found in a special list, then the
    !> transform back to grid is not very fast. This special list
    !> contains global cylindrical grids, polar stereographic grids
    !> centered at the pole and mercator grids.
    !>
    !> Only horizontal interpolation is performed. The grids are defined
    !> by their grid description sections (passed in integer form as
    !> decoded by subprogram w3fi63()).
    !>
    !> The current code recognizes the following projections:
    !> - kgds(1) = 000 equidistant cylindrical
    !> - kgds(1) = 001 mercator cylindrical
    !> - kgds(1) = 003 lambert conformal conical
    !> - kgds(1) = 004 gaussian cylindrical (spectral native)
    !> - kgds(1) = 005 polar stereographic azimuthal
    !> - kgds(1) = 203 rotated equidistant cylindrical (e-stagger)
    !> - kgds(1) = 205 rotated equidistant cylindrical (b-stagger)
    !>
    !> Where kgds could be either input kgdsi or output kgdso. As an
    !> added bonus the number of output grid points and their latitudes
    !> and longitudes are also returned. On the other hand, the output
    !> can be a set of station points if kgdso(1)<0, in which case the
    !> number of points and their latitudes and longitudes must be
    !> input. Output bitmaps will not be created.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | Iredell | Initial
    !> 2001-06-18 | Iredell | improve detection of special fast transform
    !> 2015-01-27 | Gayno | replace calls to gdswiz() with new merged version of gdswzd().
    !>
    !> @param[in] ipopt (20) interpolation options ipopt(1)=0 for
    !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
    !> number (defaults to sensible if ipopt(2)=-1).
    !> @param[in] kgdsi (200) input gds parameters as decoded by w3fi63
    !> @param[in] kgdso (200) output gds parameters (kgdso(1)<0 implies random station points)
    !> @param[in] mi skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1
    !> @param[in] mo skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1
    !> @param[in] km number of fields to interpolate
    !> @param[in] ibi (km) input bitmap flags (must be all 0)
    !> @param[in] gi (mi,km) input fields to interpolate
    !> @param[out] no number of output points (only if kgdso(1)<0)
    !> @param[out] rlat (no) output latitudes in degrees (if kgdso(1)<0)
    !> @param[out] rlon (no) output longitudes in degrees (if kgdso(1)<0)
    !> @param[out] ibo (km) output bitmap flags
    !> @param[out] lo (mo,km) output bitmaps (always output)
    !> @param[out] go (mo,km) output fields interpolated
    !> @param[out] iret return code
    !> - 0 successful interpolation
    !> - 2 unrecognized input grid or no grid overlap
    !> - 3 unrecognized output grid
    !> - 41 invalid nonglobal input grid
    !> - 42 invalid spectral method parameters
    !>
    !> @author Iredell @date 96-04-10
    subroutine polates4_grib1(ipopt,kgdsi,kgdso,mi,mo,km,ibi,gi, &
                              no,rlat,rlon,ibo,lo,go,iret)
        integer,intent(in) :: ipopt(20),kgdsi(200)
        integer,intent(in) :: kgdso(200),mi,mo
        integer,intent(in) :: ibi(km),km
        integer,intent(out) :: ibo(km),iret
        !
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: gi(mi,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: go(mo,km)
        !
        real,parameter     :: fill=-9999.
        real,parameter     :: rerth=6.3712e6
        real,parameter     :: pi=3.14159265358979
        real,parameter     :: dpr=180./pi
        !
        integer                         :: idrti,idrto,ig,jg,im,jm
        integer                         :: igo,jgo,imo,jmo
        integer                         :: iscan,jscan,nscan
        integer                         :: iscano,jscano,nscano
        integer                         :: iskipi,jskipi
        integer                         :: imaxi,jmaxi,ispec
        integer                         :: ip,iprime,iproj,iromb,k
        integer                         :: maxwv,n,ni,nj,nps,no
        !
        real                            :: de,dr,dy
        real                            :: dlat,dlon,dlato,dlono
        real                            :: go2(mo,km),h,hi,hj
        real                            :: orient
        real                            :: rlat1,rlon1,rlat2,rlon2,rlati
        real                            :: xmesh,xp,yp
        real                            :: xpts(mo),ypts(mo)

        type(grib1_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out

        desc_in=init_descriptor(kgdsi)
        desc_out=init_descriptor(kgdso)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        iret=0
        if(kgdso(1).ge.0) then
            call gdswzd(grid_out,0,mo,fill,xpts,ypts,rlon,rlat,no)
            if(no.eq.0) iret=3
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  AFFIRM APPROPRIATE INPUT GRID
        !    LAT/LON OR GAUSSIAN
        !    NO BITMAPS
        !    FULL ZONAL COVERAGE
        !    FULL MERIDIONAL COVERAGE
        idrti=kgdsi(1)
        im=kgdsi(2)
        jm=kgdsi(3)
        rlon1=kgdsi(5)*1.e-3
        rlon2=kgdsi(8)*1.e-3
        iscan=mod(kgdsi(11)/128,2)
        jscan=mod(kgdsi(11)/64,2)
        nscan=mod(kgdsi(11)/32,2)
        if(idrti.ne.0.and.idrti.ne.4) iret=41
        do k=1,km
            if(ibi(k).ne.0) iret=41
        enddo
        if(iret.eq.0) then
            if(iscan.eq.0) then
                dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
            else
                dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
            endif
            ig=nint(360/abs(dlon))
            iprime=1+mod(-nint(rlon1/dlon)+ig,ig)
            imaxi=ig
            jmaxi=jm
            if(mod(ig,2).ne.0.or.im.lt.ig) iret=41
        endif
        if(iret.eq.0.and.idrti.eq.0) then
            rlat1=kgdsi(4)*1.e-3
            rlat2=kgdsi(7)*1.e-3
            dlat=(rlat2-rlat1)/(jm-1)
            jg=nint(180/abs(dlat))
            if(jm.eq.jg) idrti=256
            if(jm.ne.jg.and.jm.ne.jg+1) iret=41
        elseif(iret.eq.0.and.idrti.eq.4) then
            jg=kgdsi(10)*2
            if(jm.ne.jg) iret=41
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        if(iret.eq.0) then
            iromb=ipopt(1)
            maxwv=ipopt(2)
            if(maxwv.eq.-1) then
                if(iromb.eq.0.and.idrti.eq.4) maxwv=(jmaxi-1)
                if(iromb.eq.1.and.idrti.eq.4) maxwv=(jmaxi-1)/2
                if(iromb.eq.0.and.idrti.eq.0) maxwv=(jmaxi-3)/2
                if(iromb.eq.1.and.idrti.eq.0) maxwv=(jmaxi-3)/4
                if(iromb.eq.0.and.idrti.eq.256) maxwv=(jmaxi-1)/2
                if(iromb.eq.1.and.idrti.eq.256) maxwv=(jmaxi-1)/4
            endif
            if((iromb.ne.0.and.iromb.ne.1).or.maxwv.lt.0) iret=42
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  INTERPOLATE
        if(iret.eq.0) then
            if(nscan.eq.0) then
                iskipi=1
                jskipi=im
            else
                iskipi=jm
                jskipi=1
            endif
            if(iscan.eq.1) iskipi=-iskipi
            if(jscan.eq.0) jskipi=-jskipi
            ispec=0
            !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
            if((kgdso(1).eq.0.or.kgdso(1).eq.4).and. &
               mod(kgdso(2),2).eq.0.and.kgdso(5).eq.0.and.kgdso(11).eq.0) then
                idrto=kgdso(1)
                imo=kgdso(2)
                jmo=kgdso(3)
                rlon2=kgdso(8)*1.e-3
                dlono=(mod(rlon2-1+3600,360.)+1)/(imo-1)
                igo=nint(360/abs(dlono))
                if(imo.eq.igo.and.idrto.eq.0) then
                    rlat1=kgdso(4)*1.e-3
                    rlat2=kgdso(7)*1.e-3
                    dlat=(rlat2-rlat1)/(jmo-1)
                    jgo=nint(180/abs(dlat))
                    if(jmo.eq.jgo) idrto=256
                    if(jmo.eq.jgo.or.jmo.eq.jgo+1) ispec=1
                elseif(imo.eq.igo.and.idrto.eq.4) then
                    jgo=kgdso(10)*2
                    if(jmo.eq.jgo) ispec=1
                endif
                if(ispec.eq.1) then
                    call sptrun(iromb,maxwv,idrti,imaxi,jmaxi,idrto,imo,jmo, &
                                km,iprime,iskipi,jskipi,mi,0,0,mo,0,gi,go)
                endif
                !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
            elseif(kgdso(1).eq.5.and. &
                   kgdso(2).eq.kgdso(3).and.mod(kgdso(2),2).eq.1.and. &
                   kgdso(8).eq.kgdso(9).and.kgdso(11).eq.64) then
                nps=kgdso(2)
                rlat1=kgdso(4)*1.e-3
                rlon1=kgdso(5)*1.e-3
                orient=kgdso(7)*1.e-3
                xmesh=kgdso(8)
                iproj=mod(kgdso(10)/128,2)
                ip=(nps+1)/2
                h=(-1.)**iproj
                de=(1.+sin(60./dpr))*rerth
                dr=de*cos(rlat1/dpr)/(1+h*sin(rlat1/dpr))
                xp=1-h*sin((rlon1-orient)/dpr)*dr/xmesh
                yp=1+cos((rlon1-orient)/dpr)*dr/xmesh
                if(nint(xp).eq.ip.and.nint(yp).eq.ip) then
                    if(iproj.eq.0) then
                        call sptruns(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                     iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                     60.,xmesh,orient,gi,go,go2)
                    else
                        call sptruns(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                     iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                     60.,xmesh,orient,gi,go2,go)
                    endif
                    ispec=1
                endif
                !  SPECIAL CASE OF MERCATOR GRID
            elseif(kgdso(1).eq.1) then
                ni=kgdso(2)
                nj=kgdso(3)
                rlat1=kgdso(4)*1.e-3
                rlon1=kgdso(5)*1.e-3
                rlon2=kgdso(8)*1.e-3
                rlati=kgdso(9)*1.e-3
                iscano=mod(kgdso(11)/128,2)
                jscano=mod(kgdso(11)/64,2)
                nscano=mod(kgdso(11)/32,2)
                dy=kgdso(13)
                hi=(-1.)**iscano
                hj=(-1.)**(1-jscano)
                dlono=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(ni-1)
                dlato=hj*dy/(rerth*cos(rlati/dpr))*dpr
                if(nscano.eq.0) then
                    call sptrunm(iromb,maxwv,idrti,imaxi,jmaxi,km,ni,nj, &
                                 iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                 rlat1,rlon1,dlato,dlono,gi,go)
                    ispec=1
                endif
            endif
            !  GENERAL SLOW CASE
            if(ispec.eq.0) then
                call sptrung(iromb,maxwv,idrti,imaxi,jmaxi,km,no, &
                             iprime,iskipi,jskipi,mi,mo,0,0,0,rlat,rlon,gi,go)
            endif
            do k=1,km
                ibo(k)=0
                do n=1,no
                    lo(n,k)=.true.
                enddo
            enddo
        else
            do k=1,km
                ibo(k)=1
                do n=1,no
                    lo(n,k)=.false.
                    go(n,k)=0.
                enddo
            enddo
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    endsubroutine polates4_grib1

    !> Interpolate vector fields (spectral).
    !>
    !> This subprogram performs spectral interpolation from any grid to
    !> any grid for vector fields. It requires that the input fields be
    !> uniformly global. Options allow choices between triangular shape
    !> (ipopt(1)=0) and rhomboidal shape (ipopt(1)=1) which has no
    !> default; a second option is the truncation (ipopt(2)) which
    !> defaults to a sensible truncation for the input grid (if
    !> opt(2)=-1).
    !>
    !> @note If the output grid is not found in a special list, then the
    !> transform back to grid is not very fast.  This special list
    !> contains global cylindrical grids, polar stereographic grids
    !> centered at the pole and mercator grids. Only horizontal
    !> interpolation is performed.
    !>
    !> The input and output grids are defined by their grib 2 grid
    !> definition template as decoded by the ncep g2 library. The code
    !> recognizes the following projections, where "igdtnumi/o" is the
    !> grib 2 grid defintion template number for the input and output
    !> grids, respectively:
    !> - igdtnumi/o=00 equidistant cylindrical
    !> - igdtnumo  =01 rotated equidistant cylindrical. "e" and non-"e" staggered
    !> - igdtnumo  =10 mercator cylindrical
    !> - igdtnumo  =20 polar stereographic azimuthal
    !> - igdtnumo  =30 lambert conformal conical
    !> - igdtnumi/o=40 gaussian cylindrical
    !>
    !> The input and output vectors are rotated so that they are either
    !> resolved relative to the defined grid in the direction of
    !> increasing x and y coordinates or resolved relative to easterly
    !> and northerly directions, as designated by their respective grid
    !> description sections.
    !>
    !> As an added bonus the number of output grid points and their
    !> latitudes and longitudes are also returned along with their
    !> vector rotation parameters.  On the other hand, the output can be
    !> a set of station points if igdtnumo<0, in which case the number
    !> of points and their latitudes and longitudes must be input along
    !> with their vector rotation parameters.
    !>
    !> Output bitmaps will only be created when the output grid extends
    !> outside of the domain of the input grid.  the output field is set
    !> to 0 where the output bitmap is off.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | iredell | initial
    !> 2001-06-18 | iredell | improve detection of special fast transform
    !> 2015-01-27 | gayno | replace calls to gdswiz() with new merged routine gdswzd().
    !> 2015-07-13 | gayno | convert to grib 2. replace grib 1 kgds arrays with grib 2 grid definition template arrays.
    !>
    !> @param[in] ipopt (20) interpolation options ipopt(1)=0 for
    !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
    !> number (defaults to sensible if ipopt(2)=-1).
    !> @param[in] igdtnumi grid definition template number - input
    !> grid. Corresponds to the gfld%igdtnum component of the ncep g2
    !> library gridmod data structure:
    !> - 00 equidistant cylindrical
    !> - 01 rotated equidistant cylindrical.  "e" and non-"e" staggered
    !> - 10 mercator cyclindrical
    !> - 20 polar stereographic azimuthal
    !> - 30 lambert conformal conical
    !> - 40 gaussian equidistant cyclindrical
    !> @param[in] igdtmpli (igdtleni) grid definition template array - input
    !> grid. corresponds to the gfld%igdtmpl component of the ncep g2
    !> library gridmod data structure (section 3 info).  see comments in
    !> routine ipolatev for complete definition.
    !> @param[in] igdtleni number of elements of the grid definition
    !> template array - input grid.  corresponds to the gfld%igdtlen
    !> component of the ncep g2 library gridmod data structure.
    !> @param[in] igdtnumo grid definition template number - output
    !> grid. Corresponds to the gfld%igdtnum component of the ncep g2
    !> library gridmod data structure. igdtnumo<0 means interpolate to
    !> random station points. Otherwise, same definition as "igdtnumi".
    !> @param[in] igdtmplo (igdtleno) grid definition template array -
    !> output grid. corresponds to the gfld%igdtmpl component of the
    !> ncep g2 library gridmod data structure (section 3 info).  see
    !> comments in routine ipolatev() for complete definition.
    !> @param[in] igdtleno number of elements of the grid definition
    !> template array - output grid. Corresponds to the gfld%igdtlen
    !> component of the ncep g2 library gridmod data structure.
    !> @param[in] mi skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1.
    !> @param[in] mo skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1.
    !> @param[in] km number of fields to interpolate
    !> @param[in] ibi (km) input bitmap flags (must be all 0)
    !> @param[in] ui (mi,km) input u-component fields to interpolate
    !> @param[in] vi (mi,km) input v-component fields to interpolate
    !> @param[out] no number of output points (only if igdtnumo>=0)
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
    !> - 0 successful interpolation
    !> - 2 unrecognized input grid or no grid overlap
    !> - 3 unrecognized output grid
    !> - 41 invalid nonglobal input grid
    !> - 42 invalid spectral method parameters
    !>
    !> @author IREDELL @date 96-04-10
    subroutine polatev4_grib2(ipopt,igdtnumi,igdtmpli,igdtleni, &
                              igdtnumo,igdtmplo,igdtleno, &
                              mi,mo,km,ibi,ui,vi, &
                              no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
        integer,intent(in) :: ipopt(20),ibi(km)
        integer,intent(in) :: km,mi,mo
        integer,intent(out) :: iret,ibo(km),no
        integer,intent(in) :: igdtnumi,igdtleni
        integer,intent(in) :: igdtmpli(igdtleni)
        integer,intent(in) :: igdtnumo,igdtleno
        integer,intent(in) :: igdtmplo(igdtleno)
        !
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: ui(mi,km),vi(mi,km)
        real,intent(out) :: uo(mo,km),vo(mo,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: crot(mo),srot(mo)
        !
        real,parameter :: fill=-9999.
        real,parameter :: pi=3.14159265358979
        real,parameter :: dpr=180./pi
        !
        integer                         :: idrto,iromb,iskipi,ispec
        integer                         :: idrti,imaxi,jmaxi,im,jm
        integer                         :: iprime,ig,imo,jmo,igo,jgo
        integer                         :: iscan,jscan,nscan
        integer                         :: iscano,jscano,nscano
        integer                         :: iscale,ip,iproj,jskipi,jg
        integer                         :: k,maxwv,n,ni,nj,nps
        !
        real                            :: dlat,dlon,dlato,dlono,de,dr,dy
        real                            :: dum,e2,h,hi,hj,dumm(1)
        real                            :: orient,rerth,slat
        real                            :: rlat1,rlon1,rlat2,rlon2,rlati
        real                            :: urot,vrot,uo2(mo,km),vo2(mo,km)
        real                            :: xmesh,xp,yp,xpts(mo),ypts(mo)

        type(grib2_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out

        desc_in=init_descriptor(igdtnumi,igdtleni,igdtmpli)
        desc_out=init_descriptor(igdtnumo,igdtleno,igdtmplo)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        iret=0
        if(igdtnumo.ge.0) then
            call gdswzd(grid_out,0,mo,fill,xpts,ypts, &
                        rlon,rlat,no,crot,srot)
            if(no.eq.0) iret=3
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  AFFIRM APPROPRIATE INPUT GRID
        !    LAT/LON OR GAUSSIAN
        !    NO BITMAPS
        !    FULL ZONAL COVERAGE
        !    FULL MERIDIONAL COVERAGE
        idrti=igdtnumi
        if(idrti.eq.40) idrti=4
        if(idrti.eq.0.or.idrti.eq.4) then
            im=igdtmpli(8)
            jm=igdtmpli(9)
            iscale=igdtmpli(10)*igdtmpli(11)
            if(iscale.eq.0) iscale=10**6
            rlon1=float(igdtmpli(13))/float(iscale)
            rlon2=float(igdtmpli(16))/float(iscale)
            iscan=mod(igdtmpli(19)/128,2)
            jscan=mod(igdtmpli(19)/64,2)
            nscan=mod(igdtmpli(19)/32,2)
        else
            iret=41
        endif
        do k=1,km
            if(ibi(k).ne.0) iret=41
        enddo
        if(iret.eq.0) then
            if(iscan.eq.0) then
                dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
            else
                dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
            endif
            ig=nint(360/abs(dlon))
            iprime=1+mod(-nint(rlon1/dlon)+ig,ig)
            imaxi=ig
            jmaxi=jm
            if(mod(ig,2).ne.0.or.im.lt.ig) iret=41
        endif
        if(iret.eq.0.and.idrti.eq.0) then
            iscale=igdtmpli(10)*igdtmpli(11)
            if(iscale.eq.0) iscale=10**6
            rlat1=float(igdtmpli(12))/float(iscale)
            rlat2=float(igdtmpli(15))/float(iscale)
            dlat=(rlat2-rlat1)/(jm-1)
            jg=nint(180/abs(dlat))
            if(jm.eq.jg) idrti=256
            if(jm.ne.jg.and.jm.ne.jg+1) iret=41
        elseif(iret.eq.0.and.idrti.eq.4) then
            jg=igdtmpli(18)*2
            if(jm.ne.jg) iret=41
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        if(iret.eq.0) then
            iromb=ipopt(1)
            maxwv=ipopt(2)
            if(maxwv.eq.-1) then
                if(iromb.eq.0.and.idrti.eq.4) maxwv=(jmaxi-1)
                if(iromb.eq.1.and.idrti.eq.4) maxwv=(jmaxi-1)/2
                if(iromb.eq.0.and.idrti.eq.0) maxwv=(jmaxi-3)/2
                if(iromb.eq.1.and.idrti.eq.0) maxwv=(jmaxi-3)/4
                if(iromb.eq.0.and.idrti.eq.256) maxwv=(jmaxi-1)/2
                if(iromb.eq.1.and.idrti.eq.256) maxwv=(jmaxi-1)/4
            endif
            if((iromb.ne.0.and.iromb.ne.1).or.maxwv.lt.0) iret=42
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  INTERPOLATE
        if(iret.eq.0) then
            if(nscan.eq.0) then
                iskipi=1
                jskipi=im
            else
                iskipi=jm
                jskipi=1
            endif
            if(iscan.eq.1) iskipi=-iskipi
            if(jscan.eq.0) jskipi=-jskipi
            ispec=0
            !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
            if((igdtnumo.eq.0.or.igdtnumo.eq.40).and. &
               mod(igdtmplo(8),2).eq.0.and.igdtmplo(13).eq.0.and. &
               igdtmplo(19).eq.0) then
                idrto=igdtnumo
                if(idrto.eq.40) idrto=4
                imo=igdtmplo(8)
                jmo=igdtmplo(9)
                iscale=igdtmplo(10)*igdtmplo(11)
                if(iscale.eq.0) iscale=10**6
                rlon2=float(igdtmplo(16))/float(iscale)
                dlono=(mod(rlon2-1+3600,360.)+1)/(imo-1)
                igo=nint(360/abs(dlono))
                if(imo.eq.igo.and.idrto.eq.0) then
                    rlat1=float(igdtmplo(12))/float(iscale)
                    rlat2=float(igdtmplo(15))/float(iscale)
                    dlat=(rlat2-rlat1)/(jmo-1)
                    jgo=nint(180/abs(dlat))
                    if(jmo.eq.jgo) idrto=256
                    if(jmo.eq.jgo.or.jmo.eq.jgo+1) ispec=1
                elseif(imo.eq.igo.and.idrto.eq.4) then
                    jgo=igdtmplo(18)*2
                    if(jmo.eq.jgo) ispec=1
                endif
                if(ispec.eq.1) then
                    call sptrunv(iromb,maxwv,idrti,imaxi,jmaxi,idrto,imo,jmo, &
                                 km,iprime,iskipi,jskipi,mi,0,0,mo,0,ui,vi, &
                                 .true.,uo,vo,.false.,dumm,dumm,.false.,dumm,dumm)
                endif
                !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
            elseif(igdtnumo.eq.20.and. &
                   igdtmplo(8).eq.igdtmplo(9).and.mod(igdtmplo(8),2).eq.1.and. &
                   igdtmplo(15).eq.igdtmplo(16).and.igdtmplo(18).eq.64.and. &
                   mod(igdtmplo(12)/8,2).eq.1) then
                nps=igdtmplo(8)
                rlat1=float(igdtmplo(10))*1.e-6
                rlon1=float(igdtmplo(11))*1.e-6
                orient=float(igdtmplo(14))*1.e-6
                xmesh=float(igdtmplo(15))*1.e-3
                iproj=mod(igdtmplo(17)/128,2)
                ip=(nps+1)/2
                h=(-1.)**iproj
                slat=float(abs(igdtmplo(13)))*1.e-6
                call earth_radius(igdtmplo,igdtleno,rerth,e2)
                de=(1.+sin(slat/dpr))*rerth
                dr=de*cos(rlat1/dpr)/(1+h*sin(rlat1/dpr))
                xp=1-h*sin((rlon1-orient)/dpr)*dr/xmesh
                yp=1+cos((rlon1-orient)/dpr)*dr/xmesh
                if(nint(xp).eq.ip.and.nint(yp).eq.ip) then
                    if(iproj.eq.0) then
                        call sptrunsv(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                      iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                      slat,xmesh,orient,ui,vi,.true.,uo,vo,uo2,vo2, &
                                      .false.,dumm,dumm,dumm,dumm, &
                                      .false.,dumm,dumm,dumm,dumm)
                    else
                        call sptrunsv(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                      iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                      slat,xmesh,orient,ui,vi,.true.,uo2,vo2,uo,vo, &
                                      .false.,dumm,dumm,dumm,dumm, &
                                      .false.,dumm,dumm,dumm,dumm)
                    endif
                    ispec=1
                endif
                !  SPECIAL CASE OF MERCATOR GRID
            elseif(igdtnumo.eq.10) then
                ni=igdtmplo(8)
                nj=igdtmplo(9)
                rlat1=float(igdtmplo(10))*1.0e-6
                rlon1=float(igdtmplo(11))*1.0e-6
                rlon2=float(igdtmplo(15))*1.0e-6
                rlati=float(igdtmplo(13))*1.0e-6
                iscano=mod(igdtmplo(16)/128,2)
                jscano=mod(igdtmplo(16)/64,2)
                nscano=mod(igdtmplo(16)/32,2)
                dy=float(igdtmplo(19))*1.0e-3
                hi=(-1.)**iscano
                hj=(-1.)**(1-jscano)
                call earth_radius(igdtmplo,igdtleno,rerth,e2)
                dlono=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(ni-1)
                dlato=hj*dy/(rerth*cos(rlati/dpr))*dpr
                if(nscano.eq.0) then
                    call sptrunmv(iromb,maxwv,idrti,imaxi,jmaxi,km,ni,nj, &
                                  iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                  rlat1,rlon1,dlato,dlono,ui,vi, &
                                  .true.,uo,vo,.false.,dumm,dumm,.false.,dumm,dumm)
                    ispec=1
                endif
            endif
            !  GENERAL SLOW CASE
            if(ispec.eq.0) then
                call sptrungv(iromb,maxwv,idrti,imaxi,jmaxi,km,no, &
                              iprime,iskipi,jskipi,mi,mo,0,0,0,rlat,rlon, &
                              ui,vi,.true.,uo,vo,.false.,dumm,dumm,.false.,dumm,dumm)
                do k=1,km
                    ibo(k)=0
                    do n=1,no
                        lo(n,k)=.true.
                        urot=crot(n)*uo(n,k)-srot(n)*vo(n,k)
                        vrot=srot(n)*uo(n,k)+crot(n)*vo(n,k)
                        uo(n,k)=urot
                        vo(n,k)=vrot
                    enddo
                enddo
            endif
        else
            do k=1,km
                ibo(k)=1
                do n=1,no
                    lo(n,k)=.false.
                    uo(n,k)=0.
                    vo(n,k)=0.
                enddo
            enddo
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    endsubroutine polatev4_grib2

    !> Interpolate vector fields (spectral).
    !>
    !> This subprogram performs spectral interpolation from any grid to
    !> any grid for vector fields. It requires that the input fields be
    !> uniformly global. Options allow choices between triangular shape
    !> (ipopt(1)=0) and rhomboidal shape (ipopt(1)=1) which has no
    !> default; a second option is the truncation (ipopt(2)) which
    !> defaults to a sensible truncation for the input grid (if
    !> opt(2)=-1).
    !>
    !> @note If the output grid is not found in a special list, then the
    !> transform back to grid is not very fast.  This special list
    !> contains global cylindrical grids, polar stereographic grids
    !> centered at the pole and mercator grids.
    !>
    !> Only horizontal interpolation is performed. The grids are defined
    !> by their grid description sections (passed in integer form as
    !> decoded by subprogram w3fi63).
    !>
    !> The current code recognizes the following projections:
    !> - (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
    !> - (KGDS(1)=001) MERCATOR CYLINDRICAL
    !> - (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
    !> - (KGDS(1)=004) GAUSSIAN CYLINDRICAL (SPECTRAL NATIVE)
    !> - (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
    !> - (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL (E-STAGGER)
    !> - (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL (B-STAGGER)
    !>
    !> Where kgds could be either input kgdsi or output kgdso.
    !>
    !> The input and output vectors are rotated so that they are either
    !> resolved relative to the defined grid in the direction of
    !> increasing x and y coordinates or resolved relative to easterly
    !> and northerly directions, as designated by their respective grid
    !> description sections. As an added bonus the number of output grid
    !> points and their latitudes and longitudes are also returned along
    !> with their vector rotation parameters. On the other hand, the
    !> output can be a set of station points if kgdso(1)<0, in which
    !> case the number of points and their latitudes and longitudes must
    !> be input along with their vector rotation parameters.
    !>
    !> Output bitmaps will only be created when the output grid extends
    !> outside of the domain of the input grid. The output field is set
    !> to 0 where the output bitmap is off.
    !>
    !> ### Program History Log
    !> Date | Programmer | Comments
    !> -----|------------|---------
    !> 96-04-10 | iredell | initial.
    !> 2001-06-18 | iredell |  improve detection of special fast transform
    !> 2015-01-27 | gayno | replace calls to gdswiz() with new merged routine gdswzd().
    !>
    !> @param[in] ipopt (20) interpolation options ipopt(1)=0 for
    !> triangular, ipopt(1)=1 for rhomboidal; ipopt(2) is truncation
    !> number (defaults to sensible if ipopt(2)=-1).
    !> @param[in] kgdsi (200) input gds parameters as decoded by w3fi63.
    !> @param[in] kgdso (200) output gds parameters (kgdso(1)<0 implies
    !> random station points).
    !> @param[in] mi skip number between input grid fields if km>1 or
    !> dimension of input grid fields if km=1.
    !> @param[in] mo skip number between output grid fields if km>1 or
    !> dimension of output grid fields if km=1.
    !> @param[in] km number of fields to interpolate
    !> @param[in] ibi (km) input bitmap flags (must be all 0)
    !> @param[in] ui (mi,km) input u-component fields to interpolate
    !> @param[in] vi (mi,km) input v-component fields to interpolate
    !> @param[out] no number of output points (only if kgdso(1)<0)
    !> @param[out] rlat (no) output latitudes in degrees (if kgdso(1)<0)
    !> @param[out] rlon (no) output longitudes in degrees (if kgdso(1)<0)
    !> @param[out] crot (no) vector rotation cosines (if kgdso(1)<0)
    !> @param[out] srot (no) vector rotation sines (if kgdso(1)<0)
    !> (ugrid=crot*uearth-srot*vearth; vgrid=srot*uearth+crot*vearth)
    !> @param[out] ibo (km) output bitmap flags
    !> @param[out] lo (mo,km) output bitmaps (always output)
    !> @param[out] uo (mo,km) output u-component fields interpolated
    !> @param[out] vo (mo,km) output v-component fields interpolated
    !> @param[out] iret return code
    !> - 0 successful interpolation
    !> - 2 unrecognized input grid or no grid overlap
    !> - 3 unrecognized output grid
    !> - 41 invalid nonglobal input grid
    !> - 42 invalid spectral method parameters
    !>
    !> @author IREDELL @date 96-04-10
    subroutine polatev4_grib1(ipopt,kgdsi,kgdso,mi,mo,km,ibi,ui,vi, &
                              no,rlat,rlon,crot,srot,ibo,lo,uo,vo,iret)
        integer,intent(in) :: ipopt(20),ibi(km)
        integer,intent(in) :: km,mi,mo
        integer,intent(out) :: iret,ibo(km)
        integer,intent(in) :: kgdsi(200),kgdso(200)
        !
        logical*1,intent(out) :: lo(mo,km)
        !
        real,intent(in) :: ui(mi,km),vi(mi,km)
        real,intent(out) :: uo(mo,km),vo(mo,km)
        real,intent(inout) :: rlat(mo),rlon(mo)
        real,intent(out) :: crot(mo),srot(mo)
        !
        real,parameter :: fill=-9999.
        real,parameter :: rerth=6.3712e6
        real,parameter :: pi=3.14159265358979
        real,parameter :: dpr=180./pi
        !
        integer                         :: idrto,iromb,iskipi,ispec
        integer                         :: idrti,imaxi,jmaxi,im,jm
        integer                         :: iprime,ig,imo,jmo,igo,jgo
        integer                         :: iscan,jscan,nscan
        integer                         :: iscano,jscano,nscano
        integer                         :: ip,iproj,jskipi,jg
        integer                         :: k,maxwv,n,ni,nj,no,nps
        !
        real                            :: dlat,dlon,dlato,dlono,de,dr,dy
        real                            :: dum,h,hi,hj,dumm(1)
        real                            :: orient
        real                            :: rlat1,rlon1,rlat2,rlon2,rlati
        real                            :: urot,vrot,uo2(mo,km),vo2(mo,km)
        real                            :: xmesh,xp,yp,xpts(mo),ypts(mo)

        type(grib1_descriptor) :: desc_in,desc_out
        class(ip_grid),allocatable :: grid_in,grid_out

        desc_in=init_descriptor(kgdsi)
        desc_out=init_descriptor(kgdso)

        call init_grid(grid_in,desc_in)
        call init_grid(grid_out,desc_out)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  COMPUTE NUMBER OF OUTPUT POINTS AND THEIR LATITUDES AND LONGITUDES.
        iret=0
        if(kgdso(1).ge.0) then
            call gdswzd(grid_out,0,mo,fill,xpts,ypts,rlon,rlat,no,crot,srot)
            if(no.eq.0) iret=3
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  AFFIRM APPROPRIATE INPUT GRID
        !    LAT/LON OR GAUSSIAN
        !    NO BITMAPS
        !    FULL ZONAL COVERAGE
        !    FULL MERIDIONAL COVERAGE
        idrti=kgdsi(1)
        im=kgdsi(2)
        jm=kgdsi(3)
        rlon1=kgdsi(5)*1.e-3
        rlon2=kgdsi(8)*1.e-3
        iscan=mod(kgdsi(11)/128,2)
        jscan=mod(kgdsi(11)/64,2)
        nscan=mod(kgdsi(11)/32,2)
        if(idrti.ne.0.and.idrti.ne.4) iret=41
        do k=1,km
            if(ibi(k).ne.0) iret=41
        enddo
        if(iret.eq.0) then
            if(iscan.eq.0) then
                dlon=(mod(rlon2-rlon1-1+3600,360.)+1)/(im-1)
            else
                dlon=-(mod(rlon1-rlon2-1+3600,360.)+1)/(im-1)
            endif
            ig=nint(360/abs(dlon))
            iprime=1+mod(-nint(rlon1/dlon)+ig,ig)
            imaxi=ig
            jmaxi=jm
            if(mod(ig,2).ne.0.or.im.lt.ig) iret=41
        endif
        if(iret.eq.0.and.idrti.eq.0) then
            rlat1=kgdsi(4)*1.e-3
            rlat2=kgdsi(7)*1.e-3
            dlat=(rlat2-rlat1)/(jm-1)
            jg=nint(180/abs(dlat))
            if(jm.eq.jg) idrti=256
            if(jm.ne.jg.and.jm.ne.jg+1) iret=41
        elseif(iret.eq.0.and.idrti.eq.4) then
            jg=kgdsi(10)*2
            if(jm.ne.jg) iret=41
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  SET PARAMETERS
        if(iret.eq.0) then
            iromb=ipopt(1)
            maxwv=ipopt(2)
            if(maxwv.eq.-1) then
                if(iromb.eq.0.and.idrti.eq.4) maxwv=(jmaxi-1)
                if(iromb.eq.1.and.idrti.eq.4) maxwv=(jmaxi-1)/2
                if(iromb.eq.0.and.idrti.eq.0) maxwv=(jmaxi-3)/2
                if(iromb.eq.1.and.idrti.eq.0) maxwv=(jmaxi-3)/4
                if(iromb.eq.0.and.idrti.eq.256) maxwv=(jmaxi-1)/2
                if(iromb.eq.1.and.idrti.eq.256) maxwv=(jmaxi-1)/4
            endif
            if((iromb.ne.0.and.iromb.ne.1).or.maxwv.lt.0) iret=42
        endif
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  INTERPOLATE
        if(iret.eq.0) then
            if(nscan.eq.0) then
                iskipi=1
                jskipi=im
            else
                iskipi=jm
                jskipi=1
            endif
            if(iscan.eq.1) iskipi=-iskipi
            if(jscan.eq.0) jskipi=-jskipi
            ispec=0
            !  SPECIAL CASE OF GLOBAL CYLINDRICAL GRID
            if((kgdso(1).eq.0.or.kgdso(1).eq.4).and. &
               mod(kgdso(2),2).eq.0.and.kgdso(5).eq.0.and. &
               kgdso(11).eq.0) then
                idrto=kgdso(1)
                imo=kgdso(2)
                jmo=kgdso(3)
                rlon2=kgdso(8)*1.e-3
                dlono=(mod(rlon2-1+3600,360.)+1)/(imo-1)
                igo=nint(360/abs(dlono))
                if(imo.eq.igo.and.idrto.eq.0) then
                    rlat1=kgdso(4)*1.e-3
                    rlat2=kgdso(7)*1.e-3
                    dlat=(rlat2-rlat1)/(jmo-1)
                    jgo=nint(180/abs(dlat))
                    if(jmo.eq.jgo) idrto=256
                    if(jmo.eq.jgo.or.jmo.eq.jgo+1) ispec=1
                elseif(imo.eq.igo.and.idrto.eq.4) then
                    jgo=kgdso(10)*2
                    if(jmo.eq.jgo) ispec=1
                endif
                if(ispec.eq.1) then
                    call sptrunv(iromb,maxwv,idrti,imaxi,jmaxi,idrto,imo,jmo, &
                                 km,iprime,iskipi,jskipi,mi,0,0,mo,0,ui,vi, &
                                 .true.,uo,vo,.false.,dumm,dumm,.false.,dumm,dumm)
                endif
                !  SPECIAL CASE OF POLAR STEREOGRAPHIC GRID
            elseif(kgdso(1).eq.5.and. &
                   kgdso(2).eq.kgdso(3).and.mod(kgdso(2),2).eq.1.and. &
                   kgdso(8).eq.kgdso(9).and.kgdso(11).eq.64.and. &
                   mod(kgdso(6)/8,2).eq.1) then
                nps=kgdso(2)
                rlat1=kgdso(4)*1.e-3
                rlon1=kgdso(5)*1.e-3
                orient=kgdso(7)*1.e-3
                xmesh=kgdso(8)
                iproj=mod(kgdso(10)/128,2)
                ip=(nps+1)/2
                h=(-1.)**iproj
                de=(1.+sin(60./dpr))*rerth
                dr=de*cos(rlat1/dpr)/(1+h*sin(rlat1/dpr))
                xp=1-h*sin((rlon1-orient)/dpr)*dr/xmesh
                yp=1+cos((rlon1-orient)/dpr)*dr/xmesh
                if(nint(xp).eq.ip.and.nint(yp).eq.ip) then
                    if(iproj.eq.0) then
                        call sptrunsv(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                      iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                      60.,xmesh,orient,ui,vi,.true.,uo,vo,uo2,vo2, &
                                      .false.,dumm,dumm,dumm,dumm, &
                                      .false.,dumm,dumm,dumm,dumm)
                    else
                        call sptrunsv(iromb,maxwv,idrti,imaxi,jmaxi,km,nps, &
                                      iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                      60.,xmesh,orient,ui,vi,.true.,uo2,vo2,uo,vo, &
                                      .false.,dumm,dumm,dumm,dumm, &
                                      .false.,dumm,dumm,dumm,dumm)
                    endif
                    ispec=1
                endif
                !  SPECIAL CASE OF MERCATOR GRID
            elseif(kgdso(1).eq.1) then
                ni=kgdso(2)
                nj=kgdso(3)
                rlat1=kgdso(4)*1.e-3
                rlon1=kgdso(5)*1.e-3
                rlon2=kgdso(8)*1.e-3
                rlati=kgdso(9)*1.e-3
                iscano=mod(kgdso(11)/128,2)
                jscano=mod(kgdso(11)/64,2)
                nscano=mod(kgdso(11)/32,2)
                dy=kgdso(13)
                hi=(-1.)**iscano
                hj=(-1.)**(1-jscano)
                dlono=hi*(mod(hi*(rlon2-rlon1)-1+3600,360.)+1)/(ni-1)
                dlato=hj*dy/(rerth*cos(rlati/dpr))*dpr
                if(nscano.eq.0) then
                    call sptrunmv(iromb,maxwv,idrti,imaxi,jmaxi,km,ni,nj, &
                                  iprime,iskipi,jskipi,mi,mo,0,0,0, &
                                  rlat1,rlon1,dlato,dlono,ui,vi, &
                                  .true.,uo,vo,.false.,dumm,dumm,.false.,dumm,dumm)
                    ispec=1
                endif
            endif
            !  GENERAL SLOW CASE
            if(ispec.eq.0) then
                call sptrungv(iromb,maxwv,idrti,imaxi,jmaxi,km,no, &
                              iprime,iskipi,jskipi,mi,mo,0,0,0,rlat,rlon, &
                              ui,vi,.true.,uo,vo,.false.,dumm,dumm,.false.,dumm,dumm)
                do k=1,km
                    ibo(k)=0
                    do n=1,no
                        lo(n,k)=.true.
                        urot=crot(n)*uo(n,k)-srot(n)*vo(n,k)
                        vrot=srot(n)*uo(n,k)+crot(n)*vo(n,k)
                        uo(n,k)=urot
                        vo(n,k)=vrot
                    enddo
                enddo
            endif
        else
            do k=1,km
                ibo(k)=1
                do n=1,no
                    lo(n,k)=.false.
                    uo(n,k)=0.
                    vo(n,k)=0.
                enddo
            enddo
        endif
    endsubroutine polatev4_grib1
endmodule spectral_interp_mod
