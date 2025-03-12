!> @file
!> @brief Rotated equidistant cylindrical GRIB decoder and grid
!> coordinate transformations for Arakawa grids A through D.
!>
!> @author Mark Iredell, George Gayno, Kyle Gerheiser
!> @date July 2021

!> Rotated equidistant cylindrical GRIB decoder and grid coordinate
!> transformations for Arakawa grids A through D. (To handle the E
!> grid, see ip_rot_equid_cylind_egrid_mod).
!>
!> See more info about [Awakawa
!> grids](https://en.wikipedia.org/wiki/Arakawa_grids).
!>
!> Octet numbers refer to [GRIB2 - GRID DEFINITION TEMPLATE 3.1 Rotate
!> Latitude/Longitude](https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_temp3-1.shtml).
!>
!> @author Gayno @date 2007-NOV-15
module ip_rot_equid_cylind_grid_mod
  use iso_fortran_env, only: real64
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use ip_constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_rot_equid_cylind_grid

  integer, parameter :: kd = real64 !< Fortran kind for reals.

  type, extends(ip_grid) :: ip_rot_equid_cylind_grid
     real(kd) :: clat0 !< Cosine of the latitude of the southern pole of projection.
     real(kd) :: dlats !< 'J'-direction grid increment.
     real(kd) :: dlons !< 'I'-direction grid increment.
     real(kd) :: rlon0 !< Longitude of southern pole of projection.
     real(kd) :: slat0 !< Sine of the latitude of the southern pole of projection.
     real(kd) :: wbd !< Longitude of the western boundary of the grid before rotation.
     real(kd) :: sbd !<  Latitude of the southern boundary of the grid before rotation.
     !> Rotation flag. When '0' the u/v vector components are relative
     !> to north/east. When '1' the u/v vector components are grid
     !> relative.
     integer :: irot
   contains
     !> Initializes a Rotated equidistant cylindrical grid given a
     !> grib1_descriptor object.
     !> @return N/A @copydoc ip_rot_equid_cylind_grid_mod::init_grib1
     procedure :: init_grib1
     !> Initializes a Rotated equidistant cylindrical given a
     !> grib2_descriptor object
     !> @return N/A @copydoc ip_rot_equid_cylind_grid_mod::init_grib2
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
     !> for Gaussian grids.
     !> @return N/A @copydoc ip_rot_equid_cylind_grid_mod::gdswzd_rot_equid_cylind
     procedure :: gdswzd => gdswzd_rot_equid_cylind
  end type ip_rot_equid_cylind_grid

  INTEGER :: IROT !< Local copy of irot.
  REAL(KIND=KD) :: RERTH !< Radius of the Earth.
  REAL(KIND=KD) :: CLAT0 !< Local copy of clat0.
  REAL(KIND=KD) :: DLATS !< Local copy of dlats.
  REAL(KIND=KD) :: DLONS !< Local copy of dlons.
  REAL(KIND=KD) :: RLON0 !< Local copy of rlon0.
  REAL(KIND=KD) :: SLAT0 !< Local copy of slat0.

CONTAINS

  !> Initializes a Rotated equidistant cylindrical grid given a
  !> grib1_descriptor object.
  !>
  !> @param[inout] self The grid to initialize
  !> @param[in] g1_desc A grib1_descriptor
  !>
  !> @author Gayno @date 2007-NOV-15
  subroutine init_grib1(self, g1_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    real(kd) :: hs, hs2, slat1, slat2, slatr, clon1, clon2, clat1, clat2, clatr, clonr, rlonr, rlatr

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6_KD
      self%eccen_squared = 0d0

      RLAT1=KGDS(4)*1.E-3_KD
      RLON1=KGDS(5)*1.E-3_KD
      RLAT0=KGDS(7)*1.E-3_KD
      self%RLON0=KGDS(8)*1.E-3_KD
      RLAT2=KGDS(12)*1.E-3_KD
      RLON2=KGDS(13)*1.E-3_KD

      self%IROT=MOD(KGDS(6)/8,2)
      self%IM=KGDS(2)
      self%JM=KGDS(3)

      SLAT1=SIN(RLAT1/DPR)
      CLAT1=COS(RLAT1/DPR)
      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      HS=SIGN(1._KD,MOD(RLON1-self%RLON0+180+3600,360._KD)-180)
      CLON1=COS((RLON1-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT1-self%SLAT0*CLAT1*CLON1
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT1*CLON1+self%SLAT0*SLAT1)/CLATR
      RLATR=DPR*ASIN(SLATR)
      RLONR=HS*DPR*ACOS(CLONR)

      self%WBD=RLONR
      self%SBD=RLATR
      SLAT2=SIN(RLAT2/DPR)
      CLAT2=COS(RLAT2/DPR)
      HS2=SIGN(1._KD,MOD(RLON2-self%RLON0+180+3600,360._KD)-180)
      CLON2=COS((RLON2-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT2-self%SLAT0*CLAT2*CLON2
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT2*CLON2+self%SLAT0*SLAT2)/CLATR
      NBD=DPR*ASIN(SLATR)
      EBD=HS2*DPR*ACOS(CLONR)
      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
    end associate

  end subroutine init_grib1

  !> Initializes a Rotated equidistant cylindrical grid given a
  !> grib2_descriptor object. Call 'use_ncep_post_arakawa()' before
  !> using this subroutine to use ncep_post-compatible grid definition.
  !>
  !> @param[inout] self The grid to initialize
  !> @param[in] g2_desc A grib2_descriptor
  !>
  !> @author Alex Richert @date 2024-MAY-20
  subroutine init_grib2(self, g2_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    if (ncep_post_arakawa.and.(g2_desc%gdt_num.eq.32769)) then
      call init_grib2_ncep_post(self, g2_desc)
    else
      call init_grib2_default(self, g2_desc)
    endif

  end subroutine init_grib2

  !> Initializes a Rotated equidistant cylindrical grid given a
  !> grib2_descriptor object. Uses template definitions consistent
  !> with WMO standards.
  !>
  !> @param[inout] self The grid to initialize
  !> @param[in] g2_desc A grib2_descriptor
  !>
  !> @author Gayno @date 2007-NOV-15
  subroutine init_grib2_default(self, g2_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    integer :: iscale
    integer :: i_offset_odd, i_offset_even, j_offset

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

      CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,self%rerth,self%eccen_squared)

      I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
      J_OFFSET=MOD(IGDTMPL(19)/2,2)

      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6

      RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
      RLAT0=RLAT0+90.0_KD

      self%RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)

      RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
      RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)

      self%IROT=MOD(IGDTMPL(14)/8,2)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      self%WBD=RLON1
      IF (self%WBD > 180.0) self%WBD = self%WBD - 360.0
      self%SBD=RLAT1

      NBD=RLAT2
      EBD=RLON2

      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)

      IF(I_OFFSET_ODD==1) self%WBD=self%WBD+(0.5_KD*self%DLONS)
      IF(J_OFFSET==1) self%SBD=self%SBD+(0.5_KD*self%DLATS)

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
      self%nscan = mod(igdtmpl(19) / 32, 2)
      self%nscan_field_pos = self%nscan
    end associate
  end subroutine init_grib2_default

  !> Initializes a Rotated equidistant cylindrical grid given a
  !> grib2_descriptor object. Uses template number definitions in
  !> a manner compatible with wgrib2 and ncep_post, which is based
  !> on grib1.
  !>
  !> @param[inout] self The grid to initialize
  !> @param[in] g2_desc A grib2_descriptor
  !>
  !> @author Alex Richert, George Gayno @date 2024-MAY-20
  subroutine init_grib2_ncep_post(self, g2_desc)
    class(ip_rot_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    real(kd) :: rlat1, rlon1, rlat0, rlat2, rlon2, nbd, ebd
    integer :: iscale
    integer :: i_offset_odd, i_offset_even, j_offset
    real(kd) :: hs, hs2, slat1, slat2, slatr, clon1, clon2, clat1, clat2, clatr, clonr, rlonr, rlatr

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)

      CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,self%rerth,self%eccen_squared)

      I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
      J_OFFSET=MOD(IGDTMPL(19)/2,2)

      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6

      RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      RLAT0=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)

      self%RLON0=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)

      RLAT2=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
      RLON2=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)

      self%IROT=MOD(IGDTMPL(14)/8,2)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      SLAT1=SIN(RLAT1/DPR)
      CLAT1=COS(RLAT1/DPR)

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      HS=SIGN(1._KD,MOD(RLON1-self%RLON0+180+3600,360._KD)-180)

      CLON1=COS((RLON1-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT1-self%SLAT0*CLAT1*CLON1
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT1*CLON1+self%SLAT0*SLAT1)/CLATR
      RLATR=DPR*ASIN(SLATR)
      RLONR=HS*DPR*ACOS(CLONR)

      self%WBD=RLONR
      IF (self%WBD > 180.0) self%WBD = self%WBD - 360.0
      self%SBD=RLATR

      SLAT2=SIN(RLAT2/DPR)
      CLAT2=COS(RLAT2/DPR)
      HS2=SIGN(1._KD,MOD(RLON2-self%RLON0+180+3600,360._KD)-180)
      CLON2=COS((RLON2-self%RLON0)/DPR)
      SLATR=self%CLAT0*SLAT2-self%SLAT0*CLAT2*CLON2
      CLATR=SQRT(1-SLATR**2)
      CLONR=(self%CLAT0*CLAT2*CLON2+self%SLAT0*SLAT2)/CLATR
      NBD=DPR*ASIN(SLATR)
      EBD=HS2*DPR*ACOS(CLONR)
      self%DLATS=(NBD-self%SBD)/FLOAT(self%JM-1)
      self%DLONS=(EBD-self%WBD)/FLOAT(self%IM-1)

      IF(I_OFFSET_ODD==1) self%WBD=self%WBD+(0.5_KD*self%DLONS)
      IF(J_OFFSET==1) self%SBD=self%SBD+(0.5_KD*self%DLATS)

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%kscan = 0
      self%nscan = mod(igdtmpl(19) / 32, 2)
      self%nscan_field_pos = self%nscan
    end associate
  end subroutine init_grib2_ncep_post

  !> GDS wizard for rotated equidistant cylindrical.
  !>
  !> This subprogram decodes the grib 2 grid definition template
  !> (passed in integer form as decoded by the ncep g2 library) and
  !> returns one of the following:
  !> - (iopt=+1) earth coordinates of selected grid coordinates
  !> - (iopt=-1) grid coordinates of selected earth coordinates
  !>
  !> Works for non-"e" staggered rotated equidistant cylindrical
  !> projections. the scan mode (section 3, octet 72, bits 5-6)
  !> determine whether this is an "h" or "v" grid.
  !>
  !> If the selected coordinates are more than one gridpoint beyond
  !> the the edges of the grid domain, then the relevant output
  !> elements are set to fill values. The actual number of valid
  !> points computed is returned too.
  !>
  !> Optionally, the vector rotations, the map jacobians and the grid
  !> box areas may be returned as well.
  !>
  !> To compute the vector rotations, the optional arguments 'srot'
  !> and 'crot' must be present. To compute the map jacobians, the
  !> optional arguments 'xlon', 'xlat', 'ylon', 'ylat' must be
  !> present. To compute the grid box areas, the optional argument
  !> 'area' must be present.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2010-jan-15 | gayno | based on routines gdswzdcb and gdswzdca
  !> 2015-jan-21 | gayno | merger of gdswizcd and gdswzdcd. make crot,sort,xlon,xlat,ylon,ylat and area optional arguments. make part of a module. move vector rotation, map jacobian and grid box area computations to separate subroutines.
  !> 2015-jul-13 | gayno | convert to grib 2. replace grib 1 kgds array with grib 2 grid definition template array. rename as "gdswzd_rot_equid_cylind."
  !> 2018-07-20 | wesley | add threads.
  !>
  !> @param[in] self Module reference.
  !> @param[in] iopt integer option flag
  !> - 1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param[in] npts integer maximum number of coordinates
  !> @param[in] fill real fill value to set invalid output data
  !> (must be impossible value; suggested value: -9999.)
  !> @param[inout] xpts real (npts) grid x point coordinates if iopt>0
  !> @param[inout] ypts real (npts) grid y point coordinates if iopt>0
  !> @param[inout] rlon real (npts) earth longitudes in degrees e if iopt<0
  !> (acceptable range: -360. to 360.)
  !> @param[inout] rlat real (npts) earth latitudes in degrees n if iopt<0
  !> (acceptable range: -90. to 90.)
  !> @param[out] nret integer number of valid points computed
  !> @param[out] crot real, optional (npts) clockwise vector rotation cosines
  !> @param[out] srot real, optional (npts) clockwise vector rotation sines
  !> (ugrid=crot*uearth-srot*vearth;
  !> vgrid=srot*uearth+crot*vearth)
  !> @param[out] xlon real, optional (npts) dx/dlon in 1/degrees
  !> @param[out] xlat real, optional (npts) dx/dlat in 1/degrees
  !> @param[out] ylon real, optional (npts) dy/dlon in 1/degrees
  !> @param[out] ylat real, optional (npts) dy/dlat in 1/degrees
  !> @param[out] area real, optional (npts) area weights in m**2
  !>
  !> @author Gayno @date 2007-NOV-15
  SUBROUTINE GDSWZD_ROT_EQUID_CYLIND(self,IOPT,NPTS, &
       FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE

    class(ip_rot_equid_cylind_grid), intent(in) :: self
    INTEGER,                 INTENT(IN   ) :: IOPT, NPTS
    INTEGER,                 INTENT(  OUT) :: NRET
    !
    REAL,                    INTENT(IN   ) :: FILL
    REAL,                    INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,                    INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL,  OPTIONAL,         INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                                :: IM,JM,N
    !
    LOGICAL                                :: LROT, LMAP, LAREA
    !
    REAL(KIND=KD)                          :: HS
    REAL(KIND=KD)                          :: CLONR,CLATR,SLATR
    REAL(KIND=KD)                          :: CLAT,SLAT,CLON
    REAL(KIND=KD)                          :: RLATR,RLONR
    REAL(KIND=KD)                          :: WBD,SBD
    REAL                                   :: XMIN,XMAX,YMIN,YMAX
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IS THE EARTH RADIUS DEFINED?
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  IS THIS AN "E"-STAGGER GRID?  ROUTINE CAN'T PROCESS THOSE.
    ! I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
    ! I_OFFSET_EVEN=MOD(IGDTMPL(19)/4,2)
    ! J_OFFSET=MOD(IGDTMPL(19)/2,2)
    ! IF(I_OFFSET_ODD/=I_OFFSET_EVEN) THEN
    !    CALL ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !    RETURN
    ! ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    RLON0=self%rlon0
    IROT=self%irot

    IM=self%im
    JM=self%jm

    SLAT0=self%slat0
    CLAT0=self%clat0

    WBD=self%wbd
    SBD=self%sbd

    DLATS=self%dlats
    DLONS=self%dlons

    XMIN=0
    XMAX=IM+1
    YMIN=0
    YMAX=JM+1
    NRET=0

    rerth = self%rerth
    IF(RERTH<0.)THEN
       CALL ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
       RETURN
    ENDIF

    IF(PRESENT(CROT).AND.PRESENT(SROT))THEN
       LROT=.TRUE.
    ELSE
       LROT=.FALSE.
    ENDIF
    IF(PRESENT(XLON).AND.PRESENT(XLAT).AND.PRESENT(YLON).AND.PRESENT(YLAT))THEN
       LMAP=.TRUE.
    ELSE
       LMAP=.FALSE.
    ENDIF
    IF(PRESENT(AREA))THEN
       LAREA=.TRUE.
    ELSE
       LAREA=.FALSE.
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       !$OMP PARALLEL DO PRIVATE(N,RLONR,RLATR,HS,CLONR,SLATR,CLATR,SLAT,CLAT,CLON) &
       !$OMP& REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLONR=WBD+(XPTS(N)-1._KD)*DLONS
             RLATR=SBD+(YPTS(N)-1._KD)*DLATS
             IF(RLONR <= 0._KD) THEN
                HS=-1.0_KD
             ELSE
                HS=1.0_KD
             ENDIF
             CLONR=COS(RLONR/DPR)
             SLATR=SIN(RLATR/DPR)
             CLATR=COS(RLATR/DPR)
             SLAT=CLAT0*SLATR+SLAT0*CLATR*CLONR
             IF(SLAT.LE.-1) THEN
                CLAT=0.
                CLON=COS(RLON0/DPR)
                RLON(N)=0.
                RLAT(N)=-90.
             ELSEIF(SLAT.GE.1) THEN
                CLAT=0.
                CLON=COS(RLON0/DPR)
                RLON(N)=0.
                RLAT(N)=90.
             ELSE
                CLAT=SQRT(1-SLAT**2)
                CLON=(CLAT0*CLATR*CLONR-SLAT0*SLATR)/CLAT
                CLON=MIN(MAX(CLON,-1._KD),1._KD)
                RLON(N)=REAL(MOD(RLON0+HS*DPR*ACOS(CLON)+3600,360._KD))
                RLAT(N)=REAL(DPR*ASIN(SLAT))
             ENDIF
             NRET=NRET+1
             IF(LROT) CALL ROT_EQUID_CYLIND_VECT_ROT(RLON(N), CLATR, SLATR, &
                  CLAT, SLAT, CLON, CROT(N), SROT(N))
             IF(LMAP) CALL ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON(N), CLATR, &
                  CLAT, SLAT, CLON, XLON(N), XLAT(N), YLON(N), YLAT(N))
             IF(LAREA) CALL ROT_EQUID_CYLIND_GRID_AREA(CLATR, FILL, AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N,HS,CLON,SLAT,CLAT,SLATR,CLATR,CLONR,RLONR,RLATR) &
       !$OMP& REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             HS=SIGN(1._KD,MOD(RLON(N)-RLON0+180+3600,360._KD)-180)
             CLON=COS((RLON(N)-RLON0)/DPR)
             SLAT=SIN(RLAT(N)/DPR)
             CLAT=COS(RLAT(N)/DPR)
             SLATR=CLAT0*SLAT-SLAT0*CLAT*CLON
             IF(SLATR.LE.-1) THEN
                CLATR=0._KD
                RLONR=0.
                RLATR=-90.
             ELSEIF(SLATR.GE.1) THEN
                CLATR=0._KD
                RLONR=0.
                RLATR=90.
             ELSE
                CLATR=SQRT(1-SLATR**2)
                CLONR=(CLAT0*CLAT*CLON+SLAT0*SLAT)/CLATR
                CLONR=MIN(MAX(CLONR,-1._KD),1._KD)
                RLONR=HS*DPR*ACOS(CLONR)
                RLATR=DPR*ASIN(SLATR)
             ENDIF
             XPTS(N)=REAL((RLONR-WBD)/DLONS+1._KD)
             YPTS(N)=REAL((RLATR-SBD)/DLATS+1._KD)
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT) CALL ROT_EQUID_CYLIND_VECT_ROT(RLON(N), CLATR, SLATR, &
                     CLAT, SLAT, CLON, CROT(N), SROT(N))
                IF(LMAP) CALL ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON(N), CLATR, &
                     CLAT, SLAT, CLON, XLON(N), XLAT(N), YLON(N), YLAT(N))
                IF(LAREA) CALL ROT_EQUID_CYLIND_GRID_AREA(CLATR, FILL, AREA(N))
             ELSE
                XPTS(N)=FILL
                YPTS(N)=FILL
             ENDIF
          ELSE
             XPTS(N)=FILL
             YPTS(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_ROT_EQUID_CYLIND

  !> Error handler.
  !>
  !> Upon an error, this subprogram assigns a "fill" value to the
  !> output fields.
  !>
  !> @param[in] iopt integer option flag
  !> - +1 to compute earth coords of selected grid coords
  !> - -1 to compute grid coords of selected earth coords
  !> @param[in] fill real fill value to set invalid output data
  !> (must be impossible value; suggested value: -9999.)
  !> @param[out] rlat real (npts) earth latitudes in degrees n if iopt<0
  !> @param[out] rlon real (npts) earth longitudes in degrees e if iopt<0
  !> @param[out] xpts real (npts) grid x point coordinates if iopt>0
  !> @param[out] ypts real (npts) grid y point coordinates if iopt>0
  !> @param[in] npts integer maximum number of coordinates
  !>
  !> @author Gayno @date 2015-07-13
  SUBROUTINE ROT_EQUID_CYLIND_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: IOPT, NPTS
    !
    REAL,    INTENT(IN   ) :: FILL
    REAL,    INTENT(  OUT) :: RLAT(NPTS),RLON(NPTS)
    REAL,    INTENT(  OUT) :: XPTS(NPTS),YPTS(NPTS)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(IOPT>=0) THEN
       RLON=FILL
       RLAT=FILL
    ENDIF
    IF(IOPT<=0) THEN
       XPTS=FILL
       YPTS=FILL
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE ROT_EQUID_CYLIND_ERROR

  !> Vector rotation fields for rotated equidistant cylindrical grids -
  !> non "e" stagger.
  !>
  !> This subprogram computes the vector rotation sines and cosines
  !> for a rotated equidistant cylindrical grid - non "e" stagger.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-07-19 | gayno | rename as "rot_equid_cylind_vect_rot."
  !> 2018-07-20 | wesley | pass in clatr, slatr, clat, slat, clon for threading.
  !>
  !> @param[in] rlon longitude in degrees (real)
  !> @param[in] clatr cosine of rotated latitude (real)
  !> @param[in] slatr sine of rotated latitude (real)
  !> @param[in] clat cosine of latitude (real)
  !> @param[in] slat sine of latitude (real)
  !> @param[in] clon cosine of longitude (real)
  !> @param[out] crot clockwise vector rotation cosines (real)
  !> @param[out] srot clockwise vector rotation sines (real)
  !> (ugrid=crot*uearth-srot*vearth;
  !> vgrid=srot*uearth+crot*vearth)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE ROT_EQUID_CYLIND_VECT_ROT(RLON, CLATR, SLATR, CLAT, SLAT, &
       CLON, CROT, SROT)
    IMPLICIT NONE

    REAL(KIND=KD),    INTENT(IN   ) :: CLAT, CLATR, CLON, SLAT, SLATR
    REAL         ,    INTENT(IN   ) :: RLON
    REAL         ,    INTENT(  OUT) :: CROT, SROT

    REAL(KIND=KD)                   :: SLON

    IF(IROT.EQ.1) THEN
       IF(CLATR.LE.0._KD) THEN
          CROT=REAL(-SIGN(1._KD,SLATR*SLAT0))
          SROT=0.
       ELSE
          SLON=SIN((RLON-RLON0)/DPR)
          CROT=REAL((CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR)
          SROT=REAL(SLAT0*SLON/CLATR)
       ENDIF
    ELSE
       CROT=1.
       SROT=0.
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_VECT_ROT
  
  !> Map jacobians for rotated equidistant cylindrical
  !> grids - non "e" stagger.
  !>
  !> This subprogram computes the map jacobians for a rotated
  !> equidistant cylindrical grid - non "e" stagger.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | gayno | initial version
  !> 2015-09-17 | gayno | rename as "rot_equid_cylind_map_jacob".
  !> 2018-07-20 | wesley | pass in clatr, clat, slat, clon to allow threading.
  !>
  !> @param[in] fill fill value for undefined points (real)
  !> @param[in] rlon longitude in degrees (real)
  !> @param[in] clatr cosine of unrotated latitude (real)
  !> @param[in] clat cosine of latitude (real)
  !> @param[in] slat sine of latitude (real)
  !> @param[in] clon cosine of latitude (real)
  !> @param[out] xlon dx/dlon in 1/degrees (real)
  !> @param[out] xlat dx/dlat in 1/degrees (real)
  !> @param[out] ylon dy/dlon in 1/degrees (real)
  !> @param[out] ylat dy/dlat in 1/degrees (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE ROT_EQUID_CYLIND_MAP_JACOB(FILL, RLON, CLATR, CLAT, &
       SLAT, CLON, XLON, XLAT, YLON, YLAT)
    IMPLICIT NONE

    REAL(KIND=KD),    INTENT(IN   ) :: CLATR, CLAT, SLAT, CLON
    REAL         ,    INTENT(IN   ) :: FILL, RLON
    REAL         ,    INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    REAL(KIND=KD)                   :: SLON, TERM1, TERM2

    IF(CLATR.LE.0._KD) THEN
       XLON=FILL
       XLAT=FILL
       YLON=FILL
       YLAT=FILL
    ELSE
       SLON=SIN((RLON-RLON0)/DPR)
       TERM1=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
       TERM2=SLAT0*SLON/CLATR
       XLON=REAL(TERM1*CLAT/(DLONS*CLATR))
       XLAT=REAL(-TERM2/(DLONS*CLATR))
       YLON=REAL(TERM2*CLAT/DLATS)
       YLAT=REAL(TERM1/DLATS)
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_MAP_JACOB

  !> Grid box area for rotated equidistant cylindrical grids - non "e"
  !> stagger.
  !>
  !> This subprogram computes the grid box area for a rotated
  !> equidistant cylindrical grid - non "e" stagger.
  !>
  !> ### Program History Log
  !> Date | Programmer | Comments
  !> -----|------------|---------
  !> 2015-01-21 | Gayno | initial version
  !> 2015-07-19 | gayno | rename as "rot_equid_cylind_grid_area."
  !> 2018-07-20 | wesley | pass in clatr for threading
  !>
  !> @param[in] clatr cosine of unrotated latitude (real)
  !> @param[in] fill fill value for undefined points (real)
  !> @param[out] area area weights in m**2 (real)
  !>
  !> @author Gayno @date 2015-01-21
  SUBROUTINE ROT_EQUID_CYLIND_GRID_AREA(CLATR, FILL, AREA)
    IMPLICIT NONE

    REAL(KIND=KD),    INTENT(IN   ) :: CLATR
    REAL,             INTENT(IN   ) :: FILL
    REAL,             INTENT(  OUT) :: AREA

    IF(CLATR.LE.0._KD) THEN
       AREA=FILL
    ELSE
       AREA=REAL(2._KD*(RERTH**2)*CLATR*(DLONS/DPR)*SIN(0.5_KD*DLATS/DPR))
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_GRID_AREA

end module ip_rot_equid_cylind_grid_mod

