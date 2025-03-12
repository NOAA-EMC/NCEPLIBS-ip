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
  use iso_fortran_env, only: real64
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use ip_constants_mod, only: DPR, PI
  use earth_radius_mod
  implicit none

  private
  public :: ip_rot_equid_cylind_egrid

  integer, parameter :: kd = real64 !< Kind of reals.

  type, extends(ip_grid) :: ip_rot_equid_cylind_egrid
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
     !> grib1_descriptor object.
     !> @return N/A @copydoc ip_rot_equid_cylind_egrid_mod::init_grib1
     procedure :: init_grib1
     !> Initializes a rotated equidistant cylindrical grid given a
     !> grib2_descriptor object.
     !> @return N/A @copydoc ip_rot_equid_cylind_egrid_mod::init_grib2
     procedure :: init_grib2
     !> Calculates Earth coordinates (iopt = 1) or grid coorindates
     !> (iopt = -1).
     !> @return N/A @copydoc ip_rot_equid_cylind_egrid_mod::gdswzd_rot_equid_cylind_egrid
     procedure :: gdswzd => gdswzd_rot_equid_cylind_egrid
  end type ip_rot_equid_cylind_egrid

  INTEGER :: IROT !< Local copy of irot.

  REAL(KIND=KD) :: CLAT !< Cosine of the latitude.
  REAL(KIND=KD) :: CLAT0 !< Local copy of clat0.
  REAL(KIND=KD) :: CLATR !< Cosine of the rotated latitude.
  REAL(KIND=KD) :: CLON !< Cosine of the difference between rlon and rlon0.
  REAL(KIND=KD) :: DLATS !< Local copy of dlats.
  REAL(KIND=KD) :: DLONS !< Local copy of dlons.
  REAL(KIND=KD) :: RERTH !< Radius of the Earth.
  REAL(KIND=KD) :: RLON0 !< Local copy of rlon0.
  REAL(KIND=KD) :: SLAT !< Sine of the latitude.
  REAL(KIND=KD) :: SLAT0 !< Local copy of slat0.
  REAL(KIND=KD) :: SLATR !< Sine of the rotated latitude.

contains

  !> Initializes a rotated equidistant cylindrical grid given a
  !> grib1_descriptor object.
  !> 
  !> @param[inout] self The grid to initialize
  !> @param[in] g1_desc A grib1_descriptor
  !>
  !> @author Kyle Gerheiser
  !> @date July 2021
  subroutine init_grib1(self, g1_desc)
    class(ip_rot_equid_cylind_egrid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan
    real(kd) :: rlat0

    real(kd) :: rlat1, rlon1, rlon0, slat1, clat1, slat0, clat0, clon1
    real(kd) :: slatr, clatr, clonr, rlatr, rlonr, dlats, dlons, hs, hi
    integer :: im, jm

    integer :: is1, kscan,  irot

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6_KD
      self%eccen_squared = 0.0

      IM=KGDS(2)
      JM=KGDS(3)

      self%nscan_field_pos = 3
      self%nscan = MOD(KGDS(11)/32,2)

      RLAT1=KGDS(4)*1.E-3_KD
      RLON1=KGDS(5)*1.E-3_KD
      RLAT0=KGDS(7)*1.E-3_KD
      RLON0=KGDS(8)*1.E-3_KD

      IROT=MOD(KGDS(6)/8,2)
      KSCAN=MOD(KGDS(11)/256,2)
      ISCAN=MOD(KGDS(11)/128,2)
      HI=(-1.)**ISCAN
      SLAT1=SIN(RLAT1/DPR)
      CLAT1=COS(RLAT1/DPR)
      SLAT0=SIN(RLAT0/DPR)
      CLAT0=COS(RLAT0/DPR)
      HS=SIGN(1._KD,MOD(RLON1-RLON0+180+3600,360._KD)-180)
      CLON1=COS((RLON1-RLON0)/DPR)
      SLATR=CLAT0*SLAT1-SLAT0*CLAT1*CLON1
      CLATR=SQRT(1-SLATR**2)
      CLONR=(CLAT0*CLAT1*CLON1+SLAT0*SLAT1)/CLATR
      RLATR=DPR*ASIN(SLATR)
      RLONR=HS*DPR*ACOS(CLONR)
      DLATS=RLATR/(-(JM-1)/2)
      DLONS=RLONR/(-((IM * 2 - 1) -1)/2)

      IF(KSCAN.EQ.0) THEN
         IS1=(JM+1)/2
      ELSE
         IS1=JM/2
      ENDIF

      self%im = im
      self%jm = jm
      self%rlon0 = rlon0
      self%rlon1 = rlon1
      self%rlat1 = rlat1
      self%clat0 = clat0
      self%slat0 = slat0
      self%dlats = dlats
      self%dlons = dlons
      self%hi = hi
      self%irot = irot
      self%kscan = kscan


    end associate

  end subroutine init_grib1


  !> Initializes a rotated equidistant cylindrical grid given a grib2_descriptor object.
  !> @param[inout] self The grid to initialize
  !> @param[in] g2_desc A grib2_descriptor
  !>
  !> @author Kyle Gerheiser
  !> @date July 2021
  subroutine init_grib2(self, g2_desc)
    class(ip_rot_equid_cylind_egrid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan
    real(kd) :: rlat0
    integer :: i_offset_odd!, i_offset_even

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      CALL EARTH_RADIUS(IGDTMPL,IGDTLEN,self%rerth,self%eccen_squared)

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

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)

      self%NSCAN=MOD(IGDTMPL(16)/32,2)
      self%nscan_field_pos = 3

      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6

      self%RLON0=FLOAT(IGDTMPL(21))/FLOAT(ISCALE)
      self%DLATS=FLOAT(IGDTMPL(18))/FLOAT(ISCALE)
      ! THE GRIB2 CONVENTION FOR "I" RESOLUTION IS TWICE WHAT THIS ROUTINE ASSUMES.
      self%DLONS=FLOAT(IGDTMPL(17))/FLOAT(ISCALE) * 0.5_KD

      self%IROT=MOD(IGDTMPL(14)/8,2)

      I_OFFSET_ODD=MOD(IGDTMPL(19)/8,2)
      self%KSCAN=I_OFFSET_ODD
      ISCAN=MOD(IGDTMPL(19)/128,2)

      self%HI=(-1.)**ISCAN

      RLAT0=FLOAT(IGDTMPL(20))/FLOAT(ISCALE)
      RLAT0=RLAT0+90.0_KD

      self%SLAT0=SIN(RLAT0/DPR)
      self%CLAT0=COS(RLAT0/DPR)

      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
    end associate
  end subroutine init_grib2


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
  SUBROUTINE GDSWZD_ROT_EQUID_CYLIND_EGRID(self,IOPT,NPTS,&
       FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE
    !
    class(ip_rot_equid_cylind_egrid), intent(in) :: self

    INTEGER,         INTENT(IN   ) :: IOPT, NPTS
    INTEGER,         INTENT(  OUT) :: NRET
    !
    REAL,            INTENT(IN   ) :: FILL
    REAL,            INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,            INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,  INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,  INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,  INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                        :: IM, JM, IS1, N
    INTEGER                        :: KSCAN
    !    INTEGER                        :: I_OFFSET_ODD, I_OFFSET_EVEN
    !
    LOGICAL                        :: LROT, LMAP, LAREA
    !
    REAL(KIND=KD)                  :: RLAT1, RLON1
    REAL(KIND=KD)                  :: CLONR
    REAL(KIND=KD)                  :: RLATR, RLONR, SBD, WBD, HS, HI
    REAL                           :: XMAX, XMIN, YMAX, YMIN, XPTF, YPTF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL

    RLON0=self%rlon0
    IROT=self%irot
    IM=self%im * 2 - 1
    JM=self%jm
    DLATS=self%dlats
    DLONS=self%dlons
    KSCAN=self%kscan
    HI=self%hi
    SLAT0=self%slat0
    CLAT0=self%clat0
    RLAT1=self%rlat1
    RLON1=self%rlon1

    rerth = self%rerth

    ! IS THE EARTH RADIUS DEFINED?
    IF(RERTH<0.)THEN
       CALL ROT_EQUID_CYLIND_EGRID_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
       RETURN
    ENDIF

    SBD=RLAT1
    WBD=RLON1

    IF (WBD > 180.0) WBD = WBD - 360.0
    IF(KSCAN.EQ.0) THEN
       IS1=(JM+1)/2
    ELSE
       IS1=JM/2
    ENDIF

    XMIN=0
    XMAX=IM+2
    YMIN=0
    YMAX=JM+1
    NRET=0

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
       DO N=1,NPTS
          XPTF=YPTS(N)+(XPTS(N)-IS1)
          YPTF=YPTS(N)-(XPTS(N)-IS1)+KSCAN
          IF(XPTF.GE.XMIN.AND.XPTF.LE.XMAX.AND. &
               YPTF.GE.YMIN.AND.YPTF.LE.YMAX) THEN
             HS=HI*SIGN(1.,XPTF-(IM+1)/2)
             select type(desc => self%descriptor)
             type is(grib1_descriptor)
                RLONR=(XPTF-(IM+1)/2)*DLONS
                RLATR=(YPTF-(JM+1)/2)*DLATS
             type is(grib2_descriptor)
                RLONR=(XPTF-1.0_KD)*DLONS + WBD
                RLATR=(YPTF-1.0_KD)*DLATS + SBD 
             end select
             CLONR=COS(RLONR/DPR)
             SLATR=SIN(RLATR/DPR)
             CLATR=COS(RLATR/DPR)
             SLAT=CLAT0*SLATR+SLAT0*CLATR*CLONR
             IF(SLAT.LE.-1) THEN
                CLAT=0.
                CLON=COS(RLON0/DPR)
                RLON(N)=0
                RLAT(N)=-90
             ELSEIF(SLAT.GE.1) THEN
                CLAT=0.
                CLON=COS(RLON0/DPR)
                RLON(N)=0
                RLAT(N)=90
             ELSE
                CLAT=SQRT(1-SLAT**2)
                CLON=(CLAT0*CLATR*CLONR-SLAT0*SLATR)/CLAT
                CLON=MIN(MAX(CLON,-1._KD),1._KD)
                RLON(N)=REAL(MOD(RLON0+HS*DPR*ACOS(CLON)+3600,360._KD))
                RLAT(N)=REAL(DPR*ASIN(SLAT))
             ENDIF
             NRET=NRET+1
             IF(LROT) CALL ROT_EQUID_CYLIND_EGRID_VECT_ROT(RLON(N),CROT(N),SROT(N))
             IF(LMAP) CALL ROT_EQUID_CYLIND_EGRID_MAP_JACOB(FILL, RLON(N), &
                  XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL ROT_EQUID_CYLIND_EGRID_GRID_AREA(FILL, AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       DO N=1,NPTS
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             HS=SIGN(1._KD,MOD(RLON(N)-RLON0+180+3600,360._KD)-180)
             CLON=COS((RLON(N)-RLON0)/DPR)
             SLAT=SIN(RLAT(N)/DPR)
             CLAT=COS(RLAT(N)/DPR)
             SLATR=CLAT0*SLAT-SLAT0*CLAT*CLON
             IF(SLATR.LE.-1) THEN
                CLATR=0.
                RLONR=0
                RLATR=-90
             ELSEIF(SLATR.GE.1) THEN
                CLATR=0.
                RLONR=0
                RLATR=90
             ELSE
                CLATR=SQRT(1-SLATR**2)
                CLONR=(CLAT0*CLAT*CLON+SLAT0*SLAT)/CLATR
                CLONR=MIN(MAX(CLONR,-1._KD),1._KD)
                RLONR=HS*DPR*ACOS(CLONR)
                RLATR=DPR*ASIN(SLATR)
             ENDIF
             select type(desc => self%descriptor)
             type is(grib1_descriptor)
                XPTF=REAL(((RLONR-WBD)/DLONS)+1.0_KD)
                YPTF=REAL(((RLATR-SBD)/DLATS)+1.0_KD)
             type is(grib2_descriptor)
                XPTF=REAL((IM+1)/2+RLONR/DLONS)
                YPTF=REAL((JM+1)/2+RLATR/DLATS)
             end select

             IF(XPTF.GE.XMIN.AND.XPTF.LE.XMAX.AND. &
                  YPTF.GE.YMIN.AND.YPTF.LE.YMAX) THEN
                XPTS(N)=IS1+(XPTF-(YPTF-KSCAN))/2
                YPTS(N)=(XPTF+(YPTF-KSCAN))/2
                NRET=NRET+1
                IF(LROT) CALL ROT_EQUID_CYLIND_EGRID_VECT_ROT(RLON(N),CROT(N),SROT(N))
                IF(LMAP) CALL ROT_EQUID_CYLIND_EGRID_MAP_JACOB(FILL, RLON(N), &
                     XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL ROT_EQUID_CYLIND_EGRID_GRID_AREA(FILL, AREA(N))
             ELSE
                XPTS(N)=FILL
                YPTS(N)=FILL
             ENDIF
          ELSE
             XPTS(N)=FILL
             YPTS(N)=FILL
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE GDSWZD_ROT_EQUID_CYLIND_EGRID

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
  SUBROUTINE ROT_EQUID_CYLIND_EGRID_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN   ) :: IOPT, NPTS
    !
    REAL,    INTENT(IN   ) :: FILL
    REAL,    INTENT(  OUT) :: RLAT(NPTS),RLON(NPTS)
    REAL,    INTENT(  OUT) :: XPTS(NPTS),YPTS(NPTS)

    IF(IOPT>=0) THEN
       RLON=FILL
       RLAT=FILL
    ENDIF
    IF(IOPT<=0) THEN
       XPTS=FILL
       YPTS=FILL
    ENDIF
  END SUBROUTINE ROT_EQUID_CYLIND_EGRID_ERROR

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
  SUBROUTINE ROT_EQUID_CYLIND_EGRID_VECT_ROT(RLON, CROT, SROT)
    IMPLICIT NONE

    REAL         ,    INTENT(IN   ) :: RLON
    REAL         ,    INTENT(  OUT) :: CROT, SROT

    REAL(KIND=KD)                   :: SLON

    IF(IROT.EQ.1) THEN
       IF(CLATR.LE.0) THEN
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

  END SUBROUTINE ROT_EQUID_CYLIND_EGRID_VECT_ROT

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
  SUBROUTINE ROT_EQUID_CYLIND_EGRID_MAP_JACOB(FILL, RLON, &
       XLON, XLAT, YLON, YLAT)
    IMPLICIT NONE

    REAL         ,    INTENT(IN   ) :: FILL, RLON
    REAL         ,    INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    REAL(KIND=KD)                   :: SLON, TERM1, TERM2
    REAL(KIND=KD)                   :: XLATF, XLONF, YLATF, YLONF

    IF(CLATR.LE.0._KD) THEN
       XLON=FILL
       XLAT=FILL
       YLON=FILL
       YLAT=FILL
    ELSE
       SLON=SIN((RLON-RLON0)/DPR)
       TERM1=(CLAT0*CLAT+SLAT0*SLAT*CLON)/CLATR
       TERM2=SLAT0*SLON/CLATR
       XLONF=TERM1*CLAT/(DLONS*CLATR)
       XLATF=-TERM2/(DLONS*CLATR)
       YLONF=TERM2*CLAT/DLATS
       YLATF=TERM1/DLATS
       XLON=REAL(XLONF-YLONF)
       XLAT=REAL(XLATF-YLATF)
       YLON=REAL(XLONF+YLONF)
       YLAT=REAL(XLATF+YLATF)
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_EGRID_MAP_JACOB

  !> Computes the grid box area for a rotated equidistant cylindrical grid.
  !>
  !> @param[in] fill Fill value for undefined points.
  !> @param[out] area Area weights in m^2.
  !>
  !> @author George Gayno
  !> @date Jan 2015
  SUBROUTINE ROT_EQUID_CYLIND_EGRID_GRID_AREA(FILL, AREA)
    IMPLICIT NONE

    REAL,             INTENT(IN   ) :: FILL
    REAL,             INTENT(  OUT) :: AREA

    IF(CLATR.LE.0._KD) THEN
       AREA=FILL
    ELSE
       AREA=REAL(RERTH**2*CLATR*DLATS*DLONS)*2/DPR**2
    ENDIF

  END SUBROUTINE ROT_EQUID_CYLIND_EGRID_GRID_AREA

end module ip_rot_equid_cylind_egrid_mod


