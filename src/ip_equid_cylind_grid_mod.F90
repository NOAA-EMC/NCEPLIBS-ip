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

  type, extends(ip_grid) :: ip_equid_cylind_grid
     real :: hi !< Scan mode in the 'i' direction. GRIB2, Section 3, octet 72.
     real :: rlat1 !< Latitude of first grid point. GRIB2, Section 3, octets 47-50.
     real :: rlon1 !< Longitude of first grid point. GRIB2, Section 3, octets 51-54.
     real :: rlat2 !< Latitude of last grid point. GRIB2, Section 3, octets 56-59.
     real :: rlon2 !< Longitude of last grid point. GRIB2, Section 3, octets 60-63.
     real :: dlat !< Di — i direction increment. GRIB2, Section 3, octets 64-67.
     real :: dlon !< Dj — j direction increment. GRIB2, Section 3, octets 68-71.
   contains
     !> @return N/A @copydoc ip_equid_cylind_grid_mod::init_grib1
     procedure :: init_grib1 !< Init GRIB1.
     !> @return N/A @copydoc ip_equid_cylind_grid_mod::init_grib2
     procedure :: init_grib2 !< Init GRIB2.
     !> @return N/A @copydoc ip_equid_cylind_grid_mod::gdswzd_equid_cylind
     procedure :: gdswzd => gdswzd_equid_cylind
  end type ip_equid_cylind_grid

  REAL :: DLAT !< Grid resolution in degrees n/s direction.
  REAL :: DLON !< Grid resolution in degrees e/w direction.
  REAL :: RERTH !< Radius of the Earth.

contains

  !> Initializes an equidistant cylindrical grid given a grib1_descriptor object.
  !! 
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  subroutine init_grib1(self, g1_desc)
    class(ip_equid_cylind_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan

    associate(kgds => g1_desc%gds)
      self%IM=KGDS(2)
      self%JM=KGDS(3)
      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLAT2=KGDS(7)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      ISCAN=MOD(KGDS(11)/128,2)
      self%HI=(-1.)**ISCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DLAT=(self%RLAT2-self%RLAT1)/(self%JM-1)

      ! defaults
      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0

      self%iwrap = nint(360/abs(self%dlon))

      if(self%im < self%iwrap) self%iwrap=0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap,2) == 0) then
         if(abs(self%rlat1) > 90-0.25*self%dlat) then
            self%jwrap1 = 2
         elseif(abs(self%rlat1) > 90-0.75*self%dlat) then
            self%jwrap1 = 1
         endif
         if(abs(self%rlat2) > 90-0.25*self%dlat) then
            self%jwrap2 = 2 * self%jm
         elseif(abs(self%rlat2) > 90-0.75*self%dlat) then
            self%jwrap2 = 2 * self%jm+1
         endif
      endif

      self%rerth = 6.3712E6
      self%eccen_squared = 0.0
    end associate

  end subroutine init_grib1
  
  !> Initializes an equidistant cylindrical grid given a grib2_descriptor object.
  !! @param[inout] self The grid to initialize
  !! @param[in] g2_desc A grib2_descriptor
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  subroutine init_grib2(self, g2_desc)
    class(ip_equid_cylind_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)
      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6
      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      self%RLAT2=FLOAT(IGDTMPL(15))/FLOAT(ISCALE)
      self%RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
      ISCAN=MOD(IGDTMPL(19)/128,2)
      self%HI=(-1.)**ISCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)
      self%DLAT=(self%RLAT2-self%RLAT1)/(self%JM-1)

      self%nscan = MOD(IGDTMPL(19)/32,2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
      self%iwrap = NINT(360/ABS(self%DLON))

      if(self%im.lt.self%iwrap) self%iwrap=0
      self%jwrap1=0
      self%jwrap2=0

      if(self%im < self%iwrap) self%iwrap=0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap,2) == 0) then
         if(abs(self%rlat1) > 90-0.25*self%dlat) then
            self%jwrap1 = 2
         elseif(abs(self%rlat1) > 90-0.75*self%dlat) then
            self%jwrap1 = 1
         endif
         if(abs(self%rlat2) > 90-0.25*self%dlat) then
            self%jwrap2 = 2 * self%jm
         elseif(abs(self%rlat2) > 90-0.75*self%dlat) then
            self%jwrap2 = 2 * self%jm+1
         endif
      endif

      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

    end associate
  end subroutine init_grib2

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
  SUBROUTINE GDSWZD_EQUID_CYLIND(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET,  &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE
    !
    class(ip_equid_cylind_grid), intent(in) :: self
    INTEGER,             INTENT(IN   ) :: IOPT, NPTS
    INTEGER,             INTENT(  OUT) :: NRET
    !
    REAL,                INTENT(IN   ) :: FILL
    REAL,                INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,                INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL,      INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL,      INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL,      INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)
    !
    INTEGER                            :: IM, JM, N
    !
    LOGICAL                            :: LROT, LMAP, LAREA
    !
    REAL                               :: HI, RLAT1, RLON1, RLAT2, RLON2
    REAL                               :: XMAX, XMIN, YMAX, YMIN
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL

    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1
    RLAT2=self%rlat2
    RLON2=self%rlon2

    HI=self%hi

    rerth = self%rerth
    dlat = self%dlat
    dlon = self%dlon

    XMIN=0
    XMAX=IM+1
    IF(IM.EQ.NINT(360/ABS(DLON))) XMAX=IM+2
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
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLON(N)=MOD(RLON1+DLON*(XPTS(N)-1)+3600,360.)
             RLAT(N)=MIN(MAX(RLAT1+DLAT*(YPTS(N)-1),-90.),90.)
             NRET=NRET+1
             IF(LROT)  CALL EQUID_CYLIND_VECT_ROT(CROT(N),SROT(N))
             IF(LMAP)  CALL EQUID_CYLIND_MAP_JACOB(XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL EQUID_CYLIND_GRID_AREA(RLAT(N),AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             XPTS(N)=1+HI*MOD(HI*(RLON(N)-RLON1)+3600,360.)/DLON
             YPTS(N)=1+(RLAT(N)-RLAT1)/DLAT
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND.  &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT)  CALL EQUID_CYLIND_VECT_ROT(CROT(N),SROT(N))
                IF(LMAP)  CALL EQUID_CYLIND_MAP_JACOB(XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL EQUID_CYLIND_GRID_AREA(RLAT(N),AREA(N))
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
  END SUBROUTINE GDSWZD_EQUID_CYLIND

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
  SUBROUTINE EQUID_CYLIND_VECT_ROT(CROT,SROT)
    IMPLICIT NONE

    REAL,                INTENT(  OUT) :: CROT, SROT

    CROT=1.0
    SROT=0.0

  END SUBROUTINE EQUID_CYLIND_VECT_ROT

  !> Computes the map jacobians for a equidistant cylindrical grid.
  !!
  !! @param[out] xlon dx/dlon in 1/degrees.
  !! @param[out] xlat dx/dlat in 1/degrees.
  !! @param[out] ylon dy/dlon in 1/degrees.
  !! @param[out] ylat dy/dlat in 1/degrees.
  !!
  !! @author George Gayno
  !! @date July 2021
  SUBROUTINE EQUID_CYLIND_MAP_JACOB(XLON,XLAT,YLON,YLAT)
    REAL,                INTENT(  OUT) :: XLON,XLAT,YLON,YLAT

    XLON=1.0/DLON
    XLAT=0.
    YLON=0.
    YLAT=1.0/DLAT

  END SUBROUTINE EQUID_CYLIND_MAP_JACOB

  !> Computes the grid box area for a equidistant cylindrical grid.
  !!
  !! @param[in] rlat Latitude of grid point in degrees.
  !! @param[out] area Area weights in m^2.
  !!
  !! @author Mark Iredell, George Gayno
  !! @date July 2021
  SUBROUTINE EQUID_CYLIND_GRID_AREA(RLAT,AREA)
    IMPLICIT NONE

    REAL,                INTENT(IN   ) :: RLAT
    REAL,                INTENT(  OUT) :: AREA

    REAL,                PARAMETER     :: PI=3.14159265358979
    REAL,                PARAMETER     :: DPR=180./PI

    REAL                               :: DSLAT, RLATU, RLATD

    RLATU=MIN(MAX(RLAT+DLAT/2,-90.),90.)
    RLATD=MIN(MAX(RLAT-DLAT/2,-90.),90.)
    DSLAT=SIN(RLATU/DPR)-SIN(RLATD/DPR)
    AREA=RERTH**2*ABS(DSLAT*DLON)/DPR

  END SUBROUTINE EQUID_CYLIND_GRID_AREA

end module ip_equid_cylind_grid_mod

