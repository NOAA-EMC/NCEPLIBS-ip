!> @file
!! Gaussian grid coordinate transformations
!! @author Mark Iredell, George Gayno, Kyle Gerheiser
!! @date July 2021

!> Gaussian grid coordinate transformations
!!
!! @author George Gayno, Mark Iredell, Kyle Gerheiser
!! @date July 2021
module ip_gaussian_grid_mod
  use ip_grid_descriptor_mod
  use ip_grid_mod
  use earth_radius_mod
  use constants_mod
  implicit none

  private
  public :: ip_gaussian_grid

  type, extends(ip_grid) :: ip_gaussian_grid
     integer :: jh
     real :: dlon, rlat1, rlon1, rlon2, hi
     integer :: jg, jscan
   contains
     procedure :: init_grib1
     procedure :: init_grib2
     procedure :: gdswzd => gdswzd_gaussian
  end type ip_gaussian_grid

  INTEGER                        :: J1, JH
  REAL,            ALLOCATABLE   :: BLAT(:)
  REAL                           :: DLON, RERTH
  REAL,            ALLOCATABLE   :: YLAT_ROW(:)

contains

  !> Initializes a gaussian grid given a grib1 grib1_descriptor object.
  !! 
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Kyle Gerheiser
  subroutine init_grib1(self, g1_desc)
    class(ip_gaussian_grid), intent(inout) :: self
    type(grib1_descriptor), intent(in) :: g1_desc

    integer :: iscan, jg

    associate(kgds => g1_desc%gds)
      self%rerth = 6.3712E6
      self%eccen_squared = 0.0

      self%IM=KGDS(2)
      self%JM=KGDS(3)
      self%RLAT1=KGDS(4)*1.E-3
      self%RLON1=KGDS(5)*1.E-3
      self%RLON2=KGDS(8)*1.E-3
      self%JG=KGDS(10)*2
      ISCAN=MOD(KGDS(11)/128,2)
      self%JSCAN=MOD(KGDS(11)/64,2)
      self%HI=(-1.)**ISCAN
      self%JH=(-1)**self%JSCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)

      self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      self%nscan = mod(kgds(11) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0

      self%iwrap=nint(360 / abs(self%dlon))
      if(self%im < self%iwrap) self%iwrap = 0

      if(self%iwrap > 0 .and. mod(self%iwrap, 2) == 0) then
         jg=kgds(10)*2
         if(self%jm == self%jg) then
            self%jwrap1 = 1
            self%jwrap2 = 2 * self%jm + 1
         endif
      endif

    end associate
  end subroutine init_grib1

  !> Initializes a gaussian grid given a grib1 grib2_descriptor object.
  !! @param[inout] self The grid to initialize
  !! @param[in] g1_desc A grib1_descriptor
  !!
  !! @author Kyle Gerheiser
  !! @date July 2021
  subroutine init_grib2(self, g2_desc)
    class(ip_gaussian_grid), intent(inout) :: self
    type(grib2_descriptor), intent(in) :: g2_desc

    integer :: iscale, iscan, jg

    associate(igdtmpl => g2_desc%gdt_tmpl, igdtlen => g2_desc%gdt_len)
      call EARTH_RADIUS(igdtmpl, igdtlen, self%rerth, self%eccen_squared)

      self%IM=IGDTMPL(8)
      self%JM=IGDTMPL(9)
      ISCALE=IGDTMPL(10)*IGDTMPL(11)
      IF(ISCALE==0) ISCALE=10**6
      self%RLAT1=FLOAT(IGDTMPL(12))/FLOAT(ISCALE)
      self%RLON1=FLOAT(IGDTMPL(13))/FLOAT(ISCALE)
      self%RLON2=FLOAT(IGDTMPL(16))/FLOAT(ISCALE)
      self%JG=IGDTMPL(18)*2
      ISCAN=MOD(IGDTMPL(19)/128,2)
      self%JSCAN=MOD(IGDTMPL(19)/64,2)
      self%HI=(-1.)**ISCAN
      self%JH=(-1)**self%JSCAN
      self%DLON=self%HI*(MOD(self%HI*(self%RLON2-self%RLON1)-1+3600,360.)+1)/(self%IM-1)


      self%iwrap = nint(360 / abs(self%dlon))
      if(self%im < self%iwrap) self%iwrap = 0
      self%jwrap1 = 0
      self%jwrap2 = 0
      if(self%iwrap > 0 .and. mod(self%iwrap, 2) == 0) then
         jg = igdtmpl(18) * 2
         if(self%jm == jg) then
            self%jwrap1=1
            self%jwrap2 = 2 * self%jm + 1
         endif
      endif
      self%nscan = mod(igdtmpl(19) / 32, 2)
      self%nscan_field_pos = self%nscan
      self%kscan = 0
    end associate

  end subroutine init_grib2

  !> Calculates Earth coordinates (iopt = 1) or grid coorindates (iopt = -1)
  !! for Gaussian grids.
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
  SUBROUTINE GDSWZD_GAUSSIAN(self,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    IMPLICIT NONE
    !
    class(ip_gaussian_grid), intent(in) :: self
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
    INTEGER                        :: JSCAN, IM, JM
    INTEGER                        :: J, JA, JG
    INTEGER                        :: ISCALE, N
    !
    LOGICAL                        :: LROT, LMAP, LAREA
    !
    REAL,            ALLOCATABLE   :: ALAT(:), ALAT_JSCAN(:)
    REAL,            ALLOCATABLE   :: ALAT_TEMP(:),BLAT_TEMP(:)
    REAL                           :: HI, RLATA, RLATB, RLAT1, RLON1, RLON2
    REAL                           :: XMAX, XMIN, YMAX, YMIN, YPTSA, YPTSB
    REAL                           :: WB
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(PRESENT(CROT)) CROT=FILL
    IF(PRESENT(SROT)) SROT=FILL
    IF(PRESENT(XLON)) XLON=FILL
    IF(PRESENT(XLAT)) XLAT=FILL
    IF(PRESENT(YLON)) YLON=FILL
    IF(PRESENT(YLAT)) YLAT=FILL
    IF(PRESENT(AREA)) AREA=FILL
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

    IM=self%im
    JM=self%jm

    RLAT1=self%rlat1
    RLON1=self%rlon1
    RLON2=self%rlon2

    JG=self%jg
    JSCAN=self%jscan
    HI=self%hi

    JH=self%jh
    DLON=self%dlon
    rerth = self%rerth

    ALLOCATE(ALAT_TEMP(JG))
    ALLOCATE(BLAT_TEMP(JG))
    CALL SPLAT(4,JG,ALAT_TEMP,BLAT_TEMP)
    ALLOCATE(ALAT(0:JG+1))
    ALLOCATE(BLAT(0:JG+1))
    !$OMP PARALLEL DO PRIVATE(JA) SCHEDULE(STATIC)
    DO JA=1,JG
       ALAT(JA)=DPR*ASIN(ALAT_TEMP(JA))
       BLAT(JA)=BLAT_TEMP(JA)
    ENDDO
    !$OMP END PARALLEL DO
    DEALLOCATE(ALAT_TEMP,BLAT_TEMP)
    ALAT(0)=180.-ALAT(1)
    ALAT(JG+1)=-ALAT(0)
    BLAT(0)=-BLAT(1)
    BLAT(JG+1)=BLAT(0)
    J1=1
    DO WHILE(J1.LT.JG.AND.RLAT1.LT.(ALAT(J1)+ALAT(J1+1))/2)
       J1=J1+1
    ENDDO
    IF(LMAP)THEN
       ALLOCATE(ALAT_JSCAN(JG))
       DO JA=1,JG
          ALAT_JSCAN(J1+JH*(JA-1))=ALAT(JA)
       ENDDO
       ALLOCATE(YLAT_ROW(0:JG+1))
       DO JA=2,(JG-1)
          YLAT_ROW(JA)=2.0/(ALAT_JSCAN(JA+1)-ALAT_JSCAN(JA-1))
       ENDDO
       YLAT_ROW(1)=1.0/(ALAT_JSCAN(2)-ALAT_JSCAN(1))
       YLAT_ROW(0)=YLAT_ROW(1)
       YLAT_ROW(JG)=1.0/(ALAT_JSCAN(JG)-ALAT_JSCAN(JG-1))
       YLAT_ROW(JG+1)=YLAT_ROW(JG)
       DEALLOCATE(ALAT_JSCAN)
    ENDIF
    XMIN=0
    XMAX=IM+1
    IF(IM.EQ.NINT(360/ABS(DLON))) XMAX=IM+2
    YMIN=0.5
    YMAX=JM+0.5
    NRET=0
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  TRANSLATE GRID COORDINATES TO EARTH COORDINATES
    IF(IOPT.EQ.0.OR.IOPT.EQ.1) THEN
       !$OMP PARALLEL DO PRIVATE(N,J,WB,RLATA,RLATB) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
               YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
             RLON(N)=MOD(RLON1+DLON*(XPTS(N)-1)+3600,360.)
             J=YPTS(N)
             WB=YPTS(N)-J
             RLATA=ALAT(J1+JH*(J-1))
             RLATB=ALAT(J1+JH*J)
             RLAT(N)=RLATA+WB*(RLATB-RLATA)
             NRET=NRET+1
             IF(LROT) CALL GAUSSIAN_VECT_ROT(CROT(N),SROT(N))
             IF(LMAP) CALL GAUSSIAN_MAP_JACOB(YPTS(N),&
                  XLON(N),XLAT(N),YLON(N),YLAT(N))
             IF(LAREA) CALL GAUSSIAN_GRID_AREA(YPTS(N),AREA(N))
          ELSE
             RLON(N)=FILL
             RLAT(N)=FILL
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
       !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       !  TRANSLATE EARTH COORDINATES TO GRID COORDINATES
    ELSEIF(IOPT.EQ.-1) THEN
       !$OMP PARALLEL DO PRIVATE(N,JA,YPTSA, YPTSB, WB) REDUCTION(+:NRET) SCHEDULE(STATIC)
       DO N=1,NPTS
          XPTS(N)=FILL
          YPTS(N)=FILL
          IF(ABS(RLON(N)).LE.360.AND.ABS(RLAT(N)).LE.90) THEN
             XPTS(N)=1+HI*MOD(HI*(RLON(N)-RLON1)+3600,360.)/DLON
             JA=MIN(INT((JG+1)/180.*(90-RLAT(N))),JG)
             IF(RLAT(N).GT.ALAT(JA)) JA=MAX(JA-2,0)
             IF(RLAT(N).LT.ALAT(JA+1)) JA=MIN(JA+2,JG)
             IF(RLAT(N).GT.ALAT(JA)) JA=JA-1
             IF(RLAT(N).LT.ALAT(JA+1)) JA=JA+1
             YPTSA=1+JH*(JA-J1)
             YPTSB=1+JH*(JA+1-J1)
             WB=(ALAT(JA)-RLAT(N))/(ALAT(JA)-ALAT(JA+1))
             YPTS(N)=YPTSA+WB*(YPTSB-YPTSA)
             IF(XPTS(N).GE.XMIN.AND.XPTS(N).LE.XMAX.AND. &
                  YPTS(N).GE.YMIN.AND.YPTS(N).LE.YMAX) THEN
                NRET=NRET+1
                IF(LROT) CALL GAUSSIAN_VECT_ROT(CROT(N),SROT(N))
                IF(LMAP) CALL GAUSSIAN_MAP_JACOB(YPTS(N), &
                     XLON(N),XLAT(N),YLON(N),YLAT(N))
                IF(LAREA) CALL GAUSSIAN_GRID_AREA(YPTS(N),AREA(N))
             ELSE
                XPTS(N)=FILL
                YPTS(N)=FILL
             ENDIF
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF
    DEALLOCATE(ALAT, BLAT)
    IF (ALLOCATED(YLAT_ROW)) DEALLOCATE(YLAT_ROW)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE GDSWZD_GAUSSIAN
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE GAUSSIAN_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !$$$  SUBPROGRAM DOCUMENTATION BLOCK
    !
    ! SUBPROGRAM:  GAUSSIAN_ERROR   ERROR HANDLER
    !   PRGMMR: GAYNO       ORG: W/NMC23       DATE: 2015-07-13
    !
    ! ABSTRACT: UPON AN ERROR, THIS SUBPROGRAM ASSIGNS
    !           A "FILL" VALUE TO THE OUTPUT FIELDS.
    !
    ! PROGRAM HISTORY LOG:
    ! 2015-07-13  GAYNO     INITIAL VERSION
    !
    ! USAGE:    CALL GAUSSIAN_ERROR(IOPT,FILL,RLAT,RLON,XPTS,YPTS,NPTS)
    !
    !   INPUT ARGUMENT LIST:
    !     IOPT     - INTEGER OPTION FLAG
    !                (+1 TO COMPUTE EARTH COORDS OF SELECTED GRID COORDS)
    !                (-1 TO COMPUTE GRID COORDS OF SELECTED EARTH COORDS)
    !     NPTS     - INTEGER MAXIMUM NUMBER OF COORDINATES
    !     FILL     - REAL FILL VALUE TO SET INVALID OUTPUT DATA
    !                (MUST BE IMPOSSIBLE VALUE; SUGGESTED VALUE: -9999.)
    !   OUTPUT ARGUMENT LIST:
    !     RLON     - REAL (NPTS) EARTH LONGITUDES IN DEGREES E IF IOPT<0
    !     RLAT     - REAL (NPTS) EARTH LATITUDES IN DEGREES N IF IOPT<0
    !     XPTS     - REAL (NPTS) GRID X POINT COORDINATES IF IOPT>0
    !     YPTS     - REAL (NPTS) GRID Y POINT COORDINATES IF IOPT>0
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !
    !$$$
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
  END SUBROUTINE GAUSSIAN_ERROR

  !> This subprogram computes the vector rotation sines and
  !! cosines for a gaussian cylindrical grid.
  !!
  !! @param[out] crot Clockwise vector rotation cosines.
  !! @param[out] srot Clockwise vector rotation sines.
  !!
  !! @note
  !! ugrid=crot*uearth-srot*vearth;
  !! vgrid=srot*uearth+crot*vearth)
  !!
  !! @author George Gayno
  !! @date July 2021
  SUBROUTINE GAUSSIAN_VECT_ROT(CROT,SROT)
    IMPLICIT NONE

    REAL,                INTENT(  OUT) :: CROT, SROT

    CROT=1.0
    SROT=0.0

  END SUBROUTINE GAUSSIAN_VECT_ROT

  
  !> Computes the map jacobians for a gaussian cylindrical grid.
  !!
  !! @param[in] ypts y-index of grid point.
  !! @param[out] xlon dx/dlon in 1/degrees.
  !! @param[out] xlat dx/dlat in 1/degrees.
  !! @param[out] ylon dy/dlon in 1/degrees.
  !! @param[out] ylat dy/dlat in 1/degrees.
  !!
  !! @author George Gayno
  !! @date July 2021
  SUBROUTINE GAUSSIAN_MAP_JACOB(YPTS, XLON, XLAT, YLON, YLAT)
    IMPLICIT NONE

    REAL,                INTENT(IN   ) :: YPTS
    REAL,                INTENT(  OUT) :: XLON, XLAT, YLON, YLAT

    XLON=1/DLON
    XLAT=0.
    YLON=0.
    YLAT=YLAT_ROW(NINT(YPTS))

  END SUBROUTINE GAUSSIAN_MAP_JACOB

  !> Computes the grid box area for a gaussian cylindrical grid.
  !!
  !! @param[in] ypts y-index of grid point.
  !! @param[out] area Area weights in m^2
  !!
  !! @author Mark Iredell, George Gayno
  !! @date July 2021
  SUBROUTINE GAUSSIAN_GRID_AREA(YPTS,AREA)
    IMPLICIT NONE

    REAL,            INTENT(IN   ) :: YPTS
    REAL,            INTENT(  OUT) :: AREA

    INTEGER                        :: J

    REAL                           :: WB, WLAT, WLATA, WLATB

    J = YPTS
    WB=YPTS-J
    WLATA=BLAT(J1+JH*(J-1))
    WLATB=BLAT(J1+JH*J)
    WLAT=WLATA+WB*(WLATB-WLATA)
    AREA=RERTH**2*WLAT*DLON/DPR

  END SUBROUTINE GAUSSIAN_GRID_AREA


end module ip_gaussian_grid_mod

