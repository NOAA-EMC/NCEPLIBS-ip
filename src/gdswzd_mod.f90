!> @file
!! @brief Driver module for gdswzd routines.
!!
!! @date Jan 2015
!! @author George Gayno, Kyle Gerheiser

!> Driver module for gdswzd routines.
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
MODULE GDSWZD_MOD
  use ip_grid_descriptor_mod
  use ip_grids_mod
  use ip_grid_mod
  use ip_grid_factory_mod
  
  IMPLICIT NONE

  PRIVATE

  public :: GDSWZD_2D_ARRAY_grib1, gdswzd_grib1, gdswzd

  INTERFACE GDSWZD
     MODULE PROCEDURE GDSWZD_1D_ARRAY
     MODULE PROCEDURE GDSWZD_2D_ARRAY
     MODULE PROCEDURE GDSWZD_SCALAR
     module procedure gdswzd_grib1
     module procedure GDSWZD_2D_ARRAY_grib1
     module procedure gdswzd_grid
  END INTERFACE GDSWZD


CONTAINS


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
  !! @param[out] xlon dx/dlon in 1/degrees
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
  subroutine gdswzd_grid(grid,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

    class(ip_grid), intent(in) :: grid
    INTEGER,        INTENT(IN   ) :: IOPT, NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,           INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)

    INTEGER                       :: IS1, IM, JM, NM, KSCAN, NSCAN, N
    INTEGER                       :: IOPF, NN, I, J
    INTEGER                       :: I_OFFSET_ODD, I_OFFSET_EVEN

    !  COMPUTE GRID COORDINATES FOR ALL GRID POINTS
    IF(IOPT.EQ.0) THEN
       IOPF=1

       im = grid%im
       jm = grid%jm
       nm = im * jm
       nscan = grid%nscan
       kscan = grid%kscan

       if (nm > npts) then
          RLAT=FILL
          RLON=FILL
          XPTS=FILL
          YPTS=FILL
          return
       end if

       select type(grid)
       type is(ip_rot_equid_cylind_egrid)
          if(kscan == 0) then
             is1 = (jm + 1) / 2
          else
             is1 = jm / 2
          end if

          DO N=1,NM
             IF(NSCAN.EQ.0) THEN
                J=(N-1)/IM+1
                I=(N-IM*(J-1))*2-MOD(J+KSCAN,2)
             ELSE
                NN=(N*2)-1+KSCAN
                I = (NN-1)/JM + 1
                J = MOD(NN-1,JM) + 1
                IF (MOD(JM,2)==0.AND.MOD(I,2)==0.AND.KSCAN==0) J = J + 1
                IF (MOD(JM,2)==0.AND.MOD(I,2)==0.AND.KSCAN==1) J = J - 1
             ENDIF
             XPTS(N)=IS1+(I-(J-KSCAN))/2
             YPTS(N)=(I+(J-KSCAN))/2
          ENDDO
          class default
          DO N=1,NM
             IF(NSCAN.EQ.0) THEN
                J=(N-1)/IM+1
                I=N-IM*(J-1)
             ELSE
                I=(N-1)/JM+1
                J=N-JM*(I-1)
             ENDIF
             XPTS(N)=I
             YPTS(N)=J
          ENDDO
       end select

       DO N=NM+1,NPTS
          XPTS(N)=FILL
          YPTS(N)=FILL
       ENDDO

    ELSE  ! IOPT /= 0
       IOPF=IOPT
    ENDIF ! IOPT CHECK

    call grid%gdswzd(IOPF,NPTS,FILL,  &
            XPTS,YPTS,RLON,RLAT,NRET, &
            CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

  end subroutine gdswzd_grid


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
  !! @param[out] xlon dx/dlon in 1/degrees
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
  SUBROUTINE GDSWZD_SCALAR(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

    IMPLICIT NONE
    !
    INTEGER,        INTENT(IN   ) :: IGDTNUM, IGDTLEN
    INTEGER,        INTENT(IN   ) :: IGDTMPL(IGDTLEN)
    INTEGER,        INTENT(IN   ) :: IOPT, NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON, RLAT
    REAL,           INTENT(INOUT) :: XPTS, YPTS
    REAL, OPTIONAL, INTENT(  OUT) :: CROT, SROT
    REAL, OPTIONAL, INTENT(  OUT) :: XLON, XLAT
    REAL, OPTIONAL, INTENT(  OUT) :: YLON, YLAT, AREA

    REAL                          :: RLONA(1),RLATA(1)
    REAL                          :: XPTSA(1),YPTSA(1)
    REAL                          :: CROTA(1),SROTA(1)
    REAL                          :: XLONA(1),XLATA(1)
    REAL                          :: YLONA(1),YLATA(1),AREAA(1)

    RLONA(1) = RLON
    RLATA(1) = RLAT
    XPTSA(1) = XPTS
    YPTSA(1) = YPTS

    NRET = 0

    ! CALL WITHOUT EXTRA FIELDS.

    IF (.NOT. PRESENT(CROT) .AND. &
         .NOT. PRESENT(SROT) .AND. &
         .NOT. PRESENT(XLON) .AND. &
         .NOT. PRESENT(XLAT) .AND. &
         .NOT. PRESENT(YLON) .AND. &
         .NOT. PRESENT(YLAT) .AND. &
         .NOT. PRESENT(AREA) ) THEN

       CALL GDSWZD_1D_ARRAY(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
            XPTSA,YPTSA,RLONA,RLATA,NRET)

       RLON = RLONA(1)
       RLAT = RLATA(1)
       XPTS = XPTSA(1)
       YPTS = YPTSA(1)

    ENDIF

    ! MIMIC CALL TO OLD 'GDSWIZ' ROUTINES.

    IF (PRESENT(CROT) .AND. &
         PRESENT(SROT) .AND. &
         .NOT. PRESENT(XLON) .AND. &
         .NOT. PRESENT(XLAT) .AND. &
         .NOT. PRESENT(YLON) .AND. &
         .NOT. PRESENT(YLAT) .AND. &
         .NOT. PRESENT(AREA) ) THEN

       CALL GDSWZD_1D_ARRAY(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
            XPTSA,YPTSA,RLONA,RLATA,NRET,CROTA,SROTA)

       RLON = RLONA(1)
       RLAT = RLATA(1)
       XPTS = XPTSA(1)
       YPTS = YPTSA(1)
       CROT = CROTA(1)
       SROT = SROTA(1)

    ENDIF

    ! MIMIC CALL TO OLD 'GDSWZD' ROUTINES.

    IF (PRESENT(CROT) .AND. &
         PRESENT(SROT) .AND. &
         PRESENT(XLON) .AND. &
         PRESENT(XLAT) .AND. &
         PRESENT(YLON) .AND. &
         PRESENT(YLAT) .AND. &
         PRESENT(AREA) ) THEN

       CALL GDSWZD_1D_ARRAY(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
            XPTSA,YPTSA,RLONA,RLATA,NRET, &
            CROTA,SROTA,XLONA,XLATA,YLONA,YLATA,AREAA)

       RLON = RLONA(1)
       RLAT = RLATA(1)
       XPTS = XPTSA(1)
       YPTS = YPTSA(1)
       CROT = CROTA(1)
       SROT = SROTA(1)
       XLON = XLONA(1)
       XLAT = XLATA(1)
       YLON = YLONA(1)
       YLAT = YLATA(1)
       AREA = AREAA(1)

    ENDIF

    RETURN

  END SUBROUTINE GDSWZD_SCALAR

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
  !! @param[out] xlon dx/dlon in 1/degrees
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
  SUBROUTINE GDSWZD_2D_ARRAY(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

    IMPLICIT NONE
    !
    INTEGER,        INTENT(IN   ) :: IGDTNUM, IGDTLEN
    INTEGER,        INTENT(IN   ) :: IGDTMPL(IGDTLEN)
    INTEGER,        INTENT(IN   ) :: IOPT, NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON(:,:),RLAT(:,:)
    REAL,           INTENT(INOUT) :: XPTS(:,:),YPTS(:,:)
    REAL, OPTIONAL, INTENT(  OUT) :: CROT(:,:),SROT(:,:)
    REAL, OPTIONAL, INTENT(  OUT) :: XLON(:,:),XLAT(:,:)
    REAL, OPTIONAL, INTENT(  OUT) :: YLON(:,:),YLAT(:,:),AREA(:,:)

    CALL GDSWZD_1D_ARRAY(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
         XPTS,YPTS,RLON,RLAT,NRET, &
         CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

  END SUBROUTINE GDSWZD_2D_ARRAY

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
  !! @param[out] xlon dx/dlon in 1/degrees
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
  SUBROUTINE GDSWZD_1D_ARRAY(IGDTNUM,IGDTMPL,IGDTLEN,IOPT,NPTS,FILL, &
       XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    INTEGER,        INTENT(IN   ) :: IGDTNUM, IGDTLEN
    INTEGER,        INTENT(IN   ) :: IGDTMPL(IGDTLEN)
    INTEGER,        INTENT(IN   ) :: IOPT, NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,           INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)

    type(grib2_descriptor) :: desc
    class(ip_grid), allocatable :: grid

    desc = init_descriptor(igdtnum, igdtlen, igdtmpl)
    grid = init_grid(desc)
    
    call gdswzd_grid(grid,IOPT,NPTS,FILL, &
         XPTS,YPTS,RLON,RLAT,NRET, &
         CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
  
  END SUBROUTINE GDSWZD_1D_ARRAY

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
  !! @param[out] xlon dx/dlon in 1/degrees
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
  SUBROUTINE GDSWZD_grib1(KGDS,IOPT,NPTS,FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
    INTEGER,        INTENT(IN   ) :: IOPT, KGDS(200), NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON(NPTS),RLAT(NPTS)
    REAL,           INTENT(INOUT) :: XPTS(NPTS),YPTS(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: CROT(NPTS),SROT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: XLON(NPTS),XLAT(NPTS)
    REAL, OPTIONAL, INTENT(  OUT) :: YLON(NPTS),YLAT(NPTS),AREA(NPTS)

    
    type(grib1_descriptor) :: desc
    class(ip_grid), allocatable :: grid

    desc = init_descriptor(kgds)
    grid = init_grid(desc)
    
    call gdswzd_grid(grid,IOPT,NPTS,FILL, &
         XPTS,YPTS,RLON,RLAT,NRET, &
         CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

  END SUBROUTINE GDSWZD_grib1

  
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
  !! @param[out] xlon dx/dlon in 1/degrees
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
  SUBROUTINE GDSWZD_2d_array_grib1(KGDS,IOPT,NPTS,FILL,XPTS,YPTS,RLON,RLAT,NRET, &
       CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)
   
    !$$$
    INTEGER,        INTENT(IN   ) :: IOPT, KGDS(200), NPTS
    INTEGER,        INTENT(  OUT) :: NRET
    !
    REAL,           INTENT(IN   ) :: FILL
    REAL,           INTENT(INOUT) :: RLON(:,:),RLAT(:,:)
    REAL,           INTENT(INOUT) :: XPTS(:,:),YPTS(:,:)
    REAL, OPTIONAL, INTENT(  OUT) :: CROT(:,:),SROT(:,:)
    REAL, OPTIONAL, INTENT(  OUT) :: XLON(:,:),XLAT(:,:)
    REAL, OPTIONAL, INTENT(  OUT) :: YLON(:,:),YLAT(:,:),AREA(:,:)

    
    type(grib1_descriptor) :: desc
    class(ip_grid), allocatable :: grid

    desc = init_descriptor(kgds)
    grid = init_grid(desc)
    
    call gdswzd_grid(grid,IOPT,NPTS,FILL, &
         XPTS,YPTS,RLON,RLAT,NRET, &
         CROT,SROT,XLON,XLAT,YLON,YLAT,AREA)

  END SUBROUTINE GDSWZD_2d_array_grib1



END MODULE GDSWZD_MOD
