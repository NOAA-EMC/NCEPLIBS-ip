!> @file
!> @brief Expand or contract eta grids.
!> @author Iredell @date 96-04-10

!> Expand or contract eta grids.
!>
!> This subprogram transforms between the staggered eta grids as used
!> in the eta model and for native grid transmission and their full
!> expansion as used for general interpolation and graphics. The eta
!> grids are rotated latitude-longitude grids staggered as defined by
!> the arakawa e-grid, that is with mass data points alternating with
!> wind data points.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | Iredell | Initial
!> 2015-07-14 | Gayno | Make grib 2 compliant. Replace 4-pt interpolation with call to ipolates.
!>
!> @param[in] idir integer transform option.
!> - 0 to expand staggered fields to full fields
!> - -1 to contract full mass fields to staggered fields
!> - -2 to contract full wind fields to staggered fields
!> @param[in] igdtnumi integer grid definition template number - input
!> grid. Corresponds to the gfld%igdtnum component of the ncep g2
!> library gridmod data structure. Must be = 1 (for a rotated lat/lon
!> grid.)
!> @param[in] igdtlen integer number of elements of the grid
!> definition template array - same for input and output grids (=22)
!> which are both rotated lat/lon grids.  corresponds to the
!> gfld%igdtlen component of the ncep g2 library gridmod data
!> structure.

!> @param[in] igdtmpli integer (igdtlen) grid definition template
!> array - input grid. Corresponds to the gfld%igdtmpl component of
!> the ncep g2 library gridmod data structure (section 3 info):

!> - 1 shape of earth, octet 15
!> - 2 scale factor of spherical earth radius, octet 16
!> - 3 scaled value of radius of spherical earth, octets 17-20
!> - 4 scale factor of major axis of elliptical earth, octet 21
!> - 5 scaled value of major axis of elliptical earth, octets 22-25
!> - 6 scale factor of minor axis of elliptical earth, octet 26
!> - 7 scaled value of minor axis of elliptical earth, octets 27-30
!> - 8 number of points along a parallel, octs 31-34
!> - 9 number of points along a meridian, octs 35-38
!> - 10 basic angle of initial production domain, octets 39-42
!> - 11 subdivisions of basic angle, octets 43-46
!> - 12 latitude of first grid point, octets 47-50
!> - 13 longitude of first grid point, octets 51-54
!> - 14 resolution and component flags, octet 55
!> - 15 latitude of last grid point, octets 56-59
!> - 16 longitude of last grid point, octets 60-63
!> - 17 i-direction increment, octets 64-67
!> - 18 j-direction increment, octets 68-71
!> - 19 scanning mode, octet 72
!> - 20 latitude of southern pole of projection, octets 73-76
!> - 21 longitude of southern pole of projection, octets 77-80
!> - 22 angle of rotation of projection, octs 81-84
!> @param[in] npts_input integer number points input grid
!> @param[in] bitmap_input logical (npts_input) input grid bitmap
!> @param[in] data_input real (npts_input) input grid data
!> @param[out] igdtnumo integer grid definition template number -
!> output grid. Corresponds to the gfld%igdtnum component of the ncep
!> g2 library gridmod data structure. Same as igdtnumi (=1 for a
!> rotated lat/lon grid).
!> @param[out] igdtmplo integer (igdtlen) grid definition template
!> array - output grid. Corresponds to the gfld%igdtmpl component of
!> the ncep g2 library gridmod data structure. Array definitions same
!> as "igdtmpli".
!> @param[out] npts_output integer number points output grid. the
!> j-dimension of the input and output grids are the same. When going
!> from a staggered to a full grid the i-dimension increases to
!> idim*2-1. When going from full to staggered the i-dimension
!> decreases to (idim+1)/2.
!>  @param[out] bitmap_output logical (npts_outut) output grid bitmap
!>  @param[out] data_output real (npts_output) output grid data
!> @param[out] iret integer return code
!> - 0 successful transformation
!> - non-0 invalid grid specs or problem in ipolates().
!>
! @author Iredell @date 96-04-10
 SUBROUTINE IPXETAS(IDIR, IGDTNUMI, IGDTLEN, IGDTMPLI, NPTS_INPUT,  &
                    BITMAP_INPUT, DATA_INPUT, IGDTNUMO, IGDTMPLO, &
                    NPTS_OUTPUT, BITMAP_OUTPUT, DATA_OUTPUT, IRET)
 IMPLICIT NONE
!
 INTEGER,         INTENT(IN   )    :: IDIR
 INTEGER,         INTENT(IN   )    :: IGDTNUMI, IGDTLEN
 INTEGER,         INTENT(IN   )    :: IGDTMPLI(IGDTLEN)
 INTEGER,         INTENT(IN   )    :: NPTS_INPUT, NPTS_OUTPUT
 INTEGER,         INTENT(  OUT)    :: IGDTNUMO
 INTEGER,         INTENT(  OUT)    :: IGDTMPLO(IGDTLEN)
 INTEGER,         INTENT(  OUT)    :: IRET

 LOGICAL(KIND=1), INTENT(IN   )    :: BITMAP_INPUT(NPTS_INPUT)
 LOGICAL(KIND=1), INTENT(  OUT)    :: BITMAP_OUTPUT(NPTS_OUTPUT)

 REAL,            INTENT(IN   )    :: DATA_INPUT(NPTS_INPUT)
 REAL,            INTENT(  OUT)    :: DATA_OUTPUT(NPTS_OUTPUT)

 INTEGER                           :: SCAN_MODE, ISCALE, IP, IPOPT(20)
 INTEGER                           :: IBI(1), IBO(1), J, KM, NO

 REAL                              :: DLONS
 REAL, ALLOCATABLE                 :: OUTPUT_RLAT(:), OUTPUT_RLON(:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 IRET = 0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ROUTINE ONLY WORKS FOR ROTATED LAT/LON GRIDS.
 IF (IGDTNUMI/=1) THEN
   IRET=1
   RETURN
 ENDIF
!
 SCAN_MODE=IGDTMPLI(19)
 IF((SCAN_MODE==68.OR.SCAN_MODE==72).AND.(IDIR<-2.OR.IDIR>-1))THEN
   IGDTNUMO=IGDTNUMI
   IGDTMPLO=IGDTMPLI
   IGDTMPLO(19)=64
   IGDTMPLO(8)=IGDTMPLO(8)*2-1
   IF((IGDTMPLO(8)*IGDTMPLO(9))/=NPTS_OUTPUT)THEN
     IRET=3
     RETURN
   ENDIF
   ISCALE=IGDTMPLO(10)*IGDTMPLO(11)
   IF(ISCALE==0) ISCALE=10**6
   DLONS=FLOAT(IGDTMPLO(17))/FLOAT(ISCALE)
   DLONS=DLONS*0.5
   IGDTMPLO(17)=NINT(DLONS*FLOAT(ISCALE))
 ELSEIF(SCAN_MODE==64.AND.IDIR==-1)THEN  ! FULL TO H-GRID
   IGDTNUMO=IGDTNUMI
   IGDTMPLO=IGDTMPLI
   IGDTMPLO(19)=68
   IGDTMPLO(8)=(IGDTMPLO(8)+1)/2
   IF((IGDTMPLO(8)*IGDTMPLO(9))/=NPTS_OUTPUT)THEN
     IRET=3
     RETURN
   ENDIF
   ISCALE=IGDTMPLO(10)*IGDTMPLO(11)
   IF(ISCALE==0) ISCALE=10**6
   DLONS=FLOAT(IGDTMPLO(17))/FLOAT(ISCALE)
   DLONS=DLONS*2.0
   IGDTMPLO(17)=NINT(DLONS*FLOAT(ISCALE))
 ELSEIF(SCAN_MODE==64.AND.IDIR==-2)THEN  ! FULL TO V-GRID
   IGDTNUMO=IGDTNUMI
   IGDTMPLO=IGDTMPLI
   IGDTMPLO(19)=72
   IGDTMPLO(8)=(IGDTMPLO(8)+1)/2
   IF((IGDTMPLO(8)*IGDTMPLO(9))/=NPTS_OUTPUT)THEN
     IRET=3
     RETURN
   ENDIF
   ISCALE=IGDTMPLO(10)*IGDTMPLO(11)
   IF(ISCALE==0) ISCALE=10**6
   DLONS=FLOAT(IGDTMPLO(17))/FLOAT(ISCALE)
   DLONS=DLONS*2.0
   IGDTMPLO(17)=NINT(DLONS*FLOAT(ISCALE))
 ELSE
   IRET=2
   RETURN
 ENDIF

 KM=1
 IP=0
 IPOPT=0
 IBI=1
 IBO=0

 ALLOCATE(OUTPUT_RLAT(NPTS_OUTPUT))
 ALLOCATE(OUTPUT_RLON(NPTS_OUTPUT))

 CALL IPOLATES(IP, IPOPT, IGDTNUMI, IGDTMPLI, IGDTLEN, &
               IGDTNUMO, IGDTMPLO, IGDTLEN, &
               NPTS_INPUT, NPTS_OUTPUT, KM, IBI, BITMAP_INPUT, DATA_INPUT, &
               NO, OUTPUT_RLAT, OUTPUT_RLON, IBO, BITMAP_OUTPUT, DATA_OUTPUT, IRET)

 DEALLOCATE(OUTPUT_RLAT, OUTPUT_RLON)

 IF(IRET /= 0)THEN
   PRINT*,'- PROBLEM IN IPOLATES: ', IRET
   RETURN
 ENDIF

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! REPLACE ANY UNDEFINED POINTS ALONG THE LEFT AND RIGHT EDGES.
 DO J=1, IGDTMPLO(9)
   BITMAP_OUTPUT(J*IGDTMPLO(8))=BITMAP_OUTPUT(J*IGDTMPLO(8)-1)
   DATA_OUTPUT(J*IGDTMPLO(8))=DATA_OUTPUT(J*IGDTMPLO(8)-1)
   BITMAP_OUTPUT((J-1)*IGDTMPLO(8)+1)=BITMAP_OUTPUT((J-1)*IGDTMPLO(8)+2)
   DATA_OUTPUT((J-1)*IGDTMPLO(8)+1)=DATA_OUTPUT((J-1)*IGDTMPLO(8)+2)
 ENDDO

 RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 END SUBROUTINE IPXETAS
