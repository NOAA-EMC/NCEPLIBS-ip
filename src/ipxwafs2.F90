!> @file
!> @brief Expand or contract wafs grids using linear interpolation and
!> account for bitmapped data.
!> @author Iredell @date 96-04-10

!> Expand or contract wafs grids using linear interpolation and
!> account for bitmapped data.
!>
!> This subprogram transforms between the thinned wafs grids used for
!> transmitting to the aviation community and their full expansion as
!> used for general interpolation and graphics. The thinned wafs grids
!> are latitude-longitude grids where the number of points in each row
!> decrease toward the pole. This information is stored in the grib 2
!> grid definition template (section 3) starting at octet 73.
!>
!> The full grid counterparts have an equal number of points per row.
!>
!> The transform between the full and thinned wafs wafs grid is done
!> by linear interpolation and is not reversible. This routine works
!> with bitmapped data.
!>
!> This subroutine is similar to:
!> - ipxwafs() which uses linear interpolation.
!> - ipxwafs3() which uses neighbor interpolation and accounts for
!> - bitmapped data.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | iredell | initial version
!> 99-01-25 | gilbert | changed bitmap fields from logical to logical*1
!> 2015-jul | gayno | convert to grib 2
!>
!> @param[in] idir integer transform option
!> - 1 to expand thinned fields to full fields
!> - -1 to contract full fields to thinned fields
!> @param[in] numpts_thin integer number of grid points - thinned
!> grid. Must be 3447.
!> @param[in] numpts_full integer number of grid points - full
!> grid. Must be 5329.
!> @param[in] km integer number of fields to transform
!> @param[in] num_opt integer number of values to describe the thinned
!> grid.  must be 73.  dimension of array opt_pts.
!> @param[inout] opt_pts integer (num_opt) number of grid points per row -
!> thinned grid - if idir=+1
!> @param[in] igdtlen integer grid defintion template array length.
!> must be 19 for lat/lon grids. corresponds to the gfld%igdtlen
!> component of the ncep g2 library gridmod data structure.  same for
!> thin and full grids which are both lat/lon.
!> @param[in] igdtmpl_thin integer (igdtlen) grid definition template
!> array - thinned grid - if idir=+1. corresponds to the gfld%igdtmpl
!> component of the ncep g2 library gridmod data structure (section 3
!> info):
!> - 1 shape of earth, octet 15
!> - 2 scale factor of spherical earth radius, octet 16
!> - 3 scaled value of radius of spherical earth, octets 17-20
!> - 4 scale factor of major axis of elliptical earth, octet 21
!> - 5 scaled value of major axis of elliptical earth, octets 22-25
!> - 6 scale factor of minor axis of elliptical earth, octet 26
!> - 7 scaled value of minor axis of elliptical earth, octets 27-30
!> - 8 set to missing for thinned grid., octs 31-34
!> - 9 number of points along a meridian, octs 35-38
!> - 10 basic angle of initial production domain, octets 39-42.
!> - 11 subdivisions of basic angle, octets 43-46
!> - 12 latitude of first grid point, octets 47-50
!> - 13 longitude of first grid point, octets 51-54
!> - 14 resolution and component flags, octet 55
!> - 15 latitude of last grid point, octets 56-59
!> - 16 longitude of last grid point, octets 60-63
!> - 17 set to missing for thinned grid, octets 64-67
!> - 18 j-direction increment, octets 68-71
!> - 19 scanning mode, octet 72
!> @param[inout] data_thin real (numpts_thin,km) thinned grid fields if idir=+1
!> @param[inout] ib_thin integer (km) bitmap flags thinned grid - if idir=+1
!> @param[inout] bitmap_thin logical (numpts_thin,km) bitmap fields thin grid - if idir=+1
!> @param[inout] igdtmpl_full integer (igdtlen) grid definition template
!> array - full grid - if idir=-1. corresponds to the gfld%igdtmpl
!> component of the ncep g2 library gridmod data structure. same as
!> igdtmpl_thin except: (8): number of points along a parallel, octs
!> 31-34 (17): i-direction increment, octets 64-67
!> @param[inout] data_full real (numpts_full,km) full grid fields if idir=-1
!> @param[inout] ib_full integer (km) bitmap flags full grid - if idir=-1
!> @param[inout] bitmap_full logical (numpts_full,km) bitmap fields full grid - if idir=-1
!> @param[out] iret integer return code
!>  - 0 successful transformation
!>  - 1 improper grid specification
!>
!> @author Iredell @date 96-04-10
SUBROUTINE IPXWAFS2(IDIR, NUMPTS_THIN, NUMPTS_FULL, KM, NUM_OPT, OPT_PTS, &
     IGDTLEN, IGDTMPL_THIN, DATA_THIN, IB_THIN, BITMAP_THIN,  &
     IGDTMPL_FULL, DATA_FULL, IB_FULL, BITMAP_FULL, IRET)
  IMPLICIT NONE
  !
  INTEGER,               INTENT(IN   ) :: NUM_OPT
  INTEGER,               INTENT(INOUT) :: OPT_PTS(NUM_OPT)
  INTEGER,               INTENT(IN   ) :: IDIR, KM, NUMPTS_THIN, NUMPTS_FULL
  INTEGER,               INTENT(IN   ) :: IGDTLEN
  INTEGER,               INTENT(INOUT) :: IGDTMPL_THIN(IGDTLEN)
  INTEGER,               INTENT(INOUT) :: IGDTMPL_FULL(IGDTLEN)
  INTEGER,               INTENT(INOUT) :: IB_THIN(KM), IB_FULL(KM)
  INTEGER,               INTENT(  OUT) :: IRET
  !
  LOGICAL(KIND=1),       INTENT(INOUT) :: BITMAP_THIN(NUMPTS_THIN,KM)
  LOGICAL(KIND=1),       INTENT(INOUT) :: BITMAP_FULL(NUMPTS_FULL,KM)
  !
  REAL,                  INTENT(INOUT) :: DATA_THIN(NUMPTS_THIN,KM)
  REAL,                  INTENT(INOUT) :: DATA_FULL(NUMPTS_FULL,KM)
  !
  INTEGER,               PARAMETER     :: MISSING=-1
  !
  INTEGER                              :: SCAN_MODE, I, J, K, IDLAT, IDLON
  INTEGER                              :: IA, IB, IM, IM1, IM2, NPWAFS(73)
  INTEGER                              :: IS1, IS2, ISCAN, ISCALE
  !
  LOGICAL                              :: TEST1, TEST2
  !
  REAL                                 :: DLON, HI
  REAL                                 :: RAT1, RAT2, RLON1, RLON2
  REAL                                 :: WA, WB, X1, X2
  !
  DATA NPWAFS/ &
       73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70,&
       70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60,&
       59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43,&
       42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22,&
       20, 19, 17, 16, 14, 12, 11,  9,  8,  6,  5,  3,  2/
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  TRANSFORM GDS
  IRET=0
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  REG LAT/LON GRIDS HAVE 19 GDT ELEMENTS.
  IF (IGDTLEN /= 19 .OR. NUMPTS_THIN/=3447 .OR. NUMPTS_FULL/=5329) THEN
     IRET=1
     RETURN
  ENDIF
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  EXPAND THINNED GDS TO FULL GDS
  IF(IDIR.GT.0) THEN
     SCAN_MODE=IGDTMPL_THIN(19)
     ISCALE=IGDTMPL_THIN(10)*IGDTMPL_THIN(11)
     IF(ISCALE==0) ISCALE=10**6
     IDLAT=NINT(1.25*FLOAT(ISCALE))
     TEST1=ALL(OPT_PTS==NPWAFS)
     TEST2=ALL(OPT_PTS==NPWAFS(73:1:-1))
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     !  SOME CHECKS TO ENSURE THIS IS A WAFS GRID
     IF(SCAN_MODE==64 .AND. IGDTMPL_THIN(9)==73 .AND. &
          IDLAT==IGDTMPL_THIN(18) .AND. (TEST1 .OR. TEST2) ) THEN
        IGDTMPL_FULL=IGDTMPL_THIN
        IM=73
        IGDTMPL_FULL(8)=IM
        RLON1=FLOAT(IGDTMPL_FULL(13))/FLOAT(ISCALE)
        RLON2=FLOAT(IGDTMPL_FULL(16))/FLOAT(ISCALE)
        ISCAN=MOD(IGDTMPL_FULL(19)/128,2)
        HI=(-1.)**ISCAN
        DLON=HI*(MOD(HI*(RLON2-RLON1)-1+3600,360.)+1)/(IM-1)
        IGDTMPL_FULL(17)=NINT(DLON*FLOAT(ISCALE))
     ELSE
        IRET=1
     ENDIF
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     !  CONTRACT FULL GDS TO THINNED GDS
  ELSEIF(IDIR.LT.0) THEN
     SCAN_MODE=IGDTMPL_FULL(19)
     ISCALE=IGDTMPL_FULL(10)*IGDTMPL_FULL(11)
     IF(ISCALE==0) ISCALE=10**6
     IDLAT=NINT(1.25*FLOAT(ISCALE))
     IDLON=IDLAT
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     !  SOME CHECKS TO ENSURE THIS IS A WAFS GRID
     IF(SCAN_MODE==64 .AND. IGDTMPL_FULL(8)==73 .AND. IGDTMPL_FULL(9)==73 .AND. &
          NUM_OPT==73 .AND. IDLAT==IGDTMPL_FULL(18) .AND. IDLON==IGDTMPL_FULL(17))THEN
        IGDTMPL_THIN=IGDTMPL_FULL
        IGDTMPL_THIN(8)=MISSING
        IGDTMPL_THIN(17)=MISSING
        IF(IGDTMPL_THIN(12)==0) THEN  ! IS LATITUDE OF ROW 1 THE EQUATOR?
           OPT_PTS=NPWAFS
        ELSE
           OPT_PTS=NPWAFS(73:1:-1)
        ENDIF
     ELSE
        IRET=1
     ENDIF
  ENDIF
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !  TRANSFORM FIELDS
  IF(IRET.EQ.0) THEN
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     !  EXPAND THINNED FIELDS TO FULL FIELDS
     IF(IDIR.EQ.1) THEN
        DO K=1,KM
           IS1=0
           IS2=0
           IB_FULL(K)=0
           DO J=1,IGDTMPL_FULL(9)
              IM1=OPT_PTS(J)
              IM2=IGDTMPL_FULL(8)
              RAT1=FLOAT(IM1-1)/FLOAT(IM2-1)
              DO I=1,IM2
                 X1=(I-1)*RAT1+1
                 IA=X1
                 IA=MIN(MAX(IA,1),IM1-1)
                 IB=IA+1
                 WA=IB-X1
                 WB=X1-IA
                 IF(IB_THIN(K)==0.OR.(BITMAP_THIN(IS1+IA,K).AND.BITMAP_THIN(IS1+IB,K))) THEN
                    DATA_FULL(IS2+I,K)=WA*DATA_THIN(IS1+IA,K)+WB*DATA_THIN(IS1+IB,K)
                    BITMAP_FULL(IS2+I,K)=.TRUE.
                 ELSE
                    DATA_FULL(IS2+I,K)=0.0
                    BITMAP_FULL(IS2+I,K)=.FALSE.
                    IB_FULL(K)=1
                 ENDIF
              ENDDO
              IS1=IS1+IM1
              IS2=IS2+IM2
           ENDDO
        ENDDO
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !  CONTRACT FULL FIELDS TO THINNED FIELDS
     ELSEIF(IDIR.EQ.-1) THEN
        DO K=1,KM
           IS1=0
           IS2=0
           IB_THIN(K)=0
           DO J=1,IGDTMPL_FULL(9)
              IM1=OPT_PTS(J)
              IM2=IGDTMPL_FULL(8)
              RAT2=FLOAT(IM2-1)/FLOAT(IM1-1)
              DO I=1,IM1
                 X2=(I-1)*RAT2+1
                 IA=X2
                 IA=MIN(MAX(IA,1),IM2-1)
                 IB=IA+1
                 WA=IB-X2
                 WB=X2-IA
                 IF(IB_FULL(K)==0.OR.(BITMAP_FULL(IS2+IA,K).AND.BITMAP_FULL(IS2+IB,K))) THEN
                    DATA_THIN(IS1+I,K)=WA*DATA_FULL(IS2+IA,K)+WB*DATA_FULL(IS2+IB,K)
                    BITMAP_THIN(IS1+I,K)=.TRUE.
                 ELSE
                    DATA_THIN(IS1+I,K)=0.0
                    BITMAP_THIN(IS1+I,K)=.FALSE.
                    IB_THIN(K)=1
                 ENDIF
              ENDDO
              IS1=IS1+IM1
              IS2=IS2+IM2
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
END SUBROUTINE IPXWAFS2
