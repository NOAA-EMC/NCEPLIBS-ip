!> @file
!! IREDELL'S POLATE FOR SCALAR FIELDS
!! @author IREDELL @date 96-04-10
!
!> THIS SUBPROGRAM INTERPOLATES SCALAR FIELDS
!!           FROM ANY GRID TO ANY GRID (JOE IRWIN'S DREAM).
!!           ONLY HORIZONTAL INTERPOLATION IS PERFORMED.
!!           THE FOLLOWING INTERPOLATION METHODS ARE POSSIBLE:
!!             (IP=0) BILINEAR
!!             (IP=1) BICUBIC
!!             (IP=2) NEIGHBOR
!!             (IP=3) BUDGET
!!             (IP=4) SPECTRAL
!!             (IP=6) NEIGHBOR-BUDGET
!!           SOME OF THESE METHODS HAVE INTERPOLATION OPTIONS AND/OR
!!           RESTRICTIONS ON THE INPUT OR OUTPUT GRIDS, BOTH OF WHICH
!!           ARE DOCUMENTED MORE FULLY IN THEIR RESPECTIVE SUBPROGRAMS.
!!           THE GRIDS ARE DEFINED BY THEIR GRID DESCRIPTION SECTIONS
!!           (PASSED IN INTEGER FORM AS DECODED BY SUBPROGRAM W3FI63).
!!           THE CURRENT CODE RECOGNIZES THE FOLLOWING PROJECTIONS:
!!             (KGDS(1)=000) EQUIDISTANT CYLINDRICAL
!!             (KGDS(1)=001) MERCATOR CYLINDRICAL
!!             (KGDS(1)=003) LAMBERT CONFORMAL CONICAL
!!             (KGDS(1)=004) GAUSSIAN CYLINDRICAL
!!             (KGDS(1)=005) POLAR STEREOGRAPHIC AZIMUTHAL
!!             (KGDS(1)=203) ROTATED EQUIDISTANT CYLINDRICAL - E-STAGGER
!!             (KGDS(1)=205) ROTATED EQUIDISTANT CYLINDRICAL - B-STAGGER
!!           WHERE KGDS COULD BE EITHER INPUT KGDSI OR OUTPUT KGDSO.
!!           AS AN ADDED BONUS THE NUMBER OF OUTPUT GRID POINTS
!!           AND THEIR LATITUDES AND LONGITUDES ARE ALSO RETURNED.
!!           ON THE OTHER HAND, THE OUTPUT CAN BE A SET OF STATION POINTS
!!           IF KGDSO(1)<0, IN WHICH CASE THE NUMBER OF POINTS
!!           AND THEIR LATITUDES AND LONGITUDES MUST BE INPUT.
!!           FOR THE BUDGET APPROACH, A SUBSECTION OF THE GRID MAY
!!           BE OUTPUT BY SUBTRACTING KGDSO(1) FROM 255 AND PASSING
!!           IN THE LATITUDES AND LONGITUDES OF THE POINTS.
!!           INPUT BITMAPS WILL BE INTERPOLATED TO OUTPUT BITMAPS.
!!           OUTPUT BITMAPS WILL ALSO BE CREATED WHEN THE OUTPUT GRID
!!           EXTENDS OUTSIDE OF THE DOMAIN OF THE INPUT GRID.
!!           THE OUTPUT FIELD IS SET TO 0 WHERE THE OUTPUT BITMAP IS OFF.
!!        
!!   INPUT ARGUMENT LIST:
!! @param IP       - INTEGER INTERPOLATION METHOD
!!                (IP=0 FOR BILINEAR;
!!                 IP=1 FOR BICUBIC;
!!                 IP=2 FOR NEIGHBOR;
!!                 IP=3 FOR BUDGET;
!!                 IP=4 FOR SPECTRAL;
!!                 IP=6 FOR NEIGHBOR-BUDGET)
!! @param IPOPT    - INTEGER (20) INTERPOLATION OPTIONS
!!                (IP=0: (NO OPTIONS)
!!                 IP=1: CONSTRAINT OPTION
!!                 IP=2: (NO OPTIONS)
!!                 IP=3: NUMBER IN RADIUS, RADIUS WEIGHTS, SEARCH RADIUS
!!                 IP=4: SPECTRAL SHAPE, SPECTRAL TRUNCATION
!!                 IP=6: NUMBER IN RADIUS, RADIUS WEIGHTS ...)
!! @param KGDSI    - INTEGER (200) INPUT GDS PARAMETERS AS DECODED BY W3FI63
!! @param KGDSO    - INTEGER (200) OUTPUT GDS PARAMETERS
!! @param MI       - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS IF KM>1
!!                OR DIMENSION OF INPUT GRID FIELDS IF KM=1
!! @param MO       - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS IF KM>1
!!                OR DIMENSION OF OUTPUT GRID FIELDS IF KM=1
!! @param KM       - INTEGER NUMBER OF FIELDS TO INTERPOLATE
!! @param IBI      - INTEGER (KM) INPUT BITMAP FLAGS
!! @param LI       - LOGICAL*1 (MI,KM) INPUT BITMAPS (IF RESPECTIVE IBI(K)=1)
!! @param GI       - REAL (MI,KM) INPUT FIELDS TO INTERPOLATE
!! @param[out] NO       - INTEGER NUMBER OF OUTPUT POINTS (ONLY IF KGDSO(1)<0)
!! @param[out] RLAT     - REAL (NO) OUTPUT LATITUDES IN DEGREES (IF KGDSO(1)<0)
!! @param[out] RLON     - REAL (NO) OUTPUT LONGITUDES IN DEGREES (IF KGDSO(1)<0)
!! @param[out] IBO      - INTEGER (KM) OUTPUT BITMAP FLAGS
!! @param[out] LO       - LOGICAL*1 (MO,KM) OUTPUT BITMAPS (ALWAYS OUTPUT)
!! @param[out] GO       - REAL (MO,KM) OUTPUT FIELDS INTERPOLATED
!! @param[out] IRET     - INTEGER RETURN CODE
!!                0    SUCCESSFUL INTERPOLATION
!!                1    UNRECOGNIZED INTERPOLATION METHOD
!!                2    UNRECOGNIZED INPUT GRID OR NO GRID OVERLAP
!!                3    UNRECOGNIZED OUTPUT GRID
!!                1X   INVALID BICUBIC METHOD PARAMETERS
!!                3X   INVALID BUDGET METHOD PARAMETERS
!!                4X   INVALID SPECTRAL METHOD PARAMETERS
!!
!! SUBPROGRAMS CALLED:
!!   POLATES0     INTERPOLATE SCALAR FIELDS (BILINEAR)
!!   POLATES1     INTERPOLATE SCALAR FIELDS (BICUBIC)
!!   POLATES2     INTERPOLATE SCALAR FIELDS (NEIGHBOR)
!!   POLATES3     INTERPOLATE SCALAR FIELDS (BUDGET)
!!   POLATES4     INTERPOLATE SCALAR FIELDS (SPECTRAL)
!!
!! REMARKS: EXAMPLES DEMONSTRATING RELATIVE CPU COSTS.
!!   THIS EXAMPLE IS INTERPOLATING 12 LEVELS OF TEMPERATURES
!!   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
!!   TO THE 93 X 68 HAWAIIAN MERCATOR GRID (NCEP GRID 204).
!!   THE EXAMPLE TIMES ARE FOR THE C90.  AS A REFERENCE, THE CP TIME
!!   FOR UNPACKING THE GLOBAL 12 TEMPERATURE FIELDS IS 0.04 SECONDS.
!!
!!   METHOD      |IP  |IPOPT          |CP SECONDS
!!   --------    |--  |-------------  |----------
!!   BILINEAR    |0   |               | 0.03
!!   BICUBIC     |1   |0              | 0.07
!!   BICUBIC     |1   |1              | 0.07
!!   NEIGHBOR    |2   |               | 0.01
!!   BUDGET      |3   |-1,-1          | 0.48
!!   SPECTRAL    |4   |0,40           | 0.22
!!   SPECTRAL    |4   |1,40           | 0.24
!!   SPECTRAL    |4   |0,-1           | 0.42
!!   N-BUDGET    |6   |-1,-1          | 0.15
!!
!!   THE SPECTRAL INTERPOLATION IS FAST FOR THE MERCATOR GRID.
!!   HOWEVER, FOR SOME GRIDS THE SPECTRAL INTERPOLATION IS SLOW.
!!   THE FOLLOWING EXAMPLE IS INTERPOLATING 12 LEVELS OF TEMPERATURES
!!   FROM THE 360 X 181 GLOBAL GRID (NCEP GRID 3)
!!   TO THE 93 X 65 CONUS LAMBERT CONFORMAL GRID (NCEP GRID 211).
!!
!!   METHOD      |IP  |IPOPT          |CP SECONDS
!!   --------    |--  |-------------  |----------
!!   BILINEAR    |0   |               | 0.03
!!   BICUBIC     |1   |0              | 0.07
!!   BICUBIC     |1   |1              | 0.07
!!   NEIGHBOR    |2   |               | 0.01
!!   BUDGET      |3   |-1,-1          | 0.51
!!   SPECTRAL    |4   |0,40           | 3.94
!!   SPECTRAL    |4   |1,40           | 5.02
!!   SPECTRAL    |4   |0,-1           |11.36
!!   N-BUDGET    |6   |-1,-1          | 0.18
!!
SUBROUTINE IPOLATES(IP,IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI, &
     NO,RLAT,RLON,IBO,LO,GO,IRET)
  IMPLICIT NONE
!
 INTEGER,    INTENT(IN   ) :: IP, IPOPT(20), KM, MI, MO
 INTEGER,    INTENT(IN   ) :: IBI(KM), KGDSI(200), KGDSO(200)
 INTEGER,    INTENT(INOUT) :: NO
 INTEGER,    INTENT(  OUT) :: IRET, IBO(KM)
!
 LOGICAL*1,  INTENT(IN   ) :: LI(MI,KM)
 LOGICAL*1,  INTENT(  OUT) :: LO(MO,KM)
!
 REAL,       INTENT(IN   ) :: GI(MI,KM)
 REAL,       INTENT(INOUT) :: RLAT(MO),RLON(MO)
 REAL,       INTENT(  OUT) :: GO(MO,KM)
!
 INTEGER                   :: K, N
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  BILINEAR INTERPOLATION
 IF(IP.EQ.0) THEN
   CALL POLATES0(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  BICUBIC INTERPOLATION
 ELSEIF(IP.EQ.1) THEN
   CALL POLATES1(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  NEIGHBOR INTERPOLATION
 ELSEIF(IP.EQ.2) THEN
   CALL POLATES2(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  BUDGET INTERPOLATION
 ELSEIF(IP.EQ.3) THEN
   CALL POLATES3(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SPECTRAL INTERPOLATION
 ELSEIF(IP.EQ.4) THEN
   CALL POLATES4(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  NEIGHBOR-BUDGET INTERPOLATION
 ELSEIF(IP.EQ.6) THEN
   CALL POLATES6(IPOPT,KGDSI,KGDSO,MI,MO,KM,IBI,LI,GI,NO,RLAT,RLON,IBO,LO,GO,IRET)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  UNRECOGNIZED INTERPOLATION METHOD
 ELSE
   IF(KGDSO(1).GE.0) NO=0
   DO K=1,KM
     IBO(K)=1
     DO N=1,NO
       LO(N,K)=.FALSE.
       GO(N,K)=0.
     ENDDO
   ENDDO
   IRET=1
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 END SUBROUTINE IPOLATES
