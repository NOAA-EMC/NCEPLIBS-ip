!> @file
!> @brief Move a vector along a great circle.
!> @author Iredell @date 96-04-10

!> This subprogram provides the rotation parameters to move a vector
!> along a great circle from one position to another while conserving
!> its orientation with respect to the great circle. These rotation
!> parameters are useful for vector interpolation.
!>        
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | Iredell | Initial
!> 1999-04-08 | Iredell | generalize precision
!>
!> @param[in] flat real latitude in degrees from which to move the vector.
!> @param[in] flon real longitude in degrees from which to move the vector.
!> @param[in] tlat real latitude in degrees to which to move the vector.
!> @param[in] tlon real longitude in degrees to which to move the vector.
!> @param[out] crot real clockwise vector rotation cosine.
!> @param[out] srot real clockwise vector rotation sine.
!> (uto=crot*ufrom-srot*vfrom; vto=srot*ufrom+crot*vfrom)
!>
!> @author Iredell @date 96-04-10
 SUBROUTINE MOVECT(FLAT,FLON,TLAT,TLON,CROT,SROT)
 IMPLICIT NONE
!
 INTEGER,         PARAMETER     :: KD=SELECTED_REAL_KIND(15,45)
!
 REAL,            INTENT(IN   ) :: FLAT, FLON
 REAL,            INTENT(IN   ) :: TLAT, TLON
 REAL,            INTENT(  OUT) :: CROT, SROT
!
 REAL(KIND=KD),   PARAMETER     :: CRDLIM=0.9999999
 REAL(KIND=KD),   PARAMETER     :: PI=3.14159265358979
 REAL(KIND=KD),   PARAMETER     :: DPR=180./PI
!
 REAL(KIND=KD)                  :: CTLAT,STLAT,CFLAT,SFLAT
 REAL(KIND=KD)                  :: CDLON,SDLON,CRD
 REAL(KIND=KD)                  :: SRD2RN,STR,CTR,SFR,CFR
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE COSINE OF THE RADIAL DISTANCE BETWEEN THE POINTS.
 CTLAT=COS(TLAT/DPR)
 STLAT=SIN(TLAT/DPR)
 CFLAT=COS(FLAT/DPR)
 SFLAT=SIN(FLAT/DPR)
 CDLON=COS((FLON-TLON)/DPR)
 SDLON=SIN((FLON-TLON)/DPR)
 CRD=STLAT*SFLAT+CTLAT*CFLAT*CDLON
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE ROTATIONS AT BOTH POINTS WITH RESPECT TO THE GREAT CIRCLE
!  AND COMBINE THEM TO GIVE THE TOTAL VECTOR ROTATION PARAMETERS.
 IF(ABS(CRD).LE.CRDLIM) THEN
   SRD2RN=-1/(1-CRD**2)
   STR=CFLAT*SDLON
   CTR=CFLAT*STLAT*CDLON-SFLAT*CTLAT
   SFR=CTLAT*SDLON
   CFR=CTLAT*SFLAT*CDLON-STLAT*CFLAT
   CROT=REAL(SRD2RN*(CTR*CFR-STR*SFR))
   SROT=REAL(SRD2RN*(CTR*SFR+STR*CFR))
!  USE A DIFFERENT APPROXIMATION FOR NEARLY COINCIDENT POINTS.
!  MOVING VECTORS TO ANTIPODAL POINTS IS AMBIGUOUS ANYWAY.
 ELSE
   CROT=REAL(CDLON)
   SROT=REAL(SDLON*STLAT)
 ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 END SUBROUTINE MOVECT
