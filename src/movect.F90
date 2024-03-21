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
 subroutine movect(flat,flon,tlat,tlon,crot,srot)
   implicit none
!
   integer,parameter     :: kd=selected_real_kind(15,45)
!
   real,intent(in) :: flat,flon
   real,intent(in) :: tlat,tlon
   real,intent(out) :: crot,srot
!
   real(KIND=kd),parameter     :: crdlim=0.9999999
   real(KIND=kd),parameter     :: pi=3.14159265358979
   real(KIND=kd),parameter     :: dpr=180./pi
!
   real(KIND=kd)                  :: ctlat,stlat,cflat,sflat
   real(KIND=kd)                  :: cdlon,sdlon,crd
   real(KIND=kd)                  :: srd2rn,str,ctr,sfr,cfr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE COSINE OF THE RADIAL DISTANCE BETWEEN THE POINTS.
   ctlat=cos(tlat/dpr)
   stlat=sin(tlat/dpr)
   cflat=cos(flat/dpr)
   sflat=sin(flat/dpr)
   cdlon=cos((flon-tlon)/dpr)
   sdlon=sin((flon-tlon)/dpr)
   crd=stlat*sflat+ctlat*cflat*cdlon
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE ROTATIONS AT BOTH POINTS WITH RESPECT TO THE GREAT CIRCLE
!  AND COMBINE THEM TO GIVE THE TOTAL VECTOR ROTATION PARAMETERS.
   if(abs(crd).le.crdlim) then
     srd2rn=-1/(1-crd**2)
     str=cflat*sdlon
     ctr=cflat*stlat*cdlon-sflat*ctlat
     sfr=ctlat*sdlon
     cfr=ctlat*sflat*cdlon-stlat*cflat
     crot=real(srd2rn*(ctr*cfr-str*sfr))
     srot=real(srd2rn*(ctr*sfr+str*cfr))
!  USE A DIFFERENT APPROXIMATION FOR NEARLY COINCIDENT POINTS.
!  MOVING VECTORS TO ANTIPODAL POINTS IS AMBIGUOUS ANYWAY.
   else
     crot=real(cdlon)
     srot=real(sdlon*stlat)
   endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 endsubroutine movect
