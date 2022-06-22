/** @file
 * @brief C interface to gdswzd() function for '4' library build.
 * @author NOAA Programmer 
 */

#ifndef IPLIB
#define IPLIB

/**
   GDSWZD_grib1 in C.

   @param kgds
   @param iopt
   @param npts
   @param *fill
   @param *xpts
   @param *ypts
   @param *rlon
   @param *rlat
   @param nret
   @param *crot
   @param *srot
   @param *xlon
   @param *xlat
   @param *ylon
   @param *ylat
   @param *area
 */
void gdswzd(int *kgds, int iopt, int npts, float fill,
            float *xpts, float *ypts, float *rlon, float *rlat,
            int *nret,
            float *crot, float *srot, float *xlon, float *xlat,
            float *ylon, float *ylat, float *area);

#endif
