/** @file
 * @brief C interface to gdswzd() and gdswzd_grib1() functions for 'd'
 * library build.
 * @author NOAA Programmer 
 */

#ifndef IPLIB
#define IPLIB

/**
   GDSWZD in C.

   @param igdtnum
   @param igdtmpl
   @param igdtlen
   @param iopt
   @param npts
   @param fill
   @param *xpts
   @param *ypts
   @param *rlon
   @param *rlat
   @param *nret
   @param *crot
   @param *srot
   @param *xlon
   @param *xlat
   @param *ylon
   @param *ylat
   @param *area
 */
void gdswzd(int igdtnum, int *igdtmpl, int igdtlen, int iopt, 
            int npts, double fill, double *xpts, double *ypts,
            double *rlon, double *rlat, int *nret,
            double *crot, double *srot, double *xlon, double *xlat,
            double *ylon, double *ylat, double *area);

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
void gdswzd_grib1(int kgds, int iopt, int npts, double *fill, double *xpts,
		  double *ypts, double *rlon, double *rlat, int nret, double *crot,
		  double *srot, double *xlon, double *xlat, double *ylon,
		  double *ylat, double *area);
#endif
