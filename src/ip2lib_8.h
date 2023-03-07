/** @file
 * @brief C interface to gdswzd() and gdswzd_grib1() functions for '8'
 * library build.
 * @author Jovic @date 2016
 * @author NOAA Programmer 
 */

#ifndef IPLIB
#define IPLIB

void gdswzd(long igdtnum, long *igdtmpl, long igdtlen, long iopt,
            long npts, double fill, double *xpts, double *ypts, 
            double *rlon, double *rlat, long *nret,
            double *crot, double *srot, double *xlon, double *xlat,
            double *ylon, double *ylat, double *area);

#endif
