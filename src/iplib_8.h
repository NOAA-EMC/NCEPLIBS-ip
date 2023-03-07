/** @file
 * @brief C interface to gdswzd() function for '8' library build.
 * @author NOAA Programmer 
 */

#ifndef IPLIB
#define IPLIB

void gdswzd(long *kgds, long iopt, long npts, double fill,
            double *xpts, double *ypts, double *rlon, double *rlat,
            long *nret,
            double *crot, double *srot, double *xlon, double *xlat,
            double *ylon, double *ylat, double *area);

#endif
