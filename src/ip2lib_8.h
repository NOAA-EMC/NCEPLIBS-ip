#ifndef IPLIB
#define IPLIB

void gdswzd(long igdtnum, long *igdtmpl, long igdtlen, long iopt,
            long npts, double fill, double *xpts, double *ypts, 
            double *rlon, double *rlat, long *nret,
            double *crot, double *srot, double *xlon, double *xlat,
            double *ylon, double *ylat, double *area);

void gdswzd_grib1(long kgds, long iopt, long npts, double *fill, double *xpts,
		  double *ypts, double *rlon, double *rlat, long nret, double *crot,
		  double *srot, double *xlon, double *xlat, double *ylon,
		  double *ylat, double *area);
#endif
