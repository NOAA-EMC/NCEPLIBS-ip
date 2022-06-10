#ifndef IPLIB
#define IPLIB

void gdswzd(int igdtnum, int *igdtmpl, int igdtlen, int iopt,
            int npts, float fill, float *xpts, float *ypts, 
            float *rlon, float *rlat, int *nret,
            float *crot, float *srot, float *xlon, float *xlat,
            float *ylon, float *ylat, float *area);

void gdswzd_grib1(int kgds, int iopt, int npts, float *fill, float *xpts,
		  float *ypts, float *rlon, float *rlat, int nret, float *crot,
		  float *srot, float *xlon, float *xlat, float *ylon,
		  float *ylat, float *area);

#endif
