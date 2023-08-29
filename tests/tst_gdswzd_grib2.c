/**************************************************************
  Unit test to ensure the 'c' wrapper routine for gdswzd
  is working.
  Call gdswzd for a rotated lat/lon grid with "B" stagger
  and print out the corner point lat/lons and the number
  of valid grid points returned.
  Tests the mixed precision version of gdswzd.
**************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "ip2lib.h"

#if(LSIZE==4)
#define REALTYPE float
#else
#define REALTYPE double
#endif
#if(LSIZE==8)
#define INTTYPE long
#else
#define INTTYPE int
#endif


int main()
{
  INTTYPE *igdtmpl;
  INTTYPE igdtnum, igdtlen, iopt, npts, nret;
  REALTYPE fill;
  REALTYPE *xpts, *ypts, *rlon, *rlat;
  REALTYPE *crot, *srot, *xlon, *xlat, *ylon, *ylat, *area;

  int im = 251;
  int jm = 201;

  igdtnum = 1;

  igdtlen = 22;
  igdtmpl = (INTTYPE *) malloc(igdtlen * sizeof(INTTYPE));

  igdtmpl[0] = 6;
  igdtmpl[1] = 255;
  igdtmpl[2] = -1;
  igdtmpl[3] = 255;
  igdtmpl[4] = -1;
  igdtmpl[5] = 255;
  igdtmpl[6] = -1;
  igdtmpl[7] = im;
  igdtmpl[8] = jm;
  igdtmpl[9] = 0;
  igdtmpl[10] = -1;
  igdtmpl[11] = -45036000;
  igdtmpl[12] = 299961000;
  igdtmpl[13] = 56;
  igdtmpl[14] = 45036000;
  igdtmpl[15] = 60039000;
  igdtmpl[16] = 0;
  igdtmpl[17] = 0;
  igdtmpl[18] = 64;
  igdtmpl[19] = -36000000;
  igdtmpl[20] = 254000000;
  igdtmpl[21] = 0;

  iopt = 1;
  npts = im * jm;
  fill = -9999.0;

  xpts = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  ypts = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  rlon = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  rlat = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  crot = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  srot = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  xlon = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  xlat = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  ylon = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  ylat = (REALTYPE *) malloc(npts * sizeof(REALTYPE));
  area = (REALTYPE *) malloc(npts * sizeof(REALTYPE));


  for (int j=0; j<jm; j++) {
  for (int i=0; i<im; i++) {
     xpts[j*im+i] = i+1;
     ypts[j*im+i] = j+1;
  }
  }

  nret=0;

  gdswzd(igdtnum, igdtmpl, igdtlen, iopt, npts, fill,
         xpts, ypts, rlon, rlat, &nret,
         crot, srot, xlon, xlat, ylon, ylat, area);

  int expextedPointsReturned = 50451;
  
  printf("Points returned from gdswzd = %d \n", nret);
  printf(" Expected points returned    = %d \n\n", expextedPointsReturned);
  
  if (nret != expextedPointsReturned) {
          exit(1);
  }
  
  REALTYPE expectedLastCornerLat = 44.539;
  REALTYPE expectedLastCornerLon = 14.802;

  REALTYPE actualFirstCornerLat = rlat[0];
  REALTYPE actualFirstCornerLon = rlon[0] - 360.0;


  REALTYPE expectedFirstCornerLat = -7.491;
  REALTYPE expectedFirstCornerLon = -144.134;

  REALTYPE actualLastCornerLat = rlat[nret-1];
  REALTYPE actualLastCornerLon = rlon[nret-1];


  REALTYPE MAX_RELATIVE_DIFF = 0.01;

  printf(" First corner point lat/lon = %f %f \n", actualFirstCornerLat, actualFirstCornerLon);
  printf(" Expected lat/lon           = -7.491 -144.134 \n\n");
  
  if ((abs(expectedFirstCornerLat - actualFirstCornerLat) / expectedFirstCornerLat > MAX_RELATIVE_DIFF)) {
          exit(1);
  }

  if (abs(expectedFirstCornerLon - actualFirstCornerLon) / expectedFirstCornerLon > MAX_RELATIVE_DIFF) {
          exit(1);
  }
  

  printf(" Last corner point lat/lon  = %f %f \n", actualLastCornerLat, actualLastCornerLon);
  printf(" Expected lat/lon           = 44.539 14.802 \n");

  if (abs(expectedLastCornerLat - actualLastCornerLat) / expectedLastCornerLat > MAX_RELATIVE_DIFF) {
          exit(1);
  }

  if (abs(expectedLastCornerLon - actualLastCornerLon) / expectedLastCornerLon > MAX_RELATIVE_DIFF) {
          exit(1);
  }

/*
  for (int n=0; n<nret; n++) {
    printf(" n = %d crot, srot %f %f \n", n, crot[n], srot[n]);
     printf(" n = %d rlon[n], rlat[n] = %f %f \n", n, rlon[n], rlat[n]);
  }
*/

  free(xpts);
  free(ypts);
  free(rlon);
  free(rlat);
  free(crot);
  free(srot);
  free(xlon);
  free(xlat);
  free(ylon);
  free(ylat);
  free(area);
  free(igdtmpl);

  return 0;
}
