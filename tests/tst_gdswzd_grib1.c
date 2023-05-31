#include <stdio.h>
#include <stdlib.h>

#include "iplib.h"

/**************************************************************
  Unit test to ensure the 'c' wrapper routine for gdswzd
  is working.
  Call gdswzd for a rotated lat/lon grid with "B" stagger
  and print out the corner point lat/lons and the number
  of valid grid points returned.
  Tests the mixed precision version of gdswzd.
**************************************************************/


int main()
{
  int kgds[200];
  int iopt, npts, nret;
  double fill;
  double *xpts, *ypts, *rlon, *rlat;
  double *crot, *srot, *xlon, *xlat, *ylon, *ylat, *area;

  int im = 251;
  int jm = 201;

  kgds[0]  = 205;
  kgds[1]  = im;
  kgds[2]  = jm;
  kgds[3]  =   -7446;
  kgds[4]  = -144139;
  kgds[5]  = 8;
  kgds[6]  =   54000;
  kgds[7]  = -106000;
  kgds[8]  = 0;
  kgds[9]  = 0;
  kgds[10] = 64;
  kgds[11] =   44560;
  kgds[12] =   14744;

  iopt = 1;
  npts = kgds[1] * kgds[2];
  fill = -9999.0;

  xpts = (double *) malloc(npts * sizeof(double));
  ypts = (double *) malloc(npts * sizeof(double));
  rlon = (double *) malloc(npts * sizeof(double));
  rlat = (double *) malloc(npts * sizeof(double));
  crot = (double *) malloc(npts * sizeof(double));
  srot = (double *) malloc(npts * sizeof(double));
  xlon = (double *) malloc(npts * sizeof(double));
  xlat = (double *) malloc(npts * sizeof(double));
  ylon = (double *) malloc(npts * sizeof(double));
  ylat = (double *) malloc(npts * sizeof(double));
  area = (double *) malloc(npts * sizeof(double));


  for (int j=0; j<jm; j++) {
    for (int i=0; i<im; i++) {
      xpts[j*im+i] = i+1;
      ypts[j*im+i] = j+1;
    }
  }

  nret=0;

  gdswzd_grib1(kgds, iopt, npts, fill,
         xpts, ypts, rlon, rlat,
         &nret,
         crot, srot, xlon, xlat, ylon, ylat, area);

  int expextedPointsReturned = 50451;
  
  printf("Points returned from gdswzd = %d \n", nret);
  printf(" Expected points returned    = %d \n\n", expextedPointsReturned);
  
  if (nret != expextedPointsReturned) {
    exit(1);
  }
  
  double expectedLastCornerLat = 44.539;
  double expectedLastCornerLon = 14.802;

  double actualFirstCornerLat = rlat[0];
  double actualFirstCornerLon = rlon[0] - 360.0;


  double expectedFirstCornerLat = -7.491;
  double expectedFirstCornerLon = -144.134;

  double actualLastCornerLat = rlat[nret-1];
  double actualLastCornerLon = rlon[nret-1];


  double MAX_RELATIVE_DIFF = 0.01;

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

  return 0;
}
