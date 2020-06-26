#include <stdio.h>
#include <stdlib.h>

#include "iplib.h"

/**************************************************************
  Unit test to ensure the 'c' wrapper routine for gdswzd
  is working.

  Call gdswzd for a rotated lat/lon grid with "B" stagger
  and print out the corner point lat/lons and the number
  of valid grid points returned.

  Tests the double precision version of gdswzd.
**************************************************************/

int main()
{
  long kgds[200];
  long iopt, npts, nret;
  double fill;
  double *xpts, *ypts, *rlon, *rlat;
  double *crot, *srot, *xlon, *xlat, *ylon, *ylat, *area;

  long im = 251;
  long jm = 201;

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

  gdswzd(kgds, iopt, npts, fill,
         xpts, ypts, rlon, rlat,
         &nret,
         crot, srot, xlon, xlat, ylon, ylat, area);

  printf(" Points returned from gdswzd = %li \n", nret);
  printf(" Expected points returned    = 50451 \n\n");

  printf(" First corner point lat/lon = %f %f \n", rlat[0], rlon[0]-360);
  printf(" Expected lat/lon           = -7.446 -144.139 \n\n");

  printf(" Last corner point lat/lon  = %f %f \n", rlat[nret-1], rlon[nret-1]);
  printf(" Expected lat/lon           = 44.56 14.744 \n");

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
