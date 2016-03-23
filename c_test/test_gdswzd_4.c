#include <stdio.h>
#include <stdlib.h>

#include "iplib.h"

int main()
{
  int kgds[200];
  int iopt, npts, nret;
  float fill;
  float *xpts, *ypts, *rlon, *rlat;
  float *crot, *srot, *xlon, *xlat, *ylon, *ylat, *area;

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

  xpts = (float *) malloc(npts * sizeof(float));
  ypts = (float *) malloc(npts * sizeof(float));
  rlon = (float *) malloc(npts * sizeof(float));
  rlat = (float *) malloc(npts * sizeof(float));
  crot = (float *) malloc(npts * sizeof(float));
  srot = (float *) malloc(npts * sizeof(float));
  xlon = (float *) malloc(npts * sizeof(float));
  xlat = (float *) malloc(npts * sizeof(float));
  ylon = (float *) malloc(npts * sizeof(float));
  ylat = (float *) malloc(npts * sizeof(float));
  area = (float *) malloc(npts * sizeof(float));


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

  printf(" nret = %d \n", nret);
  printf(" rlon[0], rlat[0] = %f %f \n", rlon[0]-360, rlat[0]);
  printf(" rlon[npts], rlat[npts] = %f %f \n", rlon[nret-1], rlat[nret-1]);

  for (int n=0; n<nret; n++) {
    printf(" n = %d crot, srot %f %f \n", n, crot[n], srot[n]);
     printf(" n = %d rlon[n], rlat[n] = %f %f \n", n, rlon[n], rlat[n]);
  }

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
