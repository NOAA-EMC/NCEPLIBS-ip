C> @file
C> @brief Perform a simple vector spherical transform
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
C> between spectral coefficients of divergence and curl
C> and a vector field on a global cylindrical grid.
C> The wave-space can be either triangular or rhomboidal.
C>
C> The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a Gaussian grid.
C>
C> The wave field is in sequential 'IBM order'.
C>
C> The grid field is indexed east to west, then north to south.
C>
C> For more flexibility and efficiency, call SPTRAN().
C>
C> Subprogram can be called from a multiprocessing environment.
C>
C> Minimum grid dimensions for unaliased transforms to spectral:
C> Dimension                    |Linear              |Quadratic
C> -----------------------      |---------           |-------------
C> IMAX                         |2*MAXWV+2           |3*MAXWV/2*2+2
C> JMAX (IDRT=4,IROMB=0)        |1*MAXWV+1           |3*MAXWV/2+1
C> JMAX (IDRT=4,IROMB=1)        |2*MAXWV+1           |5*MAXWV/2+1
C> JMAX (IDRT=0,IROMB=0)        |2*MAXWV+3           |3*MAXWV/2*2+3
C> JMAX (IDRT=0,IROMB=1)        |4*MAXWV+3           |5*MAXWV/2*2+3
C> JMAX (IDRT=256,IROMB=0)      |2*MAXWV+1           |3*MAXWV/2*2+1
C> JMAX (IDRT=256,IROMB=1)      |4*MAXWV+1           |5*MAXWV/2*2+1
C>
C> @param IROMB Spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV Spectral truncation
C> @param IDRT Grid identifier
C> - IDRT=4 for Gaussian grid
C> - IDRT=0 for equally-spaced grid including poles
C> - IDRT=256 for equally-spaced grid excluding poles
C> @param IMAX Even number of longitudes
C> @param JMAX Number of latitudes
C> @param[out] WAVED Wave divergence field if IDIR>0
C> where MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
C> @param[out] WAVEZ Wave vorticity field if IDIR>0
C> where MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
C> @param[out] GRIDU Grid u-wind (E->W,N->S) if IDIR<0
C> @param[out] GRIDV Grid v-wind (E->W,N->S) if IDIR<0
C> @param IDIR Transform flag
C> (IDIR>0 for wave to grid, IDIR<0 for grid to wave)
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTEZV(IROMB,MAXWV,IDRT,IMAX,JMAX,
     &                  WAVED,WAVEZ,GRIDU,GRIDV,IDIR)

      REAL WAVED((MAXWV+1)*((IROMB+1)*MAXWV+2))
      REAL WAVEZ((MAXWV+1)*((IROMB+1)*MAXWV+2))
      REAL GRIDU(IMAX,JMAX)
      REAL GRIDV(IMAX,JMAX)

      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      IP=1
      IS=1
      JN=IMAX
      JS=-JN
      KW=2*MX
      KG=IMAX*JMAX
      JB=1
      JE=(JMAX+1)/2
      JC=NCPUS()
      IF(IDIR.LT.0) WAVED=0
      IF(IDIR.LT.0) WAVEZ=0

      CALL SPTRANFV(IROMB,MAXWV,IDRT,IMAX,JMAX,1,
     &              IP,IS,JN,JS,KW,KG,JB,JE,JC,
     &              WAVED,WAVEZ,
     &              GRIDU,GRIDU(1,JMAX),GRIDV,GRIDV(1,JMAX),IDIR)
      END
