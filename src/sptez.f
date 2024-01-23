C> @file
C> @brief Perform a simple scalar spherical transform.
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
C> between spectral coefficients of a scalar quantity
C> and a field on a global cylindrical grid.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a Gaussian grid.
C>
C> The wave field is in sequential 'IBM ORDER'.
C>
C> The grid field is indexed East to West, then North to South.
C>
C> For more flexibility and efficiency, call sptran().
C>
C> Subprogram can be called from a multiprocessing environment.
C>
C> Minimum grid dimensions for unaliased transforms to spectral:
C> DIMENSION                    |LINEAR              |QUADRATIC
C> -----------------------      |---------           |-------------
C> IMAX                         |2*MAXWV+2           |3*MAXWV/2*2+2
C> JMAX (IDRT=4,IROMB=0)        |1*MAXWV+1           |3*MAXWV/2+1
C> JMAX (IDRT=4,IROMB=1)        |2*MAXWV+1           |5*MAXWV/2+1
C> JMAX (IDRT=0,IROMB=0)        |2*MAXWV+3           |3*MAXWV/2*2+3
C> JMAX (IDRT=0,IROMB=1)        |4*MAXWV+3           |5*MAXWV/2*2+3
C> JMAX (IDRT=256,IROMB=0)      |2*MAXWV+1           |3*MAXWV/2*2+1
C> JMAX (IDRT=256,IROMB=1)      |4*MAXWV+1           |5*MAXWV/2*2+1
C>
C> @param IROMB spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param IDRT grid identifier
C> - IDRT=4 for Gaussian grid,
C> - IDRT=0 for equally-spaced grid including poles,
C> - IDRT=256 for equally-spaced grid excluding poles
C> @param IMAX even number of longitudes.
C> @param JMAX number of latitudes.
C> @param[out] WAVE wave field if IDIR>0
C> where MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
C> @param[out] GRID grid field (E->W,N->S) if IDIR<0
C> @param IDIR transform flag
C> (IDIR>0 for wave to grid, IDIR<0 for grid to wave).
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTEZ(IROMB,MAXWV,IDRT,IMAX,JMAX,WAVE,GRID,IDIR)

      REAL WAVE((MAXWV+1)*((IROMB+1)*MAXWV+2))
      REAL GRID(IMAX,JMAX)

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
!	print *, " EM: SPTEZ:::JJJJJJJJJJJJJJJJJJJCCCCCCCCCCC=" ,JC	
      IF(IDIR.LT.0) WAVE=0

      CALL SPTRANF(IROMB,MAXWV,IDRT,IMAX,JMAX,1,
     &             IP,IS,JN,JS,KW,KG,JB,JE,JC,
     &             WAVE,GRID,GRID(1,JMAX),IDIR)

      END
