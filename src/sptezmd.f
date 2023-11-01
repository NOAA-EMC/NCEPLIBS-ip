C> @file
C> @brief Perform simple gradient spherical transforms.
C> @author Iredell @date 96-02-29

C> This subprogram performs spherical transforms
C> between spectral coefficients of scalar fields
C> and their means and gradients on a global cylindrical grid.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a gaussian grid.
C>
C> The wave fields are in sequential 'IBM ORDER'.
C>
C> The grid fields are indexed East to West, then North to South.
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
C> - IDRT=4 for Gaussian grid
C> - IDRT=0 for equally-spaced grid including poles
C> - IDRT=256 for equally-spaced grid excluding poles
C> @param IMAX even number of longitudes.
C> @param JMAX number of latitudes.
C> @param KMAX number
C> @param[out] WAVE wave field if IDIR>0
C> where MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)
C> @param[out] GRIDMN global mean if IDIR<0
C> @param[out] GRIDX grid x-gradients (E->W,N->S) if IDIR<0
C> @param[out] GRIDY grid y-gradients (E->W,N->S) if IDIR<0
C> @param IDIR transform flag
C> (IDIR>0 for wave to grid, IDIR<0 for grid to wave).
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTEZMD(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &                   WAVE,GRIDMN,GRIDX,GRIDY,IDIR)

      REAL WAVE((MAXWV+1)*((IROMB+1)*MAXWV+2),KMAX)
      REAL GRIDMN(KMAX),GRIDX(IMAX,JMAX,KMAX),GRIDY(IMAX,JMAX,KMAX)

      JC=NCPUS()
      CALL SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &             0,0,0,0,0,0,0,0,JC,
     &             WAVE,GRIDMN,
     &             GRIDX,GRIDX(1,JMAX,1),GRIDY,GRIDY(1,JMAX,1),IDIR)
      END
