C> @file
C> @brief Perform a scalar spherical transform.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | IREDELL | Initial
C> 1998-12-15 | IREDELL | Generic fft used, openmp directives inserted
C>
C> @author IREDELL @date 96-02-29

C> This subprogram performs a spherical transform between spectral
C> coefficients of scalar quantities and fields on a global
C> cylindrical grid.
C>
C> The wave-space can be either triangular or
C> rhomboidal.
C>
C> The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a Gaussian grid.
C>
C> The wave and grid fields may have general indexing,
C> but each wave field is in sequential 'IBM order',
C> i.e. with zonal wavenumber as the slower index.
C>
C> Transforms are done in latitude pairs for efficiency;
C> thus grid arrays for each hemisphere must be passed.
C> If so requested, just a subset of the latitude pairs
C> may be transformed in each invocation of the subprogram.
C>
C> The transforms are all multiprocessed over latitude except
C> the transform from Fourier to spectral is multiprocessed
C> over zonal wavenumber to ensure reproducibility.
C>
C> Transform several fields at a time to improve vectorization.
C> Subprogram can be called from a multiprocessing environment.
C>
C> Minimum grid dimensions for unaliased transforms to spectral:
C> DIMENSION                    |LINEAR              |QUADRATIC
C> -----------------------      |---------           |-------------
C> IMAX                         | 2*MAXWV+2          | 3*MAXWV/2*2+2
C> JMAX (IDRT=4,IROMB=0)        | 1*MAXWV+1          | 3*MAXWV/2+1
C> JMAX (IDRT=4,IROMB=1)        | 2*MAXWV+1          | 5*MAXWV/2+1
C> JMAX (IDRT=0,IROMB=0)        | 2*MAXWV+3          | 3*MAXWV/2*2+3
C> JMAX (IDRT=0,IROMB=1)        | 4*MAXWV+3          | 5*MAXWV/2*2+3
C> JMAX (IDRT=256,IROMB=0)      | 2*MAXWV+1          | 3*MAXWV/2*2+1
C> JMAX (IDRT=256,IROMB=1)      | 4*MAXWV+1          | 5*MAXWV/2*2+1
C>      
C> @param IROMB spectral domain shape
c> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param IDRT grid identifier
C> - IDRT=4 for Gaussian grid,
C> - IDRT=0 for equally-spaced grid including poles,
C> - IDRT=256 for equally-spaced grid excluding poles
C> @param IMAX even number of longitudes.
C> @param JMAX number of latitudes.
C> @param KMAX number of fields to transform.
C> @param IPRIME longitude index for the prime meridian.
C> (defaults to 1 if IPRIME=0)
C> @param ISKIP skip number between longitudes
C> (defaults to 1 if ISKIP=0)
C> @param JNSKIP skip number between n.h. latitudes from north
C> (defaults to imax if JNSKIP=0)
C> @param JSSKIP skip number between s.h. latitudes from south
c> (defaults to -imax if JSSKIP=0)
C> @param KWSKIP skip number between wave fields
c> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C> @param KGSKIP skip number between grid fields
c> (defaults to IMAX*JMAX IF KGSKIP=0)
C> @param JBEG latitude index (from pole) to begin transform
c> (defaults to 1 if JBEG=0)
C> (if JBEG=0 and IDIR<0, wave is zeroed before transform)
C> @param JEND latitude index (from pole) to end transform
c> (defaults to (JMAX+1)/2 IF JEND=0)
C> @param JCPU number of cpus over which to multiprocess
C> @param[out] WAVE wave fields if IDIR>0
C> @param[out] gridn n.h. grid fields (starting at jbeg) if IDIR<0
C> @param[out] grids s.h. grid fields (starting at jbeg) if IDIR<0
C> @param IDIR transform flag
C> (idir>0 for wave to grid, idir<0 for grid to wave)
C>
C> @author IREDELL @date 96-02-29
      SUBROUTINE SPTRAN(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &                  IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,
     &                  JBEG,JEND,JCPU,
     &                  WAVE,GRIDN,GRIDS,IDIR)

      REAL WAVE(*),GRIDN(*),GRIDS(*)

      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      IP=IPRIME
      IS=ISKIP
      JN=JNSKIP
      JS=JSSKIP
      KW=KWSKIP
      KG=KGSKIP
      JB=JBEG
      JE=JEND
      JC=JCPU
      IF(IP.EQ.0) IP=1
      IF(IS.EQ.0) IS=1
      IF(JN.EQ.0) JN=IMAX
      IF(JS.EQ.0) JS=-JN
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=IMAX*JMAX
      IF(JB.EQ.0) JB=1
      IF(JE.EQ.0) JE=(JMAX+1)/2
      IF(JC.EQ.0) JC=NCPUS()

      IF(IDIR.LT.0.AND.JBEG.EQ.0) THEN
        DO K=1,KMAX
          KWS=(K-1)*KW
          WAVE(KWS+1:KWS+2*MX)=0
        ENDDO
      ENDIF

      CALL SPTRANF(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &             IP,IS,JN,JS,KW,KG,JB,JE,JC,
     &             WAVE,GRIDN,GRIDS,IDIR)

      END
