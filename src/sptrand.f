C> @file
C> @brief Perform a gradient spherical transform.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | IREDELL | Initial
C> 1998-12-15 | IREDELL | openmp directives inserted
C>
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
C> between spectral coefficients of scalar fields
C> and their means and gradients on a global cylindrical grid.
C>
C> The wave-space can be either triangular or rhomboidal.
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
C> if so requested, just a subset of the latitude pairs
C> may be transformed in each invocation of the subprogram.
C>
C> The transforms are all multiprocessed over latitude except
C> the transform from Fourier to spectral is multiprocessed
C> over zonal wavenumber to ensure reproducibility.
C>
C> Transform several fields at a time to improve vectorization.
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
C> @param KMAX number of fields to transform.
C> @param IPRIME longitude index for the prime meridian.
C> (defaults to 1 if IPRIME=0)
C> @param ISKIP skip number between longitudes
C> (defaults to 1 if ISKIP=0)
C> @param JNSKIP skip number between n.h. latitudes from north
C> (defaults to IMAX if JNSKIP=0)
C> @param JSSKIP skip number between s.h. latitudes from south
C> (defaults to -IMAX if JSSKIP=0)
C> @param KWSKIP skip number between wave fields
C> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) if KWSKIP=0)
C> @param KGSKIP skip number between grid fields
C> (defaults to IMAX*JMAX if KGSKIP=0)
C> @param JBEG latitude index (from pole) to begin transform
C> (defaults to 1 if JBEG=0). If JBEG=0 and IDIR<0, wave is zeroed before transform.
C> @param JEND latitude index (from pole) to end transform
C> (defaults to (JMAX+1)/2 if JEND=0)
C> @param JCPU number of cpus over which to multiprocess
C> @param[out] WAVE wave fields if IDIR>0
C> @param[out] GRIDMN global means if IDIR<0
C> @param[out] GRIDXN n.h. x-gradients (starting at JBEG) if IDIR<0
C> @param[out] GRIDXS s.h. x-gradients (starting at JBEG) if IDIR<0
C> [GRIDX=(D(WAVE)/DLAM)/(CLAT*RERTH)]      
C> @param[out] GRIDYN n.h. y-gradients (starting at JBEG) if IDIR<0
C> @param[out] GRIDYS s.h. y-gradients (starting at JBEG) if IDIR<0
C> [GRIDY=(D(WAVE)/DPHI)/RERTH]
C> @param IDIR transform flag
C> (IDIR>0 for wave to grid, IDIR<0 for grid to wave)
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTRAND(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &                   IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,
     &                   JBEG,JEND,JCPU,
     &                   WAVE,GRIDMN,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

      REAL WAVE(*),GRIDMN(KMAX),GRIDXN(*),GRIDXS(*),GRIDYN(*),GRIDYS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)

C  SET PARAMETERS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      KW=KWSKIP
      IF(KW.EQ.0) KW=2*MX

C  TRANSFORM WAVE TO GRID
      IF(IDIR.GT.0) THEN
C$OMP PARALLEL DO PRIVATE(KWS)
        DO K=1,KMAX
          KWS=(K-1)*KW
          GRIDMN(K)=WAVE(KWS+1)/SQRT(2.)
          CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),1)
          WZ(1:2*MX,K)=0.
        ENDDO
        CALL SPTRANV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &               IPRIME,ISKIP,JNSKIP,JSSKIP,MDIM,KGSKIP,
     &               JBEG,JEND,JCPU,
     &               WD,WZ,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)

C  TRANSFORM GRID TO WAVE
      ELSE
C$OMP PARALLEL DO
        DO K=1,KMAX
          WD(1:2*MX,K)=0.
          WZ(1:2*MX,K)=0.
        ENDDO
        CALL SPTRANV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &               IPRIME,ISKIP,JNSKIP,JSSKIP,MDIM,KGSKIP,
     &               JBEG,JEND,JCPU,
     &               WD,WZ,GRIDXN,GRIDXS,GRIDYN,GRIDYS,IDIR)
        IF(JBEG.EQ.0) THEN
C$OMP PARALLEL DO PRIVATE(KWS)
          DO K=1,KMAX
            KWS=(K-1)*KW
            CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),-1)
            WAVE(KWS+1)=GRIDMN(K)*SQRT(2.)
          ENDDO
        ELSE
C$OMP PARALLEL DO PRIVATE(KWS)
          DO K=1,KMAX
            KWS=(K-1)*KW
            CALL SPLAPLAC(IROMB,MAXWV,ENN1,WZ(1,K),WD(1,K),-1)
            WAVE(KWS+1:KWS+2*MX)=WAVE(KWS+1:KWS+2*MX)+WZ(1:2*MX,K)
            WAVE(KWS+1)=GRIDMN(K)*SQRT(2.)
          ENDDO
        ENDIF
      ENDIF
      END
