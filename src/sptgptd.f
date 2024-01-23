C> @file
C> @brief Transform spectral to station point gradients.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | Iredell | Initial.
C> 1998-12-15 | Iredell | Openmp directives inserted.
C>
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
c> from spectral coefficients of scalar fields
c> to specified sets of station point gradients on the globe.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The wave and point fields may have general indexing,
c> but each wave field is in sequential 'IBM order',
c> i.e. with zonal wavenumber as the slower index.
C>
C> The transforms are all multiprocessed over stations.
C>
C> Transform several fields at a time to improve vectorization.
C>
C> Subprogram can be called from a multiprocessing environment.
C>
C> @param IROMB spectral domain shape
c> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param KMAX number of fields to transform.
C> @param NMAX number of station points to return
C> @param KWSKIP skip number between wave fields
C> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) if KWSKIP=0)
C> @param KGSKIP skip number between station point sets
C> (defaults to NMAX if KGSKIP=0)
C> @param NRSKIP skip number between station lats and lons
C> (defaults to 1 if NRSKIP=0)
C> @param NGSKIP skip number between station points
c> (defaults to 1 if NGSKIP=0)
C> @param RLAT station latitudes in degrees
C> @param RLON station longitudes in degrees
C> @param WAVE wave fields
C> @param XP station point x-gradient sets
C> @param YP station point y-gradient sets
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTGPTD(IROMB,MAXWV,KMAX,NMAX,
     &                   KWSKIP,KGSKIP,NRSKIP,NGSKIP,
     &                   RLAT,RLON,WAVE,XP,YP)

      REAL RLAT(*),RLON(*),WAVE(*),XP(*),YP(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      KW=KWSKIP
      IF(KW.EQ.0) KW=2*MX
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE STATION FIELDS
C$OMP PARALLEL DO PRIVATE(KWS)
      DO K=1,KMAX
        KWS=(K-1)*KW
        CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),1)
        WZ(1:2*MX,K)=0.
      ENDDO
      CALL SPTGPTV(IROMB,MAXWV,KMAX,NMAX,MDIM,KGSKIP,NRSKIP,NGSKIP,
     &             RLAT,RLON,WD,WZ,XP,YP)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
