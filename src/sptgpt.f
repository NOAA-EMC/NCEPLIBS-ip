C> @file
C> @brief Transform spectral scalar to station points.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | Iredell | Initial.
C> 1998-12-15 | Iredell | Openmp directives inserted.
C> 2003-06-30 | Iredell | Use spfftpt().
C>
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
c> from spectral coefficients of scalar quantities
c> to specified sets of station points on the globe.
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
C> (0 for triangular, 1 for rhomboidal)
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
C> (defaults to 1 if NGSKIP=0)
C> @param RLAT station latitudes in degrees
C> @param RLON station longitudes in degrees
C> @param WAVE wave fields
C> @param GP station point sets
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTGPT(IROMB,MAXWV,KMAX,NMAX,
     &                  KWSKIP,KGSKIP,NRSKIP,NGSKIP,
     &                  RLAT,RLON,WAVE,GP)

      REAL RLAT(*),RLON(*),WAVE(*),GP(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(KMAX)
      REAL WTOP(2*(MAXWV+1),KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+3,2,KMAX)
      PARAMETER(PI=3.14159265358979)

C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MXTOP=MAXWV+1
      IDIM=2*MAXWV+3
      KW=KWSKIP
      KG=KGSKIP
      NR=NRSKIP
      NG=NGSKIP
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=NMAX
      IF(NR.EQ.0) NR=1
      IF(NG.EQ.0) NG=1
      MP=0
C$OMP PARALLEL DO
      DO K=1,KMAX
        WTOP(1:2*MXTOP,K)=0
      ENDDO

C  CALCULATE STATION FIELDS
C$OMP PARALLEL DO PRIVATE(RADLAT,SLAT1,CLAT1)
C$OMP&            PRIVATE(PLN,PLNTOP,F,NK)
      DO N=1,NMAX
        RADLAT=PI/180*RLAT((N-1)*NR+1)
        IF(RLAT((N-1)*NR+1).GE.89.9995) THEN
          SLAT1=1.
          CLAT1=0.
        ELSEIF(RLAT((N-1)*NR+1).LE.-89.9995) THEN
          SLAT1=-1.
          CLAT1=0.
        ELSE
          SLAT1=SIN(RADLAT)
          CLAT1=COS(RADLAT)
        ENDIF
        CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,KW,2*MXTOP,KMAX,
     &               CLAT1,PLN,PLNTOP,MP,WAVE,WTOP,F)
        CALL SPFFTPT(MAXWV,1,2*MAXWV+3,KG,KMAX,RLON((N-1)*NR+1),
     &               F,GP((N-1)*NG+1))
      ENDDO
      END
