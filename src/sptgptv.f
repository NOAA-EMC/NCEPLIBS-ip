C> @file
C> @brief Transform spectral vector to station points.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | IREDELL | Initial
C> 1998-12-15 | IREDELL | Openmp directives inserted
C> 1999-08-18 | IREDELL | Openmp directive typo fixed 
C> 2003-06-30 | IREDELL | use spfftpt()
C>
C> @author IREDELL @date 96-02-29

C> This subprogram performs a spherical transform
C> from spectral coefficients of divergences and curls
C> to specified sets of station point vectors on the globe.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The wave and point fields may have general indexing,
C> but each wave field is in sequential 'IBM order',
C> i.e. with zonal wavenumber as the slower index.
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
c> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C> @param KGSKIP skip number between station point sets
c> (defaults to NMAX IF KGSKIP=0)
C> @param NRSKIP skip number between station lats and lons
c> (defaults to 1 if NRSKIP=0)
C> @param NGSKIP skip number between station points
c> (defaults to 1 if NGSKIP=0)
C> @param RLAT station latitudes in degrees
C> @param RLON station longitudes in degrees
C> @param WAVED wave divergence fields
C> @param WAVEZ wave vorticity fields
C> @param UP station point u-wind sets
C> @param VP station point v-wind sets
C>
C> @author IREDELL @date 96-02-29
      SUBROUTINE SPTGPTV(IROMB,MAXWV,KMAX,NMAX,
     &                   KWSKIP,KGSKIP,NRSKIP,NGSKIP,
     &                   RLAT,RLON,WAVED,WAVEZ,UP,VP)

      REAL RLAT(*),RLON(*),WAVED(*),WAVEZ(*),UP(*),VP(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(2*KMAX)
      REAL W((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,2*KMAX)
      REAL WTOP(2*(MAXWV+1),2*KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+3,2,2*KMAX)
      REAL G(2*KMAX)
      PARAMETER(PI=3.14159265358979)

C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MXTOP=MAXWV+1
      MDIM=2*MX+1
      IDIM=2*MAXWV+3
      KW=KWSKIP
      KG=KGSKIP
      NR=NRSKIP
      NG=NGSKIP
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=NMAX
      IF(NR.EQ.0) NR=1
      IF(NG.EQ.0) NG=1
      MP=1

C  CALCULATE SPECTRAL WINDS
C$OMP PARALLEL DO PRIVATE(KWS)
      DO K=1,KMAX
        KWS=(K-1)*KW
        CALL SPDZ2UV(IROMB,MAXWV,ENN1,ELONN1,EON,EONTOP,
     &               WAVED(KWS+1),WAVEZ(KWS+1),
     &               W(1,K),W(1,KMAX+K),WTOP(1,K),WTOP(1,KMAX+K))
      ENDDO

C  CALCULATE STATION FIELDS
C$OMP PARALLEL DO PRIVATE(KU,KV,RADLAT,SLAT1,CLAT1)
C$OMP&            PRIVATE(PLN,PLNTOP,F,G,NK)
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
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,MDIM,2*MXTOP,2*KMAX,
     &               CLAT1,PLN,PLNTOP,MP,W,WTOP,F)
        CALL SPFFTPT(MAXWV,1,2*MAXWV+3,1,2*KMAX,RLON((N-1)*NR+1),F,G)
        DO K=1,KMAX
          KU=K
          KV=K+KMAX
          NK=(N-1)*NG+(K-1)*KG+1
          UP(NK)=G(KU)
          VP(NK)=G(KV)
        ENDDO
      ENDDO

      END
