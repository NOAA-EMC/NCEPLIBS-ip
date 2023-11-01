C> @file
C> @brief Transform spectral vector to station points.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | Iredell | Initial.
C> 1998-12-15 | Iredell | Openmp directives inserted.
C> 1999-08-18 | Iredell | Openmp directive typo fixed.
C>
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
C> from spectral coefficients of divergences and curls
C> to specified sets of station point vectors and their
C> gradients on the globe.
C>
C> <pre>
C> DP=(D(UP)/DLON+D(VP*CLAT)/DLAT)/(R*CLAT)
C> ZP=(D(VP)/DLON-D(UP*CLAT)/DLAT)/(R*CLAT)
C> UXP=D(UP*CLAT)/DLON/(R*CLAT)
C> VXP=D(VP*CLAT)/DLON/(R*CLAT)
C> UYP=D(UP*CLAT)/DLAT/R
C> VYP=D(VP*CLAT)/DLAT/R
C> </pre>
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
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param KMAX number of fields to transform.
C> @param NMAX number of station points to return
C> @param KWSKIP skip number between wave fields
C> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C> @param KGSKIP skip number between station point sets
C> (defaults to NMAX if KGSKIP=0)
C> @param NRSKIP skip number between station lats and lons
C> (defaults to 1 if NRSKIP=0)
C> @param NGSKIP skip number between station points
C> (defaults to 1 if NGSKIP=0)
C> @param RLAT station latitudes in degrees
C> @param RLON station longitudes in degrees
C> @param WAVED wave divergence fields
C> @param WAVEZ wave vorticity fields
C> @param DP station point divergence sets
C> @param ZP station point vorticity sets
C> @param UP station point u-wind sets
C> @param VP station point v-wind sets
C> @param UXP station point u-wind x-gradient sets
C> @param VXP station point v-wind x-gradient sets
C> @param UYP station point u-wind y-gradient sets
C> @param VYP station point v-wind y-gradient sets
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTGPTVD(IROMB,MAXWV,KMAX,NMAX,
     &                    KWSKIP,KGSKIP,NRSKIP,NGSKIP,
     &                    RLAT,RLON,WAVED,WAVEZ,
     &                    DP,ZP,UP,VP,UXP,VXP,UYP,VYP)

      REAL RLAT(*),RLON(*),WAVED(*),WAVEZ(*)
      REAL DP(*),ZP(*),UP(*),VP(*),UXP(*),VXP(*),UYP(*),VYP(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(4*KMAX)
      REAL W((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2,4*KMAX)
      REAL WTOP(2*(MAXWV+1),4*KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+2,2,6*KMAX),G(6*KMAX)
      PARAMETER(PI=3.14159265358979)

C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MXTOP=MAXWV+1
      MDIM=2*MX
      IDIM=2*MAXWV+2
      KW=KWSKIP
      KG=KGSKIP
      NR=NRSKIP
      NG=NGSKIP
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=NMAX
      IF(NR.EQ.0) NR=1
      IF(NG.EQ.0) NG=1
      MP(1:2*KMAX)=0
      MP(2*KMAX+1:4*KMAX)=1

C  CALCULATE SPECTRAL WINDS
C$OMP PARALLEL DO PRIVATE(KWS,KD,KZ,KU,KV)
      DO K=1,KMAX
        KWS=(K-1)*KW
        KD=0*KMAX+K
        KZ=1*KMAX+K
        KU=2*KMAX+K
        KV=3*KMAX+K
        DO I=1,2*MX
          W(I,KD)=WAVED(KWS+I)
          W(I,KZ)=WAVEZ(KWS+I)
        ENDDO
        DO I=1,2*MXTOP
          WTOP(I,KD)=0
          WTOP(I,KZ)=0
        ENDDO
        CALL SPDZ2UV(IROMB,MAXWV,ENN1,ELONN1,EON,EONTOP,
     &               WAVED(KWS+1),WAVEZ(KWS+1),
     &               W(1,KU),W(1,KV),WTOP(1,KU),WTOP(1,KV))
      ENDDO

C  CALCULATE STATION FIELDS
C$OMP PARALLEL DO PRIVATE(KD,KZ,KU,KV,KUX,KVX,SLAT1,CLAT1)
C$OMP&            PRIVATE(PLN,PLNTOP,F,G,NK)
      DO N=1,NMAX
        KU=2*KMAX+1
        KUX=4*KMAX+1
        IF(ABS(RLAT((N-1)*NR+1)).GE.89.9995) THEN
          SLAT1=SIGN(1.,RLAT((N-1)*NR+1))
          CLAT1=0.
        ELSE
          SLAT1=SIN(PI/180*RLAT((N-1)*NR+1))
          CLAT1=COS(PI/180*RLAT((N-1)*NR+1))
        ENDIF
        CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,MDIM,2*MXTOP,4*KMAX,
     &               CLAT1,PLN,PLNTOP,MP,W,WTOP,F)
        CALL SPGRADX(MAXWV,IDIM,2*KMAX,MP(2*KMAX+1),CLAT1,
     &               F(1,1,2*KMAX+1),F(1,1,4*KMAX+1))
        CALL SPFFTPT(MAXWV,1,IDIM,1,6*KMAX,RLON((N-1)*NR+1),F,G)
        DO K=1,KMAX
          KD=0*KMAX+K
          KZ=1*KMAX+K
          KU=2*KMAX+K
          KV=3*KMAX+K
          KUX=4*KMAX+K
          KVX=5*KMAX+K
          NK=(N-1)*NG+(K-1)*KG+1
          DP(NK)=G(KD)
          ZP(NK)=G(KZ)
          UP(NK)=G(KU)
          VP(NK)=G(KV)
          UXP(NK)=G(KUX)
          VXP(NK)=G(KVX)
          UYP(NK)=G(KVX)-CLAT1*G(KZ)
          VYP(NK)=CLAT1*G(KD)-G(KUX)
        ENDDO
      ENDDO
      END
