C> @file
C> @brief Transform spectral scalar to station points.
C> @author Iredell @date 96-02-29
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
C> from spectral coefficients of scalar quantities
C> to specified sets of station point values
C> and their gradients on the globe.
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
C> @param XP station point x-gradient sets
C> @param YP station point y-gradient sets
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTGPTSD(IROMB,MAXWV,KMAX,NMAX,
     &                    KWSKIP,KGSKIP,NRSKIP,NGSKIP,
     &                    RLAT,RLON,WAVE,GP,XP,YP)

      REAL RLAT(*),RLON(*),WAVE(*)
      REAL GP(*),XP(*),YP(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(2*KMAX)
      REAL W((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2,2*KMAX)
      REAL WTOP(2*(MAXWV+1),2*KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+2,2,3*KMAX),G(3*KMAX)
      PARAMETER(PI=3.14159265358979)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
      MP(1:KMAX)=10
      MP(KMAX+1:2*KMAX)=1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE SPECTRAL WINDS
C$OMP PARALLEL DO PRIVATE(KWS,KS,KY)
      DO K=1,KMAX
        KWS=(K-1)*KW
        KS=0*KMAX+K
        KY=1*KMAX+K
        DO I=1,2*MX
          W(I,KS)=WAVE(KWS+I)
        ENDDO
        DO I=1,2*MXTOP
          WTOP(I,KS)=0
        ENDDO
        CALL SPGRADY(IROMB,MAXWV,ENN1,EON,EONTOP,
     &               WAVE(KWS+1),W(1,KY),WTOP(1,KY))
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE STATION FIELDS
C$OMP PARALLEL DO PRIVATE(KS,KY,KX,SLAT1,CLAT1)
C$OMP&            PRIVATE(PLN,PLNTOP,F,G,NK)
      DO N=1,NMAX
        IF(ABS(RLAT((N-1)*NR+1)).GE.89.9995) THEN
          SLAT1=SIGN(1.,RLAT((N-1)*NR+1))
          CLAT1=0.
        ELSE
          SLAT1=SIN(PI/180*RLAT((N-1)*NR+1))
          CLAT1=COS(PI/180*RLAT((N-1)*NR+1))
        ENDIF
        CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,MDIM,2*MXTOP,2*KMAX,
     &               CLAT1,PLN,PLNTOP,MP,W,WTOP,F)
        CALL SPGRADX(MAXWV,IDIM,KMAX,MP,CLAT1,F(1,1,1),F(1,1,2*KMAX+1))
        CALL SPFFTPT(MAXWV,1,IDIM,1,3*KMAX,RLON((N-1)*NR+1),F,G)
        DO K=1,KMAX
          KS=0*KMAX+K
          KY=1*KMAX+K
          KX=2*KMAX+K
          NK=(N-1)*NG+(K-1)*KG+1
          GP(NK)=G(KS)
          XP(NK)=G(KX)
          YP(NK)=G(KY)
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
