C> @file
C> @brief Perform a vector spherical transform.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | IREDELL | Initial.
C> 1998-12-15 | IREDELL | Generic fft used, openmp directives inserted
C>
C> @author IREDELL @date 96-02-29

C> This subprogram performs a spherical transform
C> between spectral coefficients of divergences and curls
C> and vector fields on a global cylindrical grid.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a Gaussian grid.
C> the wave and grid fields may have general indexing,
C> but each wave field is in sequential 'ibm order',
C> i.e. with zonal wavenumber as the slower index.
C>
C> Transforms are done in latitude pairs for efficiency;
C> thus grid arrays for each hemisphere must be passed.
C> If so requested, just a subset of the latitude pairs
C> may be transformed in each invocation of the subprogram.
C>
C> The transforms are all multiprocessed over latitude except
C> the transform from fourier to spectral is multiprocessed
C> over zonal wavenumber to ensure reproducibility.
C>
C> Transform several fields at a time to improve vectorization.
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
C> @param IROMB SPECTRAL DOMAIN SHAPE
C> (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C> @param MAXWV SPECTRAL TRUNCATION
C> @param IDRT  GRID IDENTIFIER
C> - IDRT=4 FOR GAUSSIAN GRID,
C> - IDRT=0 FOR EQUALLY-SPACED GRID INCLUDING POLES,
C> - IDRT=256 FOR EQUALLY-SPACED GRID EXCLUDING POLES
C> @param IMAX EVEN NUMBER OF LONGITUDES.
C> @param JMAX NUMBER OF LATITUDES.
C> @param KMAX NUMBER OF FIELDS TO TRANSFORM.
C> @param IPRIME LONGITUDE INDEX FOR THE PRIME MERIDIAN.
C> (DEFAULTS TO 1 IF IPRIME=0)
C> @param ISKIP SKIP NUMBER BETWEEN LONGITUDES
C> (DEFAULTS TO 1 IF ISKIP=0)
C> @param JNSKIP SKIP NUMBER BETWEEN N.H. LATITUDES FROM NORTH
C> (DEFAULTS TO IMAX IF JNSKIP=0)
C> @param JSSKIP SKIP NUMBER BETWEEN S.H. LATITUDES FROM SOUTH
C> (DEFAULTS TO -IMAX IF JSSKIP=0)
C> @param KWSKIP SKIP NUMBER BETWEEN WAVE FIELDS
C> (DEFAULTS TO (MAXWV+1)*((IROMB+1)*MAXWV+2) IF KWSKIP=0)
C> @param KGSKIP SKIP NUMBER BETWEEN GRID FIELDS
C> (DEFAULTS TO IMAX*JMAX IF KGSKIP=0)
C> @param JBEG LATITUDE INDEX (FROM POLE) TO BEGIN TRANSFORM
C> - DEFAULTS TO 1 IF JBEG=0
C> - IF JBEG=0 AND IDIR<0, WAVE IS ZEROED BEFORE TRANSFORM
C> @param JEND LATITUDE INDEX (FROM POLE) TO END TRANSFORM
C> (DEFAULTS TO (JMAX+1)/2 IF JEND=0)
C> @param JCPU NUMBER OF CPUS OVER WHICH TO MULTIPROCESS
C> @param[out] WAVED (*) WAVE DIVERGENCE FIELDS IF IDIR>0
C> [WAVED=(D(GRIDU)/DLAM+D(CLAT*GRIDV)/DPHI)/(CLAT*RERTH)]
C> @param[out] WAVEZ (*) WAVE VORTICITY FIELDS IF IDIR>0
C> [WAVEZ=(D(GRIDV)/DLAM-D(CLAT*GRIDU)/DPHI)/(CLAT*RERTH)]
C> @param[out] GRIDUN N.H. GRID U-WINDS (STARTING AT JBEG) IF IDIR<0
C> @param[out] GRIDUS S.H. GRID U-WINDS (STARTING AT JBEG) IF IDIR<0
C> @param[out] GRIDVN N.H. GRID V-WINDS (STARTING AT JBEG) IF IDIR<0
C> @param[out] GRIDVS S.H. GRID V-WINDS (STARTING AT JBEG) IF IDIR<0
C> @param IDIR TRANSFORM FLAG
C> - IDIR>0 FOR WAVE TO GRID,
C> - IDIR<0 FOR GRID TO WAVE
C>
      SUBROUTINE SPTRANV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &                   IPRIME,ISKIP,JNSKIP,JSSKIP,KWSKIP,KGSKIP,
     &                   JBEG,JEND,JCPU,
     &                   WAVED,WAVEZ,GRIDUN,GRIDUS,GRIDVN,GRIDVS,IDIR)

      REAL WAVED(*),WAVEZ(*),GRIDUN(*),GRIDUS(*),GRIDVN(*),GRIDVS(*)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IDIR.LT.0.AND.JBEG.EQ.0) THEN
        DO K=1,KMAX
          KWS=(K-1)*KW
          WAVED(KWS+1:KWS+2*MX)=0
          WAVEZ(KWS+1:KWS+2*MX)=0
        ENDDO
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SPTRANFV(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &              IP,IS,JN,JS,KW,KG,JB,JE,JC,
     &              WAVED,WAVEZ,GRIDUN,GRIDUS,GRIDVN,GRIDVS,IDIR)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
