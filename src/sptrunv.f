C> @file
C>
C> Spectrally truncate gridded vector fields
C> @author IREDELL @date 96-02-29

C> This subprogram spectrally truncates vector fields
C> on a global cylindrical grid, returning the fields
C> to a possibly different global cylindrical grid.
C> The wave-space can be either triangular or rhomboidal.
C> Either grid-space can be either an equally-spaced grid
C> (with or without pole points) or a Gaussian grid.
C> The grid fields may have general indexing.
C> The transforms are all multiprocessed.
C> Over zonal wavenumber to ensure reproducibility.
C> Transform several fields at a time to improve vectorization.
C> Subprogram can be called from a multiprocessing environment.
C>
C> PROGRAM HISTORY LOG:
C> -  96-02-29  IREDELL
C> - 1998-12-15  IREDELL  OPENMP DIRECTIVES INSERTED
C>
C> @param IROMB    - INTEGER SPECTRAL DOMAIN SHAPE
C>                (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C> @param MAXWV    - INTEGER SPECTRAL TRUNCATION
C> @param IDRTI    - INTEGER INPUT GRID IDENTIFIER
C>                (IDRTI=4 FOR GAUSSIAN GRID,
C>                 IDRTI=0 FOR EQUALLY-SPACED GRID INCLUDING POLES,
C>                 IDRTI=256 FOR EQUALLY-SPACED GRID EXCLUDING POLES)
C> @param IMAXI    - INTEGER EVEN NUMBER OF INPUT LONGITUDES.
C> @param JMAXI    - INTEGER NUMBER OF INPUT LATITUDES.
C> @param IDRTO    - INTEGER OUTPUT GRID IDENTIFIER
C>                (IDRTO=4 FOR GAUSSIAN GRID,
C>                 IDRTO=0 FOR EQUALLY-SPACED GRID INCLUDING POLES,
C>                 IDRTO=256 FOR EQUALLY-SPACED GRID EXCLUDING POLES)
C> @param IMAXO    - INTEGER EVEN NUMBER OF OUTPUT LONGITUDES.
C> @param JMAXO    - INTEGER NUMBER OF OUTPUT LATITUDES.
C> @param KMAX     - INTEGER NUMBER OF FIELDS TO TRANSFORM.
C> @param IPRIME   - INTEGER INPUT LONGITUDE INDEX FOR THE PRIME MERIDIAN.
C>                (DEFAULTS TO 1 IF IPRIME=0)
C>                (OUTPUT LONGITUDE INDEX FOR PRIME MERIDIAN ASSUMED 1.)
C> @param ISKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT LONGITUDES
C>                (DEFAULTS TO 1 IF ISKIPI=0)
C> @param JSKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT LATITUDES FROM SOUTH
C>                (DEFAULTS TO -IMAXI IF JSKIPI=0)
C> @param KSKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS
C>                (DEFAULTS TO IMAXI*JMAXI IF KSKIPI=0)
C> @param ISKIPO   - INTEGER SKIP NUMBER BETWEEN OUTPUT LONGITUDES
C>                (DEFAULTS TO 1 IF ISKIPO=0)
C> @param JSKIPO   - INTEGER SKIP NUMBER BETWEEN OUTPUT LATITUDES FROM SOUTH
C>                (DEFAULTS TO -IMAXO IF JSKIPO=0)
C> @param KSKIPO   - INTEGER SKIP NUMBER BETWEEN OUTPUT GRID FIELDS
C>                (DEFAULTS TO IMAXO*JMAXO IF KSKIPO=0)
C> @param JCPU     - INTEGER NUMBER OF CPUS OVER WHICH TO MULTIPROCESS
C>                (DEFAULTS TO ENVIRONMENT NCPUS IF JCPU=0)
C> @param GRIDUI   - REAL (*) INPUT GRID U-WINDS
C> @param GRIDVI   - REAL (*) INPUT GRID V-WINDS
C> @param LUV      - LOGICAL FLAG WHETHER TO RETURN WINDS
C> @param LDZ      - LOGICAL FLAG WHETHER TO RETURN DIVERGENCE AND VORTICITY
C> @param LPS      - LOGICAL FLAG WHETHER TO RETURN POTENTIAL AND STREAMFCN
C> @param GRIDUO   - REAL (*) OUTPUT U-WINDS IF LUV
C>                (MAY OVERLAY INPUT FIELDS IF GRID SHAPE IS APPROPRIATE)
C> @param GRIDVO   - REAL (*) OUTPUT V-WINDS IF LUV
C>                (MAY OVERLAY INPUT FIELDS IF GRID SHAPE IS APPROPRIATE)
C> @param GRIDDO   - REAL (*) OUTPUT DIVERGENCES IF LDZ
C>                (MAY OVERLAY INPUT FIELDS IF GRID SHAPE IS APPROPRIATE)
C> @param GRIDZO   - REAL (*) OUTPUT VORTICITIES IF LDZ
C>                (MAY OVERLAY INPUT FIELDS IF GRID SHAPE IS APPROPRIATE)
C> @param GRIDPO   - REAL (*) OUTPUT POTENTIALS IF LPS
C>                (MAY OVERLAY INPUT FIELDS IF GRID SHAPE IS APPROPRIATE)
C> @param GRIDSO   - REAL (*) OUTPUT STREAMFCNS IF LPS
C>                (MAY OVERLAY INPUT FIELDS IF GRID SHAPE IS APPROPRIATE)
C>
C> SUBPROGRAMS CALLED:
C>   - SPWGET()       GET WAVE-SPACE CONSTANTS
C>   - SPLAPLAC()     COMPUTE LAPLACIAN IN SPECTRAL SPACE
C>   - SPTRAN()       PERFORM A SCALAR SPHERICAL TRANSFORM
C>   - SPTRANV()      PERFORM A VECTOR SPHERICAL TRANSFORM
C>   - NCPUS()        GETS ENVIRONMENT NUMBER OF CPUS
C>
C> REMARKS: MINIMUM GRID DIMENSIONS FOR UNALIASED TRANSFORMS TO SPECTRAL:
C>   DIMENSION                    |LINEAR             |QUADRATIC
C>   -----------------------      |---------           |-------------
C>   IMAX                         |2*MAXWV+2           |3*MAXWV/2*2+2
C>   JMAX (IDRT=4,IROMB=0)        |1*MAXWV+1           |3*MAXWV/2+1
C>   JMAX (IDRT=4,IROMB=1)        |2*MAXWV+1           |5*MAXWV/2+1
C>   JMAX (IDRT=0,IROMB=0)        |2*MAXWV+3           |3*MAXWV/2*2+3
C>   JMAX (IDRT=0,IROMB=1)        |4*MAXWV+3           |5*MAXWV/2*2+3
C>   JMAX (IDRT=256,IROMB=0)      |2*MAXWV+1           |3*MAXWV/2*2+1
C>   JMAX (IDRT=256,IROMB=1)      |4*MAXWV+1           |5*MAXWV/2*2+1
      SUBROUTINE SPTRUNV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,
     &                   IDRTO,IMAXO,JMAXO,KMAX,
     &                   IPRIME,ISKIPI,JSKIPI,KSKIPI,
     &                   ISKIPO,JSKIPO,KSKIPO,JCPU,GRIDUI,GRIDVI,
     &                   LUV,GRIDUO,GRIDVO,LDZ,GRIDDO,GRIDZO,
     &                   LPS,GRIDPO,GRIDSO)
      LOGICAL LUV,LDZ,LPS
      REAL GRIDUI(*),GRIDVI(*)
      REAL GRIDUO(*),GRIDVO(*),GRIDDO(*),GRIDZO(*),GRIDPO(*),GRIDSO(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM INPUT GRID TO WAVE
      JC=JCPU
      IF(JC.EQ.0) JC=NCPUS()
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      JN=-JSKIPI
      IF(JN.EQ.0) JN=IMAXI
      JS=-JN
      INP=(JMAXI-1)*MAX(0,-JN)+1
      ISP=(JMAXI-1)*MAX(0,-JS)+1
      CALL SPTRANV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KMAX,
     &             IPRIME,ISKIPI,JN,JS,MDIM,KSKIPI,0,0,JC,
     &             WD,WZ,
     &             GRIDUI(INP),GRIDUI(ISP),GRIDVI(INP),GRIDVI(ISP),-1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO OUTPUT WINDS
      JN=-JSKIPO
      IF(JN.EQ.0) JN=IMAXO
      JS=-JN
      INP=(JMAXO-1)*MAX(0,-JN)+1
      ISP=(JMAXO-1)*MAX(0,-JS)+1
      IF(LUV) THEN
        CALL SPTRANV(IROMB,MAXWV,IDRTO,IMAXO,JMAXO,KMAX,
     &               0,ISKIPO,JN,JS,MDIM,KSKIPO,0,0,JC,
     &               WD,WZ,
     &               GRIDUO(INP),GRIDUO(ISP),GRIDVO(INP),GRIDVO(ISP),1)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO OUTPUT DIVERGENCE AND VORTICITY
      IF(LDZ) THEN
        CALL SPTRAN(IROMB,MAXWV,IDRTO,IMAXO,JMAXO,KMAX,
     &              0,ISKIPO,JN,JS,MDIM,KSKIPO,0,0,JC,
     &              WD,GRIDDO(INP),GRIDDO(ISP),1)
        CALL SPTRAN(IROMB,MAXWV,IDRTO,IMAXO,JMAXO,KMAX,
     &              0,ISKIPO,JN,JS,MDIM,KSKIPO,0,0,JC,
     &              WZ,GRIDZO(INP),GRIDZO(ISP),1)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO OUTPUT POTENTIAL AND STREAMFUNCTION
      IF(LPS) THEN
        CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
C$OMP PARALLEL DO
        DO K=1,KMAX
          CALL SPLAPLAC(IROMB,MAXWV,ENN1,WD(1,K),WD(1,K),-1)
          CALL SPLAPLAC(IROMB,MAXWV,ENN1,WZ(1,K),WZ(1,K),-1)
          WD(1:2,K)=0.
          WZ(1:2,K)=0.
        ENDDO
        CALL SPTRAN(IROMB,MAXWV,IDRTO,IMAXO,JMAXO,KMAX,
     &              0,ISKIPO,JN,JS,MDIM,KSKIPO,0,0,JC,
     &              WD,GRIDPO(INP),GRIDPO(ISP),1)
        CALL SPTRAN(IROMB,MAXWV,IDRTO,IMAXO,JMAXO,KMAX,
     &              0,ISKIPO,JN,JS,MDIM,KSKIPO,0,0,JC,
     &              WZ,GRIDSO(INP),GRIDSO(ISP),1)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
