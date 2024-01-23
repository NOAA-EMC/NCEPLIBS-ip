C> @file
C>
C> Spectrally interpolate scalars to stations
C> @author IREDELL @date 96-02-29

C> This subprogram spectrally truncates scalar fields on a global
C> cylindrical grid, returning the fields to specified sets of
C> station points on the globe. The wave-space can be either
C> triangular or rhomboidal. The grid-space can be either an
C> equally-spaced grid (with or without pole points) or a Gaussian
C> grid. The grid and point fields may have general indexing. The
C> transforms are all multiprocessed. Transform several fields at a
C> time to improve vectorization. Subprogram can be called from a
C> multiprocessing environment.
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
C> @param KMAX     - INTEGER NUMBER OF FIELDS TO TRANSFORM.
C> @param NMAX     - INTEGER NUMBER OF STATION POINTS TO RETURN
C> @param IPRIME   - INTEGER INPUT LONGITUDE INDEX FOR THE PRIME MERIDIAN.
C>                (DEFAULTS TO 1 IF IPRIME=0)
C>                (OUTPUT LONGITUDE INDEX FOR PRIME MERIDIAN ASSUMED 1.)
C> @param ISKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT LONGITUDES
C>                (DEFAULTS TO 1 IF ISKIPI=0)
C> @param JSKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT LATITUDES FROM SOUTH
C>                (DEFAULTS TO -IMAXI IF JSKIPI=0)
C> @param KSKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS
C>                (DEFAULTS TO IMAXI*JMAXI IF KSKIPI=0)
C> @param KGSKIP   - INTEGER SKIP NUMBER BETWEEN STATION POINT SETS
C>                (DEFAULTS TO NMAX IF KGSKIP=0)
C> @param NRSKIP   - INTEGER SKIP NUMBER BETWEEN STATION LATS AND LONS
C>                (DEFAULTS TO 1 IF NRSKIP=0)
C> @param NGSKIP   - INTEGER SKIP NUMBER BETWEEN STATION POINTS
C>                (DEFAULTS TO 1 IF NGSKIP=0)
C> @param RLAT     - REAL (*) STATION LATITUDES IN DEGREES
C> @param RLON     - REAL (*) STATION LONGITUDES IN DEGREES
C> @param JCPU     - INTEGER NUMBER OF CPUS OVER WHICH TO MULTIPROCESS
C>                (DEFAULTS TO ENVIRONMENT NCPUS IF JCPU=0)
C> @param GRIDI    - REAL (*) INPUT GRID FIELDS
C> @param[out] GP  - REAL (*) STATION POINT SETS
C>
C> SUBPROGRAMS CALLED:
C> - sptran()       Perform a scalar spherical transform
C> - sptgpt()       Transform spectral scalar to station points
C> - ncpus()        Gets environment number of cpus
C>
C> Minimum grid dimensions for unaliased transforms to spectral:
C>   DIMENSION                    |LINEAR              |QUADRATIC
C>   -----------------------      |---------           |-------------
C>   IMAX                         | 2*MAXWV+2          | 3*MAXWV/2*2+2
C>   JMAX (IDRT=4,IROMB=0)        | 1*MAXWV+1          | 3*MAXWV/2+1
C>   JMAX (IDRT=4,IROMB=1)        | 2*MAXWV+1          | 5*MAXWV/2+1
C>   JMAX (IDRT=0,IROMB=0)        | 2*MAXWV+3          | 3*MAXWV/2*2+3
C>   JMAX (IDRT=0,IROMB=1)        | 4*MAXWV+3          | 5*MAXWV/2*2+3
C>   JMAX (IDRT=256,IROMB=0)      | 2*MAXWV+1          | 3*MAXWV/2*2+1
C>   JMAX (IDRT=256,IROMB=1)      | 4*MAXWV+1          | 5*MAXWV/2*2+1
C>
      SUBROUTINE SPTRUNG(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KMAX,NMAX,
     &     IPRIME,ISKIPI,JSKIPI,KSKIPI,KGSKIP,
     &     NRSKIP,NGSKIP,JCPU,RLAT,RLON,GRIDI,GP)

      REAL RLAT(*),RLON(*),GRIDI(*),GP(*)
      REAL W((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
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
      CALL SPTRAN(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KMAX,
     &     IPRIME,ISKIPI,JN,JS,MDIM,KSKIPI,0,0,JC,
     &     W,GRIDI(INP),GRIDI(ISP),-1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO OUTPUT
      CALL SPTGPT(IROMB,MAXWV,KMAX,NMAX,MDIM,KGSKIP,NRSKIP,NGSKIP,
     &     RLAT,RLON,W,GP)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
