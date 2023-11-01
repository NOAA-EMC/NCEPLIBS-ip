C> @file
C>
C> Spectrally interpolate scalars to Mercator
C> @author IREDELL @date 96-02-29

C> This subprogram spectrally truncates scalar fields on a global
C> cylindrical grid, returning the fields to a Mercator grid. The
C> wave-space can be either triangular or rhomboidal. The grid-space
C> can be either an equally-spaced grid (with or without pole
C> points) or a Gaussian grid. The grid fields may have general
C> indexing. The transforms are all multiprocessed. Transform
C> several fields at a time to improve vectorization. Subprogram can
C> be called from a multiprocessing environment.
C> 
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
C> @param MI       - INTEGER NUMBER OF POINTS IN THE FASTER ZONAL DIRECTION
C> @param MJ       - INTEGER NUMBER OF POINTS IN THE SLOWER MERID DIRECTION
C> @param IPRIME   - INTEGER INPUT LONGITUDE INDEX FOR THE PRIME MERIDIAN.
C>                (DEFAULTS TO 1 IF IPRIME=0)
C>                (OUTPUT LONGITUDE INDEX FOR PRIME MERIDIAN ASSUMED 1.)
C> @param ISKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT LONGITUDES
C>                (DEFAULTS TO 1 IF ISKIPI=0)
C> @param JSKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT LATITUDES FROM SOUTH
C>                (DEFAULTS TO -IMAXI IF JSKIPI=0)
C> @param KSKIPI   - INTEGER SKIP NUMBER BETWEEN INPUT GRID FIELDS
C>                (DEFAULTS TO IMAXI*JMAXI IF KSKIPI=0)
C> @param KGSKIP   - INTEGER SKIP NUMBER BETWEEN GRID FIELDS
C>                (DEFAULTS TO NPS*NPS IF KGSKIP=0)
C> @param NISKIP   - INTEGER SKIP NUMBER BETWEEN GRID I-POINTS
C>                (DEFAULTS TO 1 IF NISKIP=0)
C> @param NJSKIP   - INTEGER SKIP NUMBER BETWEEN GRID J-POINTS
C>                (DEFAULTS TO NPS IF NJSKIP=0)
C> @param JCPU     - INTEGER NUMBER OF CPUS OVER WHICH TO MULTIPROCESS
C>                (DEFAULTS TO ENVIRONMENT NCPUS IF JCPU=0)
C> @param RLAT1    - REAL LATITUDE OF THE FIRST GRID POINT IN DEGREES
C> @param RLON1    - REAL LONGITUDE OF THE FIRST GRID POINT IN DEGREES
C> @param DLAT     - REAL LATITUDE INCREMENT IN DEGREES SUCH THAT
C>                D(PHI)/D(J)=DLAT*COS(PHI) WHERE J IS MERIDIONAL INDEX.
C>                DLAT IS NEGATIVE FOR GRIDS INDEXED SOUTHWARD.
C>                (IN TERMS OF GRID INCREMENT DY VALID AT LATITUDE RLATI,
C>                 THE LATITUDE INCREMENT DLAT IS DETERMINED AS
C>                 DLAT=DPR*DY/(RERTH*COS(RLATI/DPR))
C>                 WHERE DPR=180/PI AND RERTH IS EARTH'S RADIUS)
C> @param DLON     - REAL LONGITUDE INCREMENT IN DEGREES SUCH THAT
C>                D(LAMBDA)/D(I)=DLON WHERE I IS ZONAL INDEX.
C>                DLON IS NEGATIVE FOR GRIDS INDEXED WESTWARD.
C> @param GRIDI    - REAL (*) INPUT GRID FIELDS
C> @param GM       - REAL (*) MERCATOR FIELDS
C>
C> SUBPROGRAMS CALLED:
C>   - sptran()       Perform a scalar spherical transform
C>   - sptgpm()       Transform spectral scalar to Mercator
C>   - ncpus()        Gets environment number of cpus
C>
C> MINIMUM GRID DIMENSIONS FOR UNALIASED TRANSFORMS TO SPECTRAL:
C> DIMENSION                    |LINEAR             |QUADRATIC
C> -----------------------      |---------          |-------------
C> IMAX                         | 2*MAXWV+2         | 3*MAXWV/2*2+2
C> JMAX (IDRT=4,IROMB=0)        | 1*MAXWV+1         | 3*MAXWV/2+1
C> JMAX (IDRT=4,IROMB=1)        | 2*MAXWV+1         | 5*MAXWV/2+1
C> JMAX (IDRT=0,IROMB=0)        | 2*MAXWV+3         | 3*MAXWV/2*2+3
C> JMAX (IDRT=0,IROMB=1)        | 4*MAXWV+3         | 5*MAXWV/2*2+3
C> JMAX (IDRT=256,IROMB=0)      | 2*MAXWV+1         | 3*MAXWV/2*2+1
C> JMAX (IDRT=256,IROMB=1)      | 4*MAXWV+1         | 5*MAXWV/2*2+1
C>
      SUBROUTINE SPTRUNM(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KMAX,MI,MJ,
     &                   IPRIME,ISKIPI,JSKIPI,KSKIPI,KGSKIP,
     &                   NISKIP,NJSKIP,JCPU,RLAT1,RLON1,DLAT,DLON,
     &                   GRIDI,GM)
      REAL GRIDI(*),GM(*)
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
     &            IPRIME,ISKIPI,JN,JS,MDIM,KSKIPI,0,0,JC,
     &            W,GRIDI(INP),GRIDI(ISP),-1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO OUTPUT
      CALL SPTGPM(IROMB,MAXWV,KMAX,MI,MJ,MDIM,KGSKIP,NISKIP,NJSKIP,
     &            RLAT1,RLON1,DLAT,DLON,W,GM)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
