C> @file
C> @brief Truncate gridded scalar fields
C> @author IREDELL @date 96-02-29

C> This subprogram spectrally truncates scalar fields on a global
C> cylindrical grid, returning the fields to a possibly different
C> global cylindrical grid. The wave-space can be either triangular
C> or rhomboidal. either grid-space can be either an equally-spaced
C> grid (with or without pole points) or a Gaussian grid. the grid
C> fields may have general indexing. the transforms are all
C> multiprocessed. Transform several fields at a time to improve
C> vectorization. Subprogram can be called from a multiprocessing
C> environment.
C>
C> Remarks: Minimum grid dimensions for unaliased transforms to spectral:
C>   Dimension                 |  Linear       |     Quadratic
C>   -----------------------   |  ---------    |     -------------
C>   IMAX                      |  2*MAXWV+2    |     3*MAXWV/2*2+2
C>   JMAX (IDRT=4,IROMB=0)     |  1*MAXWV+1    |     3*MAXWV/2+1
C>   JMAX (IDRT=4,IROMB=1)     |  2*MAXWV+1    |     5*MAXWV/2+1
C>   JMAX (IDRT=0,IROMB=0)     |  2*MAXWV+3    |     3*MAXWV/2*2+3
C>   JMAX (IDRT=0,IROMB=1)     |  4*MAXWV+3    |     5*MAXWV/2*2+3
C>   JMAX (IDRT=256,IROMB=0)   |  2*MAXWV+1    |     3*MAXWV/2*2+1
C>   JMAX (IDRT=256,IROMB=1)   |  4*MAXWV+1    |     5*MAXWV/2*2+1
C>
C> @param IROMB Spectral domain shape (0 for triangular, 1 for rhomboidal)
C> @param MAXWV Spectral truncation
C> @param IDRTI Input grid identifier
C>  - IDRTI=4 for Gaussian grid
C>  - IDRTI=0 for equally-spaced grid including poles
C>  - IDRTI=256 for equally-spaced grid excluding poles
C> @param IMAXI Even number of input longitudes
C> @param JMAXI Number of input latitudes
C> @param IDRTO Output grid identifier
C>  - IDRTO=4 for Gaussian grid
C>  - IDRTO=0 for equally-spaced grid including poles
C>  - IDRTO=256 for equally-spaced grid excluding poles
C> @param IMAXO Even number of output longitudes
C> @param JMAXO Number of output latitudes
C> @param KMAX Number of fields to transform
C> @param IPRIME Input longitude index for the prime meridian.
C>  - Defaults to 1 if IPRIME=0
C>  - Output longitude index for prime meridian assumed 1
C> @param ISKIPI Skip number between input longitudes (defaults to 1 if ISKIPI=0)
C> @param JSKIPI Skip number between input latitudes from south (defaults to -IMAXI if JSKIPI=0)
C> @param KSKIPI Skip number between input grid fields (defaults to IMAXI*JMAXI if KSKIPI=0)
C> @param ISKIPO Skip number between output longitudes (defaults to 1 if ISKIPO=0)
C> @param JSKIPO Skip number between output latitudes from south (defaults to -IMAXO if JSKIPO=0)
C> @param KSKIPO Skip number between output grid fields (defaults to IMAXO*JMAXO if KSKIPO=0)
C> @param JCPU Number of CPUs over which to multiprocess (defaults to environment NCPUS if JCPU=0)
C> @param GRIDI Input grid fields
C> @param GRIDO Output grid fields (may overlay input fields if grid shape is appropriate)
C>
C> @author IREDELL @date 96-02-29
      SUBROUTINE SPTRUN(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,IDRTO,IMAXO,JMAXO,
     &                  KMAX,IPRIME,ISKIPI,JSKIPI,KSKIPI,
     &                  ISKIPO,JSKIPO,KSKIPO,JCPU,GRIDI,GRIDO)
      REAL GRIDI(*),GRIDO(*)
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
      JN=-JSKIPO
      IF(JN.EQ.0) JN=IMAXO
      JS=-JN
      INP=(JMAXO-1)*MAX(0,-JN)+1
      ISP=(JMAXO-1)*MAX(0,-JS)+1
      CALL SPTRAN(IROMB,MAXWV,IDRTO,IMAXO,JMAXO,KMAX,
     &            0,ISKIPO,JN,JS,MDIM,KSKIPO,0,0,JC,
     &            W,GRIDO(INP),GRIDO(ISP),1)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
