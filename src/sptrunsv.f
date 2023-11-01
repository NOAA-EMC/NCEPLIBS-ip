C> @file
C> @brief Spectrally interpolate vectors to polar stereo.
C>
C> 96-02-29 | Iredell | Initial.
C> 1998-12-15 | Iredell | Openmp directives inserted.
C>
C> @author Iredell @date 96-02-29

C> This subprogram spectrally truncates vector fields
C> on a global cylindrical grid, returning the fields
C> to specific pairs of polar stereographic scalar fields.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a gaussian grid.
C>
C> The grid fields may have general indexing.
C>
C> The transforms are all multiprocessed.
C>
C> Transform several fields at a time to improve vectorization.
C>
C> Subprogram can be called from a multiprocessing environment.
C>
C> Minimum grid dimensions for unaliased transforms to spectral:
C> Dimension                    |Linear              |Quadratic
C> -----------------------      |---------           |-------------
C> IMAX                         |2*MAXWV+2           |3*MAXWV/2*2+2
C> JMAX (IDRT=4,IROMB=0)        |1*MAXWV+1           |3*MAXWV/2+1
C> JMAX (IDRT=4,IROMB=1)        |2*MAXWV+1           |5*MAXWV/2+1
C> JMAX (IDRT=0,IROMB=0)        |2*MAXWV+3           |3*MAXWV/2*2+3
C> JMAX (IDRT=0,IROMB=1)        |4*MAXWV+3           |5*MAXWV/2*2+3
C> JMAX (IDRT=256,IROMB=0)      |2*MAXWV+1           |3*MAXWV/2*2+1
C> JMAX (IDRT=256,IROMB=1)      |4*MAXWV+1           |5*MAXWV/2*2+1
C>      
C> @param IROMB integer spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV integer spectral truncation
C> @param IDRTI integer input grid identifier
C> - IDRTI=4 for Gaussian grid
C> - IDRTI=0 for equally-spaced grid including poles
C> - IDRTI=256 for equally-spaced grid excluding poles
C> @param IMAXI integer even number of input longitudes.
C> @param JMAXI integer number of input latitudes.
C> @param KMAX integer number of fields to transform.
C> @param NPS integer odd order of the polar stereographic grids
C> @param IPRIME integer input longitude index for the prime meridian.
C> (defaults to 1 if IPRIME=0)
C> (output longitude index for prime meridian assumed 1.)
C> @param ISKIPI integer skip number between input longitudes
C> (defaults to 1 if ISKIPI=0)
C> @param JSKIPI integer skip number between input latitudes from south
C> (defaults to -IMAXI if JSKIPI=0)
C> @param KSKIPI integer skip number between input grid fields
C> (defaults to IMAXI*JMAXI if KSKIPI=0)
C> @param KGSKIP integer skip number between grid fields
C> (defaults to NPS*NPS if KGSKIP=0)
C> @param NISKIP integer skip number between grid i-points
C> (defaults to 1 if NISKIP=0)
C> @param NJSKIP integer skip number between grid j-points
C> (defaults to NPS if NJSKIP=0)
C> @param JCPU integer number of cpus over which to multiprocess
C> (defaults to environment NCPUS if JCPU=0)
C> @param TRUE real latitude at which ps grid is true (usually 60.)
C> @param XMESH real grid length at true latitude (m)
C> @param ORIENT real longitude at bottom of Northern PS grid
C> (Southern PS grid will have opposite orientation.)
C> @param GRIDUI real input grid u-winds
C> @param GRIDVI real input grid v-winds
C> @param LUV logical flag whether to return winds
C> @param LDZ logical flag whether to return divergence and vorticity
C> @param LPS logical flag whether to return potential and streamfcn
C> @param UN real northern ps u-winds if luv
C> @param VN real northern ps v-winds if luv
C> @param US real southern ps u-winds if luv
C> @param VS real southern ps v-winds if luv
C> @param DN real northern divergences if ldz
C> @param ZN real northern vorticities if ldz
C> @param DS real southern divergences if ldz
C> @param ZS real southern vorticities if ldz
C> @param PN real northern potentials if lps
C> @param SN real northern streamfcns if lps
C> @param PS real southern potentials if lps
C> @param SS real southern streamfcns if lps
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTRUNSV(IROMB,MAXWV,IDRTI,IMAXI,JMAXI,KMAX,NPS,
     &                    IPRIME,ISKIPI,JSKIPI,KSKIPI,KGSKIP,
     &                    NISKIP,NJSKIP,JCPU,TRUE,XMESH,ORIENT,
     &                    GRIDUI,GRIDVI,
     &                    LUV,UN,VN,US,VS,LDZ,DN,ZN,DS,ZS,
     &                    LPS,PN,SN,PS,SS)
      LOGICAL LUV,LDZ,LPS
      REAL GRIDUI(*),GRIDVI(*)
      REAL UN(*),VN(*),US(*),VS(*),DN(*),ZN(*),DS(*),ZS(*)
      REAL PN(*),SN(*),PS(*),SS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)

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

C  TRANSFORM WAVE TO OUTPUT WINDS
      IF(LUV) THEN
        CALL SPTGPSV(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &               TRUE,XMESH,ORIENT,WD,WZ,UN,VN,US,VS)
      ENDIF

C  TRANSFORM WAVE TO OUTPUT DIVERGENCE AND VORTICITY
      IF(LDZ) THEN
        CALL SPTGPS(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &              TRUE,XMESH,ORIENT,WD,DN,DS)
        CALL SPTGPS(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &              TRUE,XMESH,ORIENT,WZ,ZN,ZS)
      ENDIF

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
        CALL SPTGPS(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &              TRUE,XMESH,ORIENT,WD,PN,PS)
        CALL SPTGPS(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &              TRUE,XMESH,ORIENT,WZ,SN,SS)
      ENDIF
      END
