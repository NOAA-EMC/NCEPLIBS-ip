C> @file
C> @brief Transform spectral to Mercator gradients.
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
C> from spectral coefficients of scalar fields
C> to gradient fields on a Mercator grid.
C>
C> The wave-space can be either triangular or rhomboidal.
C> The wave and grid fields may have general indexing,
C> but each wave field is in sequential 'IBM order',
C> i.e. with zonal wavenumber as the slower index.
C>
C> The Mercator grid is identified by the location
C> of its first point and by its respective increments.
C>
C> The transforms are all multiprocessed over sector points.
C>
C> Transform several fields at a time to improve vectorization.
C> Subprogram can be called from a multiprocessing environment.
C>
C> @param IROMB Spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV Spectral truncation
C> @param KMAX Number of fields to transform
C> @param MI Number of points in the faster zonal direction
C> @param MJ Number of points in the slower merid direction
C> @param KWSKIP Skip number between wave fields
C>                (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) if KWSKIP=0)
C> @param KGSKIP Skip number between grid fields
C>                (defaults to MI*MJ if KGSKIP=0)
C> @param NISKIP Skip number between grid i-points
C>                (defaults to 1 if NISKIP=0)
C> @param NJSKIP Skip number between grid j-points
C>                (defaults to MI if NJSKIP=0)
C> @param RLAT1 Latitude of the first grid point in degrees
C> @param RLON1 Longitude of the first grid point in degrees
C> @param DLAT Latitude increment in degrees such that
C> D(PHI)/D(J)=DLAT*COS(PHI) where J is meridional index.
C> DLAT is negative for grids indexed southward.
C> (in terms of grid increment dy valid at latitude RLATI,
C> the latitude increment DLAT is determined as
C> DLAT=DPR*DY/(RERTH*COS(RLATI/DPR))
C> where DPR=180/PI and RERTH is Earth's radius)
C> @param DLON Longitude increment in degrees such that
C> D(LAMBDA)/D(I)=DLON where I is zonal index.
C> DLON is negative for grids indexed westward.
C> @param WAVE Wave fields
C> @param XM Mercator x-gradients
C> @param YM Mercator y-gradients
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTGPMD(IROMB,MAXWV,KMAX,MI,MJ,
     &                   KWSKIP,KGSKIP,NISKIP,NJSKIP,
     &                   RLAT1,RLON1,DLAT,DLON,WAVE,XM,YM)

      REAL WAVE(*),XM(*),YM(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)

C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      KW=KWSKIP
      IF(KW.EQ.0) KW=2*MX

C  CALCULATE GRADIENTS
C$OMP PARALLEL DO PRIVATE(KWS)
      DO K=1,KMAX
        KWS=(K-1)*KW
        CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),1)
        WZ(1:2*MX,K)=0.
      ENDDO
      CALL SPTGPMV(IROMB,MAXWV,KMAX,MI,MJ,MDIM,KGSKIP,NISKIP,NJSKIP,
     &             RLAT1,RLON1,DLAT,DLON,WD,WZ,XM,YM)
      END
