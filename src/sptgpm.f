C> @file
C> @brief Transform spectral scalar to Mercator
C> ### Program history log:
C> Date       | Programmer | Comments
C> -----------|------------|---------
C>   96-02-29 | IREDELL | Initial.
C> 1998-12-15 | IREDELL | OpenMP directives inserted.
C> @author IREDELL @date 96-02-29

C> This subprogram performs a spherical transform
C> from spectral coefficients of scalar quantities
C> to scalar fields on a Mercator grid.
C> The wave-space can be either triangular or rhomboidal.
C> The wave and grid fields may have general indexing,
C> but each wave field is in sequential 'ibm order',
C> i.e. with zonal wavenumber as the slower index.
C> The Mercator grid is identified by the location
C> of its first point and by its respective increments.
C> The transforms are all multiprocessed over sector points.
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
C> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) if KWSKIP=0)
C> @param KGSKIP Skip number between grid fields
C> (defaults to MI*MJ if KGSKIP=0)
C> @param NISKIP Skip number between grid i-points
C> (defaults to 1 if NISKIP=0)
C> @param NJSKIP Skip number between grid j-points
C> (defaults to MI if NJSKIP=0)
C> @param RLAT1 Latitude of the first grid point in degrees
C> @param RLON1 Longitude of the first grid point in degrees
C> @param DLAT Latitude increment in degrees such that
C> D(PHI)/D(J)=DLAT*COS(PHI) where J is meridional index.
C> DLAT is negative for grids indexed southward.
C> (in terms of grid increment DY valid at latitude RLATI,
C> the latitude increment DLAT is determined as
C> DLAT=DPR*DY/(RERTH*COS(RLATI/DPR))
C> where DPR=180/PI and RERTH is earth's radius)
C> @param DLON Longitude increment in degrees such that
C> D(LAMBDA)/D(I)=DLON where I is zonal index.
C> DLON is negative for grids indexed westward.
C> @param WAVE Wave fields
C> @param GM Mercator fields
C>
C> @author IREDELL @date 96-02-29
      SUBROUTINE SPTGPM(IROMB,MAXWV,KMAX,MI,MJ,
     &                  KWSKIP,KGSKIP,NISKIP,NJSKIP,
     &                  RLAT1,RLON1,DLAT,DLON,WAVE,GM)

      REAL WAVE(*),GM(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(KMAX)
      REAL WTOP(2*(MAXWV+1),KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+3,2,KMAX)
      REAL CLAT(MJ),SLAT(MJ),CLON(MAXWV,MI),SLON(MAXWV,MI)
      PARAMETER(RERTH=6.3712E6)
      PARAMETER(PI=3.14159265358979,DPR=180./PI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MXTOP=MAXWV+1
      IDIM=2*MAXWV+3
      KW=KWSKIP
      KG=KGSKIP
      NI=NISKIP
      NJ=NJSKIP
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=MI*MJ
      IF(NI.EQ.0) NI=1
      IF(NJ.EQ.0) NJ=MI
      DO I=1,MI
        RLON=MOD(RLON1+DLON*(I-1)+3600,360.)
        DO L=1,MAXWV
          CLON(L,I)=COS(L*RLON/DPR)
          SLON(L,I)=SIN(L*RLON/DPR)
        ENDDO
      ENDDO
      YE=1-LOG(TAN((RLAT1+90)/2/DPR))*DPR/DLAT
      DO J=1,MJ
        RLAT=ATAN(EXP(DLAT/DPR*(J-YE)))*2*DPR-90
        CLAT(J)=COS(RLAT/DPR)
        SLAT(J)=SIN(RLAT/DPR)
      ENDDO
      MP=0
C$OMP PARALLEL DO
      DO K=1,KMAX
        WTOP(1:2*MXTOP,K)=0
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM TO GRID
C$OMP PARALLEL DO PRIVATE(PLN,PLNTOP,F,IJK)
      DO J=1,MJ
        CALL SPLEGEND(IROMB,MAXWV,SLAT(J),CLAT(J),EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,KW,2*MXTOP,KMAX,
     &               CLAT(J),PLN,PLNTOP,MP,WAVE,WTOP,F)
        DO K=1,KMAX
          DO I=1,MI
            IJK=(I-1)*NI+(J-1)*NJ+(K-1)*KG+1
            GM(IJK)=F(1,1,K)
          ENDDO
          DO L=1,MAXWV
            DO I=1,MI
              IJK=(I-1)*NI+(J-1)*NJ+(K-1)*KG+1
              GM(IJK)=GM(IJK)+2.*(F(2*L+1,1,K)*CLON(L,I)
     &                           -F(2*L+2,1,K)*SLON(L,I))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
