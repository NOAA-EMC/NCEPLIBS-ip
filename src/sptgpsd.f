C> @file
C> @brief Transform spectral to polar stereographic gradients
C> ### Program history log:
C> Date | Programmer | Comments
C> -----|------------|----------
C>  96-02-29 | IREDELL | Initial.
C> 1998-12-15 | IREDELL | OpenMP directives inserted.
C> @author IREDELL @date 96-02-29

C> This subprogram performs a spherical transform
C> from spectral coefficients of scalar fields
C> to gradient fields on a pair of polar stereographic grids.
C> The wave-space can be either triangular or rhomboidal.
C> The wave and grid fields may have general indexing,
C> but each wave field is in sequential 'ibm order',
C> i.e., with zonal wavenumber as the slower index.
C> The two square polar stereographic grids are centered
C> on the respective poles, with the orientation longitude
C> of the southern hemisphere grid 180 degrees opposite
C> that of the northern hemisphere grid.
C> The vectors are automatically rotated to be resolved
C> relative to the respective polar stereographic grids.
C>
C> The transform is made efficient by combining points in eight
C> sectors of each polar stereographic grid, numbered as in the
C> following diagram. The pole and the sector boundaries are
C> treated specially in the code. Unfortunately, this approach
C> induces some hairy indexing and code loquacity, for which
C> the developer apologizes.
C>
C> \verbatim
C>   \ 4 | 5 /
C>    \  |  /
C>   3 \ | / 6
C>      \|/
C>   ----+----
C>      /|\
C>   2 / | \ 7
C>    /  |  \
C>   / 1 | 8 \
C> \endverbatim
C>
C> The transforms are all multiprocessed over sector points.
C> transform several fields at a time to improve vectorization.
C> Subprogram can be called from a multiprocessing environment.
C>
C> @param IROMB Spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV Spectral truncation
C> @param KMAX Number of fields to transform
C> @param NPS Odd order of the polar stereographic grids
C> @param KWSKIP Skip number between wave fields
C> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) if KWSKIP=0)
C> @param KGSKIP Skip number between grid fields
C> (defaults to NPS*NPS if KGSKIP=0)
C> @param NISKIP Skip number between grid i-points
C> (defaults to 1 if NISKIP=0)
C> @param NJSKIP Skip number between grid j-points
C> (defaults to NPS if NJSKIP=0)
C> @param TRUE Latitude at which PS grid is true (usually 60.)
C> @param XMESH Grid length at true latitude (M)
C> @param ORIENT Longitude at bottom of northern PS grid
C> (southern PS grid will have opposite orientation.)
C> @param WAVE Wave fields
C> @param XN Northern polar stereographic x-gradients
C> @param YN Northern polar stereographic y-gradients
C> @param XS Southern polar stereographic x-gradients
C> @param YS Southern polar stereographic y-gradients
C>
C> @author IREDELL @date 96-02-29
      SUBROUTINE SPTGPSD(IROMB,MAXWV,KMAX,NPS,
     &                   KWSKIP,KGSKIP,NISKIP,NJSKIP,
     &                   TRUE,XMESH,ORIENT,WAVE,XN,YN,XS,YS)

      REAL WAVE(*),XN(*),YN(*),XS(*),YS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL WD((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
      REAL WZ((MAXWV+1)*((IROMB+1)*MAXWV+2)/2*2+1,KMAX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MDIM=2*MX+1
      KW=KWSKIP
      IF(KW.EQ.0) KW=2*MX
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  CALCULATE GRADIENTS
C$OMP PARALLEL DO PRIVATE(KWS)
      DO K=1,KMAX
        KWS=(K-1)*KW
        CALL SPLAPLAC(IROMB,MAXWV,ENN1,WAVE(KWS+1),WD(1,K),1)
        WZ(1:2*MX,K)=0.
      ENDDO
      CALL SPTGPSV(IROMB,MAXWV,KMAX,NPS,MDIM,KGSKIP,NISKIP,NJSKIP,
     &             TRUE,XMESH,ORIENT,WD,WZ,XN,YN,XS,YS)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
