C> @file
C> @brief Perform a scalar spherical transform
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | Iredell | Initial.
C> 1998-12-15 | Iredell | Generic fft used, openmp directives inserted
C> 2013-01-16 | Iredell, Mirvis | Fixing afft negative sharing effect
C>
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform between spectral
C> coefficients of scalar quantities and fields on a global
C> cylindrical grid.
C>
C> The wave-space can be either triangular or
C> rhomboidal.  The grid-space can be either an equally-spaced grid
C> (with or without pole points) or a Gaussian grid.
C>
C> The wave and grid fields may have general indexing,
C> but each wave field is in sequential 'ibm order',
C> i.e. with zonal wavenumber as the slower index.
C>
C> Transforms are done in latitude pairs for efficiency;
C> thus grid arrays for each hemisphere must be passed.
C>
C> If so requested, just a subset of the latitude pairs
C> may be transformed in each invocation of the subprogram.
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
C> IMAX                         | 2*MAXWV+2          | 3*MAXWV/2*2+2
C> JMAX (IDRT=4,IROMB=0)        | 1*MAXWV+1          | 3*MAXWV/2+1
C> JMAX (IDRT=4,IROMB=1)        | 2*MAXWV+1          | 5*MAXWV/2+1
C> JMAX (IDRT=0,IROMB=0)        | 2*MAXWV+3          | 3*MAXWV/2*2+3
C> JMAX (IDRT=0,IROMB=1)        | 4*MAXWV+3          | 5*MAXWV/2*2+3
C> JMAX (IDRT=256,IROMB=0)      | 2*MAXWV+1          | 3*MAXWV/2*2+1
C> JMAX (IDRT=256,IROMB=1)      | 4*MAXWV+1          | 5*MAXWV/2*2+1
C>
C> @param IROMB spectral domain shape
c> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param IDRT grid identifier
C> - IDRT=4 for Gaussian grid,
C> - IDRT=0 for equally-spaced grid including poles
C> - IDRT=256 for equally-spaced grid excluding poles
C> @param IMAX even number of longitudes.
C> @param JMAX number of latitudes.
C> @param KMAX number of fields to transform.
C> @param IP longitude index for the prime meridian
C> @param IS skip number between longitudes
C> @param JN skip number between n.h. latitudes from north
C> @param JS skip number between s.h. latitudes from south
C> @param KW skip number between wave fields
C> @param KG skip number between grid fields
C> @param JB latitude index (from pole) to begin transform
C> @param JE latitude index (from pole) to end transform
C> @param JC number of cpus over which to multiprocess
C> @param[out] WAVE wave fields if IDIR>0
C> @param[out] GRIDN n.h. grid fields (starting at JB) if IDIR<0
C> @param[out] GRIDS s.h. grid fields (starting at JB) if IDIR<0
C> @param IDIR transform flag
C> (IDIR>0 for wave to grid, IDIR<0 for grid to wave)
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTRANF(IROMB,MAXWV,IDRT,IMAX,JMAX,KMAX,
     &                   IP,IS,JN,JS,KW,KG,JB,JE,JC,
     &                   WAVE,GRIDN,GRIDS,IDIR)

      REAL WAVE(*),GRIDN(*),GRIDS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL(8) AFFT(50000+4*IMAX), AFFT_TMP(50000+4*IMAX)
      REAL CLAT(JB:JE),SLAT(JB:JE),WLAT(JB:JE)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2,JB:JE)
      REAL PLNTOP(MAXWV+1,JB:JE)
      REAL WTOP(2*(MAXWV+1))
      REAL G(IMAX,2)
!      write(0,*) 'sptranf top'

C  SET PARAMETERS
      MP=0
      CALL SPTRANF0(IROMB,MAXWV,IDRT,IMAX,JMAX,JB,JE,
     &              EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP,
     &              AFFT,CLAT,SLAT,WLAT,PLN,PLNTOP)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  TRANSFORM WAVE TO GRID
      IF(IDIR.GT.0) THEN
C$OMP PARALLEL DO PRIVATE(AFFT_TMP,KWS,WTOP,G,IJKN,IJKS)
        DO K=1,KMAX
			  AFFT_TMP=AFFT
          KWS=(K-1)*KW
          WTOP=0
          DO J=JB,JE
            CALL SPTRANF1(IROMB,MAXWV,IDRT,IMAX,JMAX,J,J,
     &                    EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP,
     &                    AFFT_TMP,CLAT(J),SLAT(J),WLAT(J),
     &                    PLN(1,J),PLNTOP(1,J),MP,
     &                    WAVE(KWS+1),WTOP,G,IDIR)
            IF(IP.EQ.1.AND.IS.EQ.1) THEN
              DO I=1,IMAX
                IJKN=I+(J-JB)*JN+(K-1)*KG
                IJKS=I+(J-JB)*JS+(K-1)*KG
                GRIDN(IJKN)=G(I,1)
                GRIDS(IJKS)=G(I,2)
              ENDDO
            ELSE
              DO I=1,IMAX
                IJKN=MOD(I+IP-2,IMAX)*IS+(J-JB)*JN+(K-1)*KG+1
                IJKS=MOD(I+IP-2,IMAX)*IS+(J-JB)*JS+(K-1)*KG+1
                GRIDN(IJKN)=G(I,1)
                GRIDS(IJKS)=G(I,2)
              ENDDO
            ENDIF
          ENDDO
        ENDDO

C  TRANSFORM GRID TO WAVE
      ELSE
C$OMP PARALLEL DO PRIVATE(AFFT_TMP,KWS,WTOP,G,IJKN,IJKS)
        DO K=1,KMAX
			  AFFT_TMP=AFFT
          KWS=(K-1)*KW
          WTOP=0
          DO J=JB,JE
            IF(WLAT(J).GT.0.) THEN
              IF(IP.EQ.1.AND.IS.EQ.1) THEN
                DO I=1,IMAX
                  IJKN=I+(J-JB)*JN+(K-1)*KG
                  IJKS=I+(J-JB)*JS+(K-1)*KG
                  G(I,1)=GRIDN(IJKN)
                  G(I,2)=GRIDS(IJKS)
                ENDDO
              ELSE
                DO I=1,IMAX
                  IJKN=MOD(I+IP-2,IMAX)*IS+(J-JB)*JN+(K-1)*KG+1
                  IJKS=MOD(I+IP-2,IMAX)*IS+(J-JB)*JS+(K-1)*KG+1
                  G(I,1)=GRIDN(IJKN)
                  G(I,2)=GRIDS(IJKS)
                ENDDO
              ENDIF
              CALL SPTRANF1(IROMB,MAXWV,IDRT,IMAX,JMAX,J,J,
     &                      EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP,
     &                      AFFT_TMP,CLAT(J),SLAT(J),WLAT(J),
     &                      PLN(1,J),PLNTOP(1,J),MP,
     &                      WAVE(KWS+1),WTOP,G,IDIR)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      END
