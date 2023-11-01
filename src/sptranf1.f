C> @file
C> @brief Sptranf spectral transform.
C> @author Iredell @date 96-02-29

C> This subprogram performs an single latitude transform for
C> subprogram sptranf(). Use this subprogram outside
C> the sptranf() family context at your own risk.
C>
C> @param IROMB spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param IDRT grid identifier
C> - IDRT=4 for Gaussian grid,
C> - IDRT=0 for equally-spaced grid including poles,
C> - IDRT=256 for equally-spaced grid excluding poles
C> @param IMAX even number of longitudes
C> @param JMAX number of latitudes
C> @param JB latitude index (from pole) to begin transform
C> @param JE latitude index (from pole) to end transform
C> @param EPS
C> @param EPSTOP
C> @param ENN1
C> @param ELONN1
C> @param EON
C> @param EONTOP
C> @param CLAT cosines of latitude
C> @param SLAT sines of latitude
C> @param WLAT Gaussian weights
C> @param AFFT auxiliary array if IDIR=0
C> @param PLN Legendre polynomials
C> @param PLNTOP Legendre polynomial over top
C> @param MP identifier (0 for scalar, 1 for vector)
C> @param[out] W wave field if IDIR>0
C> @param[out] WTOP wave field over top if IDIR>0
C> @param[out] G grid field if IDIR<0
C> @param IDIR transform flag
C> (IDIR>0 for wave to grid, IDIR<0 for grid to wave)
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTRANF1(IROMB,MAXWV,IDRT,IMAX,JMAX,JB,JE,
     &                    EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP,
     &                    AFFT,CLAT,SLAT,WLAT,PLN,PLNTOP,MP,
     &                    W,WTOP,G,IDIR)

      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL(8) AFFT(50000+4*IMAX)
      REAL CLAT(JB:JE),SLAT(JB:JE),WLAT(JB:JE)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2,JB:JE)
      REAL PLNTOP(MAXWV+1,JB:JE)
      REAL W((MAXWV+1)*((IROMB+1)*MAXWV+2))
      REAL WTOP(2*(MAXWV+1))
      REAL G(IMAX,2,JB:JE)
      REAL F(IMAX+2,2)

      KW=(MAXWV+1)*((IROMB+1)*MAXWV+2)
      KWTOP=2*(MAXWV+1)
      IF(IDIR.GT.0) THEN
        DO J=JB,JE
          CALL SPSYNTH(IROMB,MAXWV,IMAX,IMAX+2,KW,KWTOP,1,
     &                 CLAT(J),PLN(1,J),PLNTOP(1,J),MP,
     &                 W,WTOP,F)
          CALL SPFFTE(IMAX,(IMAX+2)/2,IMAX,2,F,G(1,1,J),+1,AFFT)
        ENDDO
      ELSE
        DO J=JB,JE
          CALL SPFFTE(IMAX,(IMAX+2)/2,IMAX,2,F,G(1,1,J),-1,AFFT)
          CALL SPANALY(IROMB,MAXWV,IMAX,IMAX+2,KW,KWTOP,1,
     &                 WLAT(J),CLAT(J),PLN(1,J),PLNTOP(1,J),MP,
     &                 F,W,WTOP)
        ENDDO
      ENDIF

      END
