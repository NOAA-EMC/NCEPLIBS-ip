C> @file
C> @brief Sptranf spectral initialization.
C> @author IREDELL @date 96-02-29

C> This subprogram performs an initialization for
C> subprogram sptranf(). Use this subprogram outside
C> the sptranf() family context at your own risk.
C>
C> @param IROMB spectral domain shape
c> (0 for triangular, 1 for rhomboidal)
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
C> @param AFFT auxiliary array if IDIR=0
C> @param CLAT cosines of latitude
C> @param SLAT sines of latitude
C> @param WLAT Gaussian weights
C> @param PLN Legendre polynomials
C> @param PLNTOP Legendre polynomial over top
C>
C> @author IREDELL @date 96-02-29
      SUBROUTINE SPTRANF0(IROMB,MAXWV,IDRT,IMAX,JMAX,JB,JE,
     &                    EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP,
     &                    AFFT,CLAT,SLAT,WLAT,PLN,PLNTOP)

      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      REAL(8) AFFT(50000+4*IMAX)
      REAL CLAT(JB:JE),SLAT(JB:JE),WLAT(JB:JE)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2,JB:JE)
      REAL PLNTOP(MAXWV+1,JB:JE)
      REAL SLATX(JMAX),WLATX(JMAX)

      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      CALL SPFFTE(IMAX,(IMAX+2)/2,IMAX,2,0.,0.,0,AFFT)
      CALL SPLAT(IDRT,JMAX,SLATX,WLATX)
      JHE=(JMAX+1)/2
      IF(JHE.GT.JMAX/2) WLATX(JHE)=WLATX(JHE)/2
      DO J=JB,JE
        CLAT(J)=SQRT(1.-SLATX(J)**2)
        SLAT(J)=SLATX(J)
        WLAT(J)=WLATX(J)
      ENDDO
C$OMP PARALLEL DO
      DO J=JB,JE
        CALL SPLEGEND(IROMB,MAXWV,SLAT(J),CLAT(J),EPS,EPSTOP,
     &                PLN(1,J),PLNTOP(1,J))
      ENDDO

      END
