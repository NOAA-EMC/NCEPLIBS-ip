C> @file
C>
C> Compute Legendre polynomials
C> @author IREDELL @date 92-10-31

C> Evaluates the orthonormal associated Legendre polynomials in the
C> spectral domain at a given latitude. Subprogram splegend should
C> be called already. If l is the zonal wavenumber, N is the total
C> wavenumber, and EPS(L,N)=SQRT((N**2-L**2)/(4*N**2-1)) then the
C> following bootstrapping formulas are used:
C>
C> <pre>
C> PLN(0,0)=SQRT(0.5)
C> PLN(L,L)=PLN(L-1,L-1)*CLAT*SQRT(FLOAT(2*L+1)/FLOAT(2*L))
C> PLN(L,N)=(SLAT*PLN(L,N-1)-EPS(L,N-1)*PLN(L,N-2))/EPS(L,N)
C> </pre>
C>
C> Synthesis at the pole needs only two zonal wavenumbers. Scalar
C> fields are synthesized with zonal wavenumber 0 while vector
C> fields are synthesized with zonal wavenumber 1. (Thus polar
C> vector fields are implicitly divided by clat.) The following
C> bootstrapping formulas are used at the pole:
C>
C> <pre>
C> PLN(0,0)=SQRT(0.5)
C> PLN(1,1)=SQRT(0.75)
C> PLN(L,N)=(PLN(L,N-1)-EPS(L,N-1)*PLN(L,N-2))/EPS(L,N)
C> </pre>
C>
C> PROGRAM HISTORY LOG:
C> -  91-10-31  MARK IREDELL
C> - 98-06-10  MARK IREDELL  GENERALIZE PRECISION
C>
C> @param I        - INTEGER SPECTRAL DOMAIN SHAPE
C>                (0 FOR TRIANGULAR, 1 FOR RHOMBOIDAL)
C> @param M        - INTEGER SPECTRAL TRUNCATION
C> @param SLAT     - REAL SINE OF LATITUDE
C> @param CLAT     - REAL COSINE OF LATITUDE
C> @param EPS      - REAL ((M+1)*((I+1)*M+2)/2) SQRT((N**2-L**2)/(4*N**2-1))
C> @param EPSTOP   - REAL (M+1) SQRT((N**2-L**2)/(4*N**2-1)) OVER TOP
C> @param[out] PLN - REAL ((M+1)*((I+1)*M+2)/2) LEGENDRE POLYNOMIAL
C> @param[out] PLNTOP - REAL (M+1) LEGENDRE POLYNOMIAL OVER TOP
C>
      SUBROUTINE SPLEGEND(I,M,SLAT,CLAT,EPS,EPSTOP,PLN,PLNTOP)

CFPP$ NOCONCUR R
      REAL EPS((M+1)*((I+1)*M+2)/2),EPSTOP(M+1)
      REAL PLN((M+1)*((I+1)*M+2)/2),PLNTOP(M+1)
      REAL(KIND=SELECTED_REAL_KIND(15,45)):: DLN((M+1)*((I+1)*M+2)/2)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ITERATIVELY COMPUTE PLN WITHIN SPECTRAL DOMAIN AT POLE
      M1=M+1
      M2=2*M+I+1
      MX=(M+1)*((I+1)*M+2)/2
      IF(CLAT.EQ.0.) THEN
        DLN(1)=SQRT(0.5)
        IF(M.GT.0) THEN
          DLN(M1+1)=SQRT(0.75)
          DLN(2)=SLAT*DLN(1)/EPS(2)
        ENDIF
        IF(M.GT.1) THEN
          DLN(M1+2)=SLAT*DLN(M1+1)/EPS(M1+2)
          DLN(3)=(SLAT*DLN(2)-EPS(2)*DLN(1))/EPS(3)
          DO N=3,M
            K=1+N
            DLN(K)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPS(K)
            K=M1+N
            DLN(K)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPS(K)
          ENDDO
          IF(I.EQ.1) THEN
            K=M2
            DLN(K)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPS(K)
          ENDIF
          DO K=M2+1,MX
            DLN(K)=0.
          ENDDO
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE POLYNOMIALS OVER TOP OF SPECTRAL DOMAIN
        K=M1+1
        PLNTOP(1)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPSTOP(1)
        IF(M.GT.0) THEN
          K=M2+1
          PLNTOP(2)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPSTOP(2)
          DO L=2,M
            PLNTOP(L+1)=0.
          ENDDO
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  ITERATIVELY COMPUTE PLN(L,L) (BOTTOM HYPOTENUSE OF DOMAIN)
      ELSE
        NML=0
        K=1
        DLN(K)=SQRT(0.5)
        DO L=1,M+(I-1)*NML
          KP=K
          K=L*(2*M+(I-1)*(L-1))/2+L+NML+1
          DLN(K)=DLN(KP)*CLAT*SQRT(FLOAT(2*L+1)/FLOAT(2*L))
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE PLN(L,L+1) (DIAGONAL NEXT TO BOTTOM HYPOTENUSE OF DOMAIN)
        NML=1
CDIR$ IVDEP
        DO L=0,M+(I-1)*NML
          K=L*(2*M+(I-1)*(L-1))/2+L+NML+1
          DLN(K)=SLAT*DLN(K-1)/EPS(K)
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE REMAINING PLN IN SPECTRAL DOMAIN
        DO NML=2,M
CDIR$ IVDEP
          DO L=0,M+(I-1)*NML
            K=L*(2*M+(I-1)*(L-1))/2+L+NML+1
            DLN(K)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPS(K)
          ENDDO
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE POLYNOMIALS OVER TOP OF SPECTRAL DOMAIN
        DO L=0,M
          NML=M+1+(I-1)*L
          K=L*(2*M+(I-1)*(L-1))/2+L+NML+1
          PLNTOP(L+1)=(SLAT*DLN(K-1)-EPS(K-1)*DLN(K-2))/EPSTOP(L+1)
        ENDDO
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  RETURN VALUES
      DO K=1,MX
        PLN(K)=DLN(K)
      ENDDO
      RETURN
      END
