C> @file
C> @brief Synthesize Fourier coefficients from spectral coefficients.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 91-10-31 | Mark Iredell | Initial.
C> 1998-12-18 | Mark Iredell | Include scalar and gradient option.
C>
C> @author Iredell @date 92-10-31

C> Synthesizes Fourier coefficients from spectral coefficients
C> for a latitude pair (Northern and Southern hemispheres).
C>
C> Vector components are divided by cosine of latitude.
C>
C> @param I spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param IM even number of Fourier coefficients
C> @param IX dimension of Fourier coefficients (IX>=IM+2)
C> @param NC dimension of spectral coefficients
C> (NC>=(M+1)*((I+1)*M+2))
C> @param NCTOP dimension of spectral coefficients over top
C> (NCTOP>=2*(M+1))
C> @param KM number of fields
C> @param CLAT cosine of latitude
C> @param PLN ((M+1)*((I+1)*M+2)/2) Legendre polynomial
C> @param PLNTOP Legendre polynomial over top
C> @param SPC spectral coefficients
C> @param SPCTOP spectral coefficients over top
C> @param MP identifiers (0 for scalar, 1 for vector,
C> or 10 for scalar and gradient)
C> @param F Fourier coefficients for latitude pair
C>
C> @author Iredell @date 92-10-31
      SUBROUTINE SPSYNTH(I,M,IM,IX,NC,NCTOP,KM,CLAT,PLN,PLNTOP,MP,
     &                   SPC,SPCTOP,F)

      REAL PLN((M+1)*((I+1)*M+2)/2),PLNTOP(M+1)
      INTEGER MP(KM)
      REAL SPC(NC,KM),SPCTOP(NCTOP,KM)
      REAL F(IX,2,KM)

C  ZERO OUT FOURIER COEFFICIENTS.
      DO K=1,KM
        DO L=0,IM/2
          F(2*L+1,1,K)=0.
          F(2*L+2,1,K)=0.
          F(2*L+1,2,K)=0.
          F(2*L+2,2,K)=0.
        ENDDO
      ENDDO

C  SYNTHESIS OVER POLE.
C  INITIALIZE FOURIER COEFFICIENTS WITH TERMS OVER TOP OF THE SPECTRUM.
C  INITIALIZE EVEN AND ODD POLYNOMIALS SEPARATELY.
      IF(CLAT.EQ.0) THEN
        LTOPE=MOD(M+1+I,2)
!C$OMP PARALLEL DO PRIVATE(LB,LE,L,KS,KP,N,F1R,F1I)
        DO K=1,KM
          LB=MP(K)
          LE=MP(K)
          IF(MP(K).EQ.10) THEN
            LB=0
            LE=1
          ENDIF
          L=LB
          IF(L.EQ.1) THEN
            IF(L.EQ.LTOPE) THEN
              F(2*L+1,1,K)=PLNTOP(L+1)*SPCTOP(2*L+1,K)
              F(2*L+2,1,K)=PLNTOP(L+1)*SPCTOP(2*L+2,K)
            ELSE
              F(2*L+1,2,K)=PLNTOP(L+1)*SPCTOP(2*L+1,K)
              F(2*L+2,2,K)=PLNTOP(L+1)*SPCTOP(2*L+2,K)
            ENDIF
          ENDIF
C  FOR EACH ZONAL WAVENUMBER, SYNTHESIZE TERMS OVER TOTAL WAVENUMBER.
C  SYNTHESIZE EVEN AND ODD POLYNOMIALS SEPARATELY.
          DO L=LB,LE
            KS=L*(2*M+(I-1)*(L-1))
            KP=KS/2+1
            DO N=L,I*L+M,2
              F(2*L+1,1,K)=F(2*L+1,1,K)+PLN(KP+N)*SPC(KS+2*N+1,K)
              F(2*L+2,1,K)=F(2*L+2,1,K)+PLN(KP+N)*SPC(KS+2*N+2,K)
            ENDDO
            DO N=L+1,I*L+M,2
              F(2*L+1,2,K)=F(2*L+1,2,K)+PLN(KP+N)*SPC(KS+2*N+1,K)
              F(2*L+2,2,K)=F(2*L+2,2,K)+PLN(KP+N)*SPC(KS+2*N+2,K)
            ENDDO
C  SEPARATE FOURIER COEFFICIENTS FROM EACH HEMISPHERE.
C  ODD POLYNOMIALS CONTRIBUTE NEGATIVELY TO THE SOUTHERN HEMISPHERE.
            F1R=F(2*L+1,1,K)
            F1I=F(2*L+2,1,K)
            F(2*L+1,1,K)=F1R+F(2*L+1,2,K)
            F(2*L+2,1,K)=F1I+F(2*L+2,2,K)
            F(2*L+1,2,K)=F1R-F(2*L+1,2,K)
            F(2*L+2,2,K)=F1I-F(2*L+2,2,K)
          ENDDO
        ENDDO

C  SYNTHESIS OVER FINITE LATITUDE.
C  INITIALIZE FOURIER COEFFICIENTS WITH TERMS OVER TOP OF THE SPECTRUM.
C  INITIALIZE EVEN AND ODD POLYNOMIALS SEPARATELY.
      ELSE
        LX=MIN(M,IM/2)
        LTOPE=MOD(M+1,2)
        LTOPO=1-LTOPE
        LE=1+I*LTOPE
        LO=2-I*LTOPO
!C$OMP PARALLEL DO PRIVATE(L,KS,KP,N,F1R,F1I)
        DO K=1,KM
          IF(MP(K).EQ.1) THEN
            DO L=LTOPE,LX,2
              F(2*L+1,LE,K)=PLNTOP(L+1)*SPCTOP(2*L+1,K)
              F(2*L+2,LE,K)=PLNTOP(L+1)*SPCTOP(2*L+2,K)
            ENDDO
            DO L=LTOPO,LX,2
              F(2*L+1,LO,K)=PLNTOP(L+1)*SPCTOP(2*L+1,K)
              F(2*L+2,LO,K)=PLNTOP(L+1)*SPCTOP(2*L+2,K)
            ENDDO
          ENDIF
C  FOR EACH ZONAL WAVENUMBER, SYNTHESIZE TERMS OVER TOTAL WAVENUMBER.
C  SYNTHESIZE EVEN AND ODD POLYNOMIALS SEPARATELY.
          DO L=0,LX
            KS=L*(2*M+(I-1)*(L-1))
            KP=KS/2+1
            DO N=L,I*L+M,2
              F(2*L+1,1,K)=F(2*L+1,1,K)+PLN(KP+N)*SPC(KS+2*N+1,K)
              F(2*L+2,1,K)=F(2*L+2,1,K)+PLN(KP+N)*SPC(KS+2*N+2,K)
            ENDDO
            DO N=L+1,I*L+M,2
              F(2*L+1,2,K)=F(2*L+1,2,K)+PLN(KP+N)*SPC(KS+2*N+1,K)
              F(2*L+2,2,K)=F(2*L+2,2,K)+PLN(KP+N)*SPC(KS+2*N+2,K)
            ENDDO
          ENDDO
C  SEPARATE FOURIER COEFFICIENTS FROM EACH HEMISPHERE.
C  ODD POLYNOMIALS CONTRIBUTE NEGATIVELY TO THE SOUTHERN HEMISPHERE.
C  DIVIDE VECTOR COMPONENTS BY COSINE LATITUDE.
          DO L=0,LX
            F1R=F(2*L+1,1,K)
            F1I=F(2*L+2,1,K)
            F(2*L+1,1,K)=F1R+F(2*L+1,2,K)
            F(2*L+2,1,K)=F1I+F(2*L+2,2,K)
            F(2*L+1,2,K)=F1R-F(2*L+1,2,K)
            F(2*L+2,2,K)=F1I-F(2*L+2,2,K)
          ENDDO
          IF(MP(K).EQ.1) THEN
            DO L=0,LX
              F(2*L+1,1,K)=F(2*L+1,1,K)/CLAT
              F(2*L+2,1,K)=F(2*L+2,1,K)/CLAT
              F(2*L+1,2,K)=F(2*L+1,2,K)/CLAT
              F(2*L+2,2,K)=F(2*L+2,2,K)/CLAT
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      END
