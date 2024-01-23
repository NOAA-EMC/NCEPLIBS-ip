C> @file
C> @brief Analyze spectral from Fourier.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 91-10-31 | Mark Iredell | Initial.
C> 94-08-01 | Mark Iredell | Moved zonal wavenumber loop inside.
C> 1998-12-15 | Iredell | Openmp directives inserted.
C>
C> @author Iredell @date 91-10-31

C> Analyzes spectral coefficients from Fourier coefficients
C> for a latitude pair (Northern and Southern hemispheres).
C>
C> Vector components are multiplied by cosine of latitude.
C>
C> @param I spectral domain shape (0 for triangular, 1 for rhomboidal)
C> @param M  spectral truncation
C> @param IM even number of Fourier coefficients
C> @param IX dimension of Fourier coefficients (IX>=IM+2)
C> @param NC dimension of spectral coefficients (NC>=(M+1)*((I+1)*M+2))
C> @param NCTOP dimension of spectral coefficients over top (NCTOP>=2*(M+1))
C> @param KM number of fields
C> @param WGT Gaussian weight
C> @param CLAT cosine of latitude
C> @param PLN Legendre polynomials
C> @param PLNTOP Legendre polynomial over top
C> @param MP identifiers (0 for scalar, 1 for vector)
C> @param F Fourier coefficients combined
C> @param SPC spectral coefficients
C> @param SPCTOP spectral coefficients over top
C>
C> @author Iredell @date 91-10-31
      SUBROUTINE SPANALY(I,M,IM,IX,NC,NCTOP,KM,WGT,CLAT,PLN,PLNTOP,MP,
     &                   F,SPC,SPCTOP)
      INTEGER MP(KM)
      REAL PLN((M+1)*((I+1)*M+2)/2),PLNTOP(M+1)
      REAL F(IX,2,KM)
      REAL SPC(NC,KM),SPCTOP(NCTOP,KM)
      REAL FW(2,2)

C FOR EACH ZONAL WAVENUMBER, ANALYZE TERMS OVER TOTAL WAVENUMBER.
C ANALYZE EVEN AND ODD POLYNOMIALS SEPARATELY.
      LX=MIN(M,IM/2)
!C$OMP PARALLEL DO PRIVATE(L,NT,KS,KP,FW)
      DO K=1,KM
        DO L=0,LX
          NT=MOD(M+1+(I-1)*L,2)+1
          KS=L*(2*M+(I-1)*(L-1))
          KP=KS/2+1
          IF(MP(K).EQ.0) THEN
            FW(1,1)=WGT*(F(2*L+1,1,K)+F(2*L+1,2,K))
            FW(2,1)=WGT*(F(2*L+2,1,K)+F(2*L+2,2,K))
            FW(1,2)=WGT*(F(2*L+1,1,K)-F(2*L+1,2,K))
            FW(2,2)=WGT*(F(2*L+2,1,K)-F(2*L+2,2,K))
          ELSE
            FW(1,1)=WGT*CLAT*(F(2*L+1,1,K)+F(2*L+1,2,K))
            FW(2,1)=WGT*CLAT*(F(2*L+2,1,K)+F(2*L+2,2,K))
            FW(1,2)=WGT*CLAT*(F(2*L+1,1,K)-F(2*L+1,2,K))
            FW(2,2)=WGT*CLAT*(F(2*L+2,1,K)-F(2*L+2,2,K))
            SPCTOP(2*L+1,K)=SPCTOP(2*L+1,K)+PLNTOP(L+1)*FW(1,NT)
            SPCTOP(2*L+2,K)=SPCTOP(2*L+2,K)+PLNTOP(L+1)*FW(2,NT)
          ENDIF
          DO N=L,I*L+M,2
            SPC(KS+2*N+1,K)=SPC(KS+2*N+1,K)+PLN(KP+N)*FW(1,1)
            SPC(KS+2*N+2,K)=SPC(KS+2*N+2,K)+PLN(KP+N)*FW(2,1)
          ENDDO
          DO N=L+1,I*L+M,2
            SPC(KS+2*N+1,K)=SPC(KS+2*N+1,K)+PLN(KP+N)*FW(1,2)
            SPC(KS+2*N+2,K)=SPC(KS+2*N+2,K)+PLN(KP+N)*FW(2,2)
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
