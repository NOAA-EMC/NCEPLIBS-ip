C> @file
C> @brief Pad or truncate a spectral field.
C> @author Iredell @date 92-10-31

C> Pad or truncate a spectral field.
C>
C> @param I1 input spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param M1 input spectral truncation
C> @param Q1 ((M+1)*((I+1)*M+2)) input field
C> @param I2 output spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param M2 output spectral truncation
C> @param Q2 ((M+1)*((I+1)*M+2)) output field
C>
C> @author Iredell @date 92-10-31
      SUBROUTINE SPPAD(I1,M1,Q1,I2,M2,Q2)

      REAL Q1((M1+1)*((I1+1)*M1+2))
      REAL Q2((M2+1)*((I2+1)*M2+2))

      DO L=0,M2
        DO N=L,I2*L+M2
          KS2=L*(2*M2+(I2-1)*(L-1))+2*N
          IF(L.LE.M1.AND.N.LE.I1*L+M1) THEN
            KS1=L*(2*M1+(I1-1)*(L-1))+2*N
            Q2(KS2+1)=Q1(KS1+1)
            Q2(KS2+2)=Q1(KS1+2)
          ELSE
            Q2(KS2+1)=0
            Q2(KS2+2)=0
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
