C> @file
C> @brief Compute variance by total wavenumber.
C> @author Iredell @date 92-10-31

C> Computes the variances by total wavenumber
C> of a scalar field in spectral space.
C>
C> @param I spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param Q ((M+1)*((I+1)*M+2)) scalar field
C> @param QVAR (0:(I+1)*M) variances
C>
C> @author Iredell @date 92-10-31
      SUBROUTINE SPVAR(I,M,Q,QVAR)
      REAL Q((M+1)*((I+1)*M+2))
      REAL QVAR(0:(I+1)*M)

      L=0
      DO N=0,M
        KS=L*(2*M+(I-1)*(L-1))+2*N
        QVAR(N)=0.5*Q(KS+1)**2
      ENDDO
      DO N=M+1,(I+1)*M
        QVAR(N)=0.
      ENDDO
      DO N=0,(I+1)*M
        DO L=MAX(1,N-M),MIN(N,M)
          KS=L*(2*M+(I-1)*(L-1))+2*N
          QVAR(N)=QVAR(N)+Q(KS+1)**2+Q(KS+2)**2
        ENDDO
      ENDDO

      RETURN
      END
