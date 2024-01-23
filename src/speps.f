C> @file
C> @brief Compute utility spectral fields.
C> @author Iredell  @date 92-10-31

C> Computes constant fields indexed in the spectral domain
C> in "IBM ORDER" (Zonal wavenumber is the slower index).
C>      
C> If L is the zonal wavenumber and N is the total wavenumber
C> and A is the earth radius, then the fields returned are:
C> - (1) normalizing factor EPSILON=SQRT((N**2-L**2)/(4*N**2-1))
C> - (2) Laplacian factor N*(N+1)/A**2
C> - (3) zonal derivative/Laplacian factor L/(N*(N+1))*A
C> - (4) Meridional derivative/Laplacian factor EPSILON/N*A
C>
C> @param I spectral domain shape (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param EPS ((M+1)*((I+1)*M+2)/2) SQRT((N**2-L**2)/(4*N**2-1))
C> @param EPSTOP (M+1) SQRT((N**2-L**2)/(4*N**2-1)) OVER TOP
C> @param ENN1 ((M+1)*((I+1)*M+2)/2) N*(N+1)/A**2
C> @param ELONN1 ((M+1)*((I+1)*M+2)/2) L/(N*(N+1))*A
C> @param EON ((M+1)*((I+1)*M+2)/2) EPSILON/N*A
C> @param EONTOP (M+1) EPSILON/N*A OVER TOP
C>
C> @author Iredell  @date 92-10-31
      SUBROUTINE SPEPS(I,M,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      REAL EPS((M+1)*((I+1)*M+2)/2),EPSTOP(M+1)
      REAL ENN1((M+1)*((I+1)*M+2)/2),ELONN1((M+1)*((I+1)*M+2)/2)
      REAL EON((M+1)*((I+1)*M+2)/2),EONTOP(M+1)
      PARAMETER(RERTH=6.3712E6,RA2=1./RERTH**2)

      DO L=0,M
        K=L*(2*M+(I-1)*(L-1))/2+L+1
        EPS(K)=0.
        ENN1(K)=RA2*L*(L+1)
        ELONN1(K)=RERTH/(L+1)
        EON(K)=0.
      ENDDO
      DO L=0,M
        DO N=L+1,I*L+M
          K=L*(2*M+(I-1)*(L-1))/2+N+1
          EPS(K)=SQRT(FLOAT(N**2-L**2)/FLOAT(4*N**2-1))
          ENN1(K)=RA2*N*(N+1)
          ELONN1(K)=RERTH*L/(N*(N+1))
          EON(K)=RERTH/N*EPS(K)
        ENDDO
      ENDDO
      DO L=0,M
        N=I*L+M+1
        EPSTOP(L+1)=SQRT(FLOAT(N**2-L**2)/FLOAT(4*N**2-1))
        EONTOP(L+1)=RERTH/N*EPSTOP(L+1)
      ENDDO
      RETURN
      END
