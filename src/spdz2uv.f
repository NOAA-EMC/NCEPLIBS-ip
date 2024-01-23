C> @file
C> @brief Compute winds from divergence and vorticity.
C> @author Iredell @date 92-10-31

C> Computes the wind components from divergence and vorticity
C> in spectral space.
C>
C> Subprogram speps() should be called already.
C>
C> If L is the zonal wavenumber, N is the total wavenumber,
C> <pre>      
C> EPS(L,N) = SQRT((N**2-L**2)/(4*N**2-1))
C> </pre>
C> and A is earth radius,
C> then the zonal wind component U is computed as
C> <pre>
C> U(L,N)=-I*L/(N*(N+1))*A*D(L,N)
C> +EPS(L,N+1)/(N+1)*A*Z(L,N+1)-EPS(L,N)/N*A*Z(L,N-1)
C> </pre>
C> and the meridional wind component V is computed as
C> <pre>
C> V(L,N)=-I*L/(N*(N+1))*A*Z(L,N)
C> -EPS(L,N+1)/(N+1)*A*D(L,N+1)+EPS(L,N)/N*A*D(L,N-1)
C> </pre>
C> where D is divergence and Z is vorticity.
C>
C> U and V are weighted by the cosine of latitude.
C>
C> Cxtra terms are computed over top of the spectral domain.
C>
C> Advantage is taken of the fact that EPS(L,L)=0
C> in order to vectorize over the entire spectral domain.
C>
C> @param I spectral domain shape (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param ENN1 ((M+1)*((I+1)*M+2)/2) N*(N+1)/A**2
C> @param ELONN1 ((M+1)*((I+1)*M+2)/2) L/(N*(N+1))*A
C> @param EON ((M+1)*((I+1)*M+2)/2) EPSILON/N*A
C> @param EONTOP (M+1) EPSILON/N*A OVER TOP
C> @param D ((M+1)*((I+1)*M+2)) divergence
C> @param Z ((M+1)*((I+1)*M+2)) vorticity
C> @param U ((M+1)*((I+1)*M+2)) zonal wind (times coslat)
C> @param V ((M+1)*((I+1)*M+2)) merid wind (times coslat)
C> @param UTOP (2*(M+1)) zonal wind (times coslat) over top
C> @param VTOP (2*(M+1)) merid wind (times coslat) over top
C>
C> @author Iredell @date 92-10-31
      SUBROUTINE SPDZ2UV(I,M,ENN1,ELONN1,EON,EONTOP,D,Z,U,V,UTOP,VTOP)
      REAL ENN1((M+1)*((I+1)*M+2)/2),ELONN1((M+1)*((I+1)*M+2)/2)
      REAL EON((M+1)*((I+1)*M+2)/2),EONTOP(M+1)
      REAL D((M+1)*((I+1)*M+2)),Z((M+1)*((I+1)*M+2))
      REAL U((M+1)*((I+1)*M+2)),V((M+1)*((I+1)*M+2))
      REAL UTOP(2*(M+1)),VTOP(2*(M+1))

C  COMPUTE WINDS IN THE SPECTRAL DOMAIN
      K=1
      U(2*K-1)=EON(K+1)*Z(2*K+1)
      U(2*K)=EON(K+1)*Z(2*K+2)
      V(2*K-1)=-EON(K+1)*D(2*K+1)
      V(2*K)=-EON(K+1)*D(2*K+2)
      DO K=2,(M+1)*((I+1)*M+2)/2-1
        U(2*K-1)=ELONN1(K)*D(2*K)+EON(K+1)*Z(2*K+1)-EON(K)*Z(2*K-3)
        U(2*K)=-ELONN1(K)*D(2*K-1)+EON(K+1)*Z(2*K+2)-EON(K)*Z(2*K-2)
        V(2*K-1)=ELONN1(K)*Z(2*K)-EON(K+1)*D(2*K+1)+EON(K)*D(2*K-3)
        V(2*K)=-ELONN1(K)*Z(2*K-1)-EON(K+1)*D(2*K+2)+EON(K)*D(2*K-2)
      ENDDO
      K=(M+1)*((I+1)*M+2)/2
      U(2*K-1)=ELONN1(K)*D(2*K)-EON(K)*Z(2*K-3)
      U(2*K)=-ELONN1(K)*D(2*K-1)-EON(K)*Z(2*K-2)
      V(2*K-1)=ELONN1(K)*Z(2*K)+EON(K)*D(2*K-3)
      V(2*K)=-ELONN1(K)*Z(2*K-1)+EON(K)*D(2*K-2)

C  COMPUTE WINDS OVER TOP OF THE SPECTRAL DOMAIN
      DO L=0,M
        K=L*(2*M+(I-1)*(L-1))/2+I*L+M+1
        UTOP(2*L+1)=-EONTOP(L+1)*Z(2*K-1)
        UTOP(2*L+2)=-EONTOP(L+1)*Z(2*K)
        VTOP(2*L+1)=EONTOP(L+1)*D(2*K-1)
        VTOP(2*L+2)=EONTOP(L+1)*D(2*K)
      ENDDO
      RETURN
      END
