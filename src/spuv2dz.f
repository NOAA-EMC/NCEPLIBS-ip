C> @file
C> @brief Compute divergence and vorticity from winds.
C> @author Iredell @date 92-10-31

C> Computes the divergence and vorticity from wind components
C> in spectral space.
C>
C> Subprogram speps() should be called already.
C>
C> If L is the zonal wavenumber, N is the total wavenumber,
C> EPS(L,N)=SQRT((N**2-L**2)/(4*N**2-1)) and A is earth radius,
C> then the divergence D is computed as:
C> <pre>
C> D(L,N)=I*L*A*U(L,N)
C> +EPS(L,N+1)*N*A*V(L,N+1)-EPS(L,N)*(N+1)*A*V(L,N-1)
C> </pre>
C>
C> and the vorticity Z is computed as:
C> <pre>
C> Z(L,N)=I*L*A*V(L,N)
C> -EPS(L,N+1)*N*A*U(L,N+1)+EPS(L,N)*(N+1)*A*U(L,N-1)
C> </pre>
C>
C> where U is the zonal wind and V is the meridional wind.
C>
C> U and V are weighted by the secant of latitude.
C>
C> Extra terms are used over top of the spectral domain.
C>
C> Advantage is taken of the fact that EPS(L,L)=0
C> in order to vectorize over the entire spectral domain.
C>
C> @param I integer spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param M INTEGER spectral truncation
C> @param ENN1 ((M+1)*((I+1)*M+2)/2) N*(N+1)/A**2
C> @param ELONN1 ((M+1)*((I+1)*M+2)/2) L/(N*(N+1))*A
C> @param EON ((M+1)*((I+1)*M+2)/2) EPSILON/N*A
C> @param EONTOP (M+1) EPSILON/N*A over top
C> @param U ((M+1)*((I+1)*M+2)) zonal wind (over coslat)
C> @param V ((M+1)*((I+1)*M+2)) merid wind (over coslat)
C> @param UTOP (2*(M+1)) zonal wind (over coslat) over top
C> @param VTOP (2*(M+1)) merid wind (over coslat) over top
C> @param D ((M+1)*((I+1)*M+2)) divergence
C> @param Z ((M+1)*((I+1)*M+2)) vorticity
C>
C> @author Iredell @date 92-10-31
      SUBROUTINE SPUV2DZ(I,M,ENN1,ELONN1,EON,EONTOP,U,V,UTOP,VTOP,D,Z)
      REAL ENN1((M+1)*((I+1)*M+2)/2),ELONN1((M+1)*((I+1)*M+2)/2)
      REAL EON((M+1)*((I+1)*M+2)/2),EONTOP(M+1)
      REAL U((M+1)*((I+1)*M+2)),V((M+1)*((I+1)*M+2))
      REAL UTOP(2*(M+1)),VTOP(2*(M+1))
      REAL D((M+1)*((I+1)*M+2)),Z((M+1)*((I+1)*M+2))

C  COMPUTE TERMS FROM THE SPECTRAL DOMAIN
      K=1
      D(2*K-1)=0.
      D(2*K)=0.
      Z(2*K-1)=0.
      Z(2*K)=0.
      DO K=2,(M+1)*((I+1)*M+2)/2-1
        D(2*K-1)=-ELONN1(K)*U(2*K)+EON(K+1)*V(2*K+1)-EON(K)*V(2*K-3)
        D(2*K)=ELONN1(K)*U(2*K-1)+EON(K+1)*V(2*K+2)-EON(K)*V(2*K-2)
        Z(2*K-1)=-ELONN1(K)*V(2*K)-EON(K+1)*U(2*K+1)+EON(K)*U(2*K-3)
        Z(2*K)=ELONN1(K)*V(2*K-1)-EON(K+1)*U(2*K+2)+EON(K)*U(2*K-2)
      ENDDO
      K=(M+1)*((I+1)*M+2)/2
      D(2*K-1)=-ELONN1(K)*U(2*K)-EON(K)*V(2*K-3)
      D(2*K)=ELONN1(K)*U(2*K-1)-EON(K)*V(2*K-2)
      Z(2*K-1)=-ELONN1(K)*V(2*K)+EON(K)*U(2*K-3)
      Z(2*K)=ELONN1(K)*V(2*K-1)+EON(K)*U(2*K-2)

C  COMPUTE TERMS FROM OVER TOP OF THE SPECTRAL DOMAIN
CDIR$ IVDEP
      DO L=0,M
        K=L*(2*M+(I-1)*(L-1))/2+I*L+M+1
        D(2*K-1)=D(2*K-1)+EONTOP(L+1)*VTOP(2*L+1)
        D(2*K)=D(2*K)+EONTOP(L+1)*VTOP(2*L+2)
        Z(2*K-1)=Z(2*K-1)-EONTOP(L+1)*UTOP(2*L+1)
        Z(2*K)=Z(2*K)-EONTOP(L+1)*UTOP(2*L+2)
      ENDDO

C  MULTIPLY BY LAPLACIAN TERM
      DO K=2,(M+1)*((I+1)*M+2)/2
        D(2*K-1)=D(2*K-1)*ENN1(K)
        D(2*K)=D(2*K)*ENN1(K)
        Z(2*K-1)=Z(2*K-1)*ENN1(K)
        Z(2*K)=Z(2*K)*ENN1(K)
      ENDDO
      RETURN
      END
