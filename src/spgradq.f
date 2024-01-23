C> @file
C> @brief Compute gradient in spectral space.
C> @author Iredell @date 92-10-31

C> Computes the horizontal vector gradient of a scalar field
C> in spectral space.
C>
C> Subprogram speps() should be called already.
C>
C> If l is the zonal wavenumber, n is the total wavenumber,
C> eps(l,n)=sqrt((n**2-l**2)/(4*n**2-1)) and a is earth radius,
C> then the zonal gradient of q(l,n) is simply i*l/a*q(l,n)
C> while the meridional gradient of q(l,n) is computed as
C> eps(l,n+1)*(n+2)/a*q(l,n+1)-eps(l,n+1)*(n-1)/a*q(l,n-1).
C>
C> Extra terms are computed over top of the spectral domain.
C>
C> Advantage is taken of the fact that eps(l,l)=0
C> in order to vectorize over the entire spectral domain.
C>
C> @param I spectral domain shape (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param ENN1 
C> @param ELONN1 
C> @param EON EPSILON/N*A
C> @param EONTOP EPSILON/N*A over top
C> @param Q scalar field
C> @param QDX zonal gradient (times coslat)
C> @param QDY merid gradient (times coslat)
C> @param QDYTOP merid gradient (times coslat) over top
C>      
C> @author IREDELL @date 92-10-31
      SUBROUTINE SPGRADQ(I,M,ENN1,ELONN1,EON,EONTOP,Q,QDX,QDY,QDYTOP)

      REAL ENN1((M+1)*((I+1)*M+2)/2),ELONN1((M+1)*((I+1)*M+2)/2)
      REAL EON((M+1)*((I+1)*M+2)/2),EONTOP(M+1)
      REAL Q((M+1)*((I+1)*M+2))
      REAL QDX((M+1)*((I+1)*M+2)),QDY((M+1)*((I+1)*M+2))
      REAL QDYTOP(2*(M+1))

C  TAKE ZONAL AND MERIDIONAL GRADIENTS
      K=1
      QDX(2*K-1)=0.
      QDX(2*K)=0.
      QDY(2*K-1)=EON(K+1)*ENN1(K+1)*Q(2*K+1)
      QDY(2*K)=EON(K+1)*ENN1(K+1)*Q(2*K+2)
      DO K=2,(M+1)*((I+1)*M+2)/2-1
        QDX(2*K-1)=-ELONN1(K)*ENN1(K)*Q(2*K)
        QDX(2*K)=ELONN1(K)*ENN1(K)*Q(2*K-1)
        QDY(2*K-1)=EON(K+1)*ENN1(K+1)*Q(2*K+1)-EON(K)*ENN1(K-1)*Q(2*K-3)
        QDY(2*K)=EON(K+1)*ENN1(K+1)*Q(2*K+2)-EON(K)*ENN1(K-1)*Q(2*K-2)
      ENDDO
      K=(M+1)*((I+1)*M+2)/2
      QDX(2*K-1)=-ELONN1(K)*ENN1(K)*Q(2*K)
      QDX(2*K)=ELONN1(K)*ENN1(K)*Q(2*K-1)
      QDY(2*K-1)=-EON(K)*ENN1(K-1)*Q(2*K-3)
      QDY(2*K)=-EON(K)*ENN1(K-1)*Q(2*K-2)

C  TAKE MERIDIONAL GRADIENT OVER TOP
      DO L=0,M
        K=L*(2*M+(I-1)*(L-1))/2+I*L+M+1
        QDYTOP(2*L+1)=-EONTOP(L+1)*ENN1(K)*Q(2*K-1)
        QDYTOP(2*L+2)=-EONTOP(L+1)*ENN1(K)*Q(2*K)
      ENDDO
      RETURN
      END
