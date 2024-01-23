C> @file
C> @brief Compute laplacian in spectral space.
C> @author Iredell @date 92-10-31

C> Computes the laplacian or the inverse laplacian
C> of a scalar field in spectral space.
C>      
C> Subprogram speps() should be called already.
C>      
C> The Laplacian of Q(L,N) is simply -N*(N+1)/A**2*Q(L,N)
C>
C> @param I spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param M spectral truncation
C> @param ENN1 N*(N+1)/A**2
C> @param[out] Q if IDIR > 0, scalar field
C> (Q(0,0) is not computed)
C> @param[out] QD2 if IDIR < 0, Laplacian
C> @param IDIR flag
C> - IDIR > 0 to take Laplacian
C> - IDIR < 0 to take inverse Laplacian
C>
C> @author Iredell @date 92-10-31
      SUBROUTINE SPLAPLAC(I,M,ENN1,Q,QD2,IDIR)

      REAL ENN1((M+1)*((I+1)*M+2)/2)
      REAL Q((M+1)*((I+1)*M+2))
      REAL QD2((M+1)*((I+1)*M+2))

C  TAKE LAPLACIAN
      IF(IDIR.GT.0) THEN
        K=1
        QD2(2*K-1)=0.
        QD2(2*K)=0.
        DO K=2,(M+1)*((I+1)*M+2)/2
          QD2(2*K-1)=Q(2*K-1)*(-ENN1(K))
          QD2(2*K)=Q(2*K)*(-ENN1(K))
        ENDDO

C  TAKE INVERSE LAPLACIAN
      ELSE
        DO K=2,(M+1)*((I+1)*M+2)/2
          Q(2*K-1)=QD2(2*K-1)/(-ENN1(K))
          Q(2*K)=QD2(2*K)/(-ENN1(K))
        ENDDO
      ENDIF

      RETURN
      END
