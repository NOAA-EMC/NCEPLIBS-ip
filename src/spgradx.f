C> @file
C> @brief Compute x-gradient in Fourier space
C> @author IREDELL @date 96-02-20

C> This subprogram computes the x-gradient of fields
C> in complex Fourier space.
C>
C> The x-gradient of a vector field W is
C> WX=CONJG(W)*L/RERTH
C> where L is the wavenumber and RERTH is the Earth radius,
C> so that the result is the x-gradient of the pseudo-vector.
C>
C> The x-gradient of a scalar field W is
C> WX=CONJG(W)*L/(RERTH*CLAT)
C> where CLAT is the cosine of latitude.
C>
C> At the pole this is undefined, so the way to get
C> the x-gradient at the pole is by passing both
C> the weighted wavenumber 0 and the unweighted wavenumber 1 
C> amplitudes at the pole and setting MP=10.
C> In this case, the wavenumber 1 amplitudes are used
C> to compute the x-gradient and then zeroed out.
C>
C> @note This subprogram is thread-safe.
C>
C> @param M Fourier wavenumber truncation
C> @param INCW first dimension of the complex amplitude array
C> (INCW >= M+1)
C> @param KMAX number of Fourier fields
C> @param MP identifiers
C> (0 or 10 for scalar, 1 for vector)
C> @param CLAT cosine of latitude
C> @param[out] W Fourier amplitudes corrected when MP=10 and CLAT=0
C> @param[out] WX complex amplitudes of x-gradients
C>
C> @author IREDELL @date 96-02-20
      SUBROUTINE SPGRADX(M,INCW,KMAX,MP,CLAT,W,WX)

        IMPLICIT NONE
        INTEGER,INTENT(IN):: M,INCW,KMAX,MP(KMAX)
        REAL,INTENT(IN):: CLAT
        REAL,INTENT(INOUT):: W(2*INCW,KMAX)
        REAL,INTENT(OUT):: WX(2*INCW,KMAX)
        INTEGER K,L
        REAL,PARAMETER:: RERTH=6.3712E6

        DO K=1,KMAX
          IF(MP(K).EQ.1) THEN
            DO L=0,M
              WX(2*L+1,K)=-W(2*L+2,K)*(L/RERTH)
              WX(2*L+2,K)=+W(2*L+1,K)*(L/RERTH)
            ENDDO
          ELSEIF(CLAT.EQ.0.) THEN
            DO L=0,M
              WX(2*L+1,K)=0
              WX(2*L+2,K)=0
            ENDDO
            IF(MP(K).EQ.10.AND.M.GE.2) THEN
              WX(3,K)=-W(4,K)/RERTH
              WX(4,K)=+W(3,K)/RERTH
              W(3,K)=0
              W(4,K)=0
            ENDIF
          ELSE
            DO L=0,M
              WX(2*L+1,K)=-W(2*L+2,K)*(L/(RERTH*CLAT))
              WX(2*L+2,K)=+W(2*L+1,K)*(L/(RERTH*CLAT))
            ENDDO
          ENDIF
        ENDDO

      END SUBROUTINE
