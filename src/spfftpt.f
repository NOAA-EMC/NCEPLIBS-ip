C> @file
C> @brief Compute fourier transform to gridpoints.
C> @author Iredell @date 96-02-20

C> This subprogram computes a slow Fourier transform
C> from Fourier space to a set of gridpoints.
C>
C> @note This subprogram is thread-safe.
C>
C> @param M Fourier wavenumber truncation
C> @param N number of gridpoints
C> @param INCW first dimension of the complex amplitude array
C> (INCW >= M+1)
C> @param INCG first dimension of the gridpoint array
C> (INCG >= N)
C> @param KMAX number of Fourier fields
C> @param RLON grid longitudes in degrees
C> @param W Fourier amplitudes
C> @param G gridpoint values
C>
C> @author Iredell @date 96-02-20
      SUBROUTINE SPFFTPT(M,N,INCW,INCG,KMAX,RLON,W,G)

        IMPLICIT NONE
        INTEGER,INTENT(IN):: M,N,INCW,INCG,KMAX
        REAL,INTENT(IN):: RLON(N)
        REAL,INTENT(IN):: W(2*INCW,KMAX)
        REAL,INTENT(OUT):: G(INCG,KMAX)
        INTEGER I,K,L
        REAL RADLON,SLON(M),CLON(M)
        REAL,PARAMETER:: PI=3.14159265358979

        DO I=1,N
          RADLON=PI/180*RLON(I)
          DO L=1,M
            SLON(L)=SIN(L*RADLON)
            CLON(L)=COS(L*RADLON)
          ENDDO
          DO K=1,KMAX
            G(I,K)=W(1,K)
          ENDDO
          DO L=1,M
            DO K=1,KMAX
              G(I,K)=G(I,K)+2.*(W(2*L+1,K)*CLON(L)-W(2*L+2,K)*SLON(L))
            ENDDO
          ENDDO
        ENDDO
      END SUBROUTINE
