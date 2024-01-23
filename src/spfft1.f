C> @file
C> @brief Perform multiple fast Fourier transforms.
C> @author Iredell @date 96-02-20

C> This subprogram performs multiple fast Fourier transforms
C> between complex amplitudes in Fourier space and real values
C> in cyclic physical space.
C>
C> Subprogram spfft1() initializes trigonometric data each call.
C> Use subprogram spfft() to save time and initialize once.
C> This version invokes the IBM ESSL FFT.
C>
C> @note The restrictions on IMAX are that it must be a multiple of 1
C> to 25 factors of two, up to 2 factors of three, and up to 1 factor of
C> five, seven and eleven.
C>
C> @note This subprogram is thread-safe.
C>
C> @param IMAX number of values in the cyclic physical space
C> (see limitations on imax in remarks below.)
C> @param INCW first dimension of the complex amplitude array
C> (INCW >= IMAX/2+1)
C> @param INCG first dimension of the real value array (INCG >= IMAX)
C> @param KMAX number of transforms to perform
C> @param[out] W complex amplitudes if IDIR>0
C> @param[out] G values if IDIR<0
C> @param IDIR direction flag
C> - IDIR>0 to transform from Fourier to physical space
C> - IDIR<0 to transform from physical to Fourier space
C>
C> @author Iredell @date 96-02-20
      SUBROUTINE SPFFT1(IMAX,INCW,INCG,KMAX,W,G,IDIR)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: IMAX,INCW,INCG,KMAX,IDIR
        COMPLEX,INTENT(INOUT):: W(INCW,KMAX)
        REAL,INTENT(INOUT):: G(INCG,KMAX)
        REAL:: AUX1(25000+INT(0.82*IMAX))
        REAL:: AUX2(20000+INT(0.57*IMAX))
        INTEGER:: NAUX1,NAUX2
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NAUX1=25000+INT(0.82*IMAX)
        NAUX2=20000+INT(0.57*IMAX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FOURIER TO PHYSICAL TRANSFORM.
        SELECT CASE(IDIR)
          CASE(1:)
            CALL SCRFT(1,W,INCW,G,INCG,IMAX,KMAX,-1,1.,
     &                 AUX1,NAUX1,AUX2,NAUX2,0.,0)
            CALL SCRFT(0,W,INCW,G,INCG,IMAX,KMAX,-1,1.,
     &                 AUX1,NAUX1,AUX2,NAUX2,0.,0)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PHYSICAL TO FOURIER TRANSFORM.
          CASE(:-1)
            CALL SRCFT(1,G,INCG,W,INCW,IMAX,KMAX,+1,1./IMAX,
     &               AUX1,NAUX1,AUX2,NAUX2,0.,0)
            CALL SRCFT(0,G,INCG,W,INCW,IMAX,KMAX,+1,1./IMAX,
     &               AUX1,NAUX1,AUX2,NAUX2,0.,0)
        END SELECT
      END SUBROUTINE
