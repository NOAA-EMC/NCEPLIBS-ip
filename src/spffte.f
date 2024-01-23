C> @file
C> @brief Perform multiple fast Fourier transforms.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 1998-12-18 | Iredell | Initial.
C> 2012-11-12 | Mirvis | fixing hard-wired types problem on Intel/Linux.
C>
C> @author Iredell @date 96-02-20

C> This subprogram performs multiple fast Fourier transforms
C> between complex amplitudes in Fourier space and real values
C> in cyclic physical space.
C>
C> This subprogram must be invoked first with IDIR=0
C> to initialize trigonemetric data. Use subprogram spfft1()
C> to perform an FFT without previous initialization.
C>
C> This version invokes the IBM ESSL FFT.
C>
C> @note The restrictions on IMAX are that it must be a multiple
C> of 1 to 25 factors of two, up to 2 factors of three,
C> and up to 1 factor of five, seven and eleven.
C>
C> If IDIR=0, then W and G need not contain any valid data.
C> The other parameters must be supplied and cannot change
C> in succeeding calls until the next time it is called with IDIR=0.
C>
C> This subprogram is thread-safe.
C>      
C> @param IMAX number of values in the cyclic physical space
C> (see limitations on imax in remarks below.)
C> @param INCW first dimension of the complex amplitude array
C> (INCW >= IMAX/2+1)
C> @param INCG first dimension of the real value array
C> (INCG >= IMAX)
C> @param KMAX number of transforms to perform
C> @param[out] W complex amplitudes if IDIR>0
C> @param[out] G real values if IDIR<0
C> @param IDIR  direction flag
C> - IDIR=0 to initialize trigonometric data
C> - IDIR>0 to transform from Fourier to physical space
C> - IDIR<0 to transform from physical to Fourier space
C> @param[out] AFFT auxiliary array if IDIR<>0
C>
C> @author Iredell @date 96-02-20
      SUBROUTINE SPFFTE(IMAX,INCW,INCG,KMAX,W,G,IDIR,AFFT)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: IMAX,INCW,INCG,KMAX,IDIR
        REAL,INTENT(INOUT):: W(2*INCW,KMAX)
        REAL,INTENT(INOUT):: G(INCG,KMAX)
        REAL(8),INTENT(INOUT):: AFFT(50000+4*IMAX)
        INTEGER:: INIT,INC2X,INC2Y,N,M,ISIGN,NAUX1,NAUX2,NAUX3
C ==EM==       ^(4)
        REAL:: SCALE
        REAL(8):: AUX2(20000+2*IMAX),AUX3
        INTEGER:: IACR,IARC

        NAUX1=25000+2*IMAX
        NAUX2=20000+2*IMAX
        NAUX3=1
        IACR=1
        IARC=1+NAUX1

C  INITIALIZATION.
C  FILL AUXILIARY ARRAYS WITH TRIGONOMETRIC DATA
        SELECT CASE(IDIR)
          CASE(0)
            INIT=1
            INC2X=INCW
            INC2Y=INCG
            N=IMAX
            M=KMAX
            ISIGN=-1
            SCALE=1.
            IF(DIGITS(1.).LT.DIGITS(1._8)) THEN
              CALL SCRFT(INIT,W,INC2X,G,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IACR),NAUX1,AUX2,NAUX2,AUX3,NAUX3)
            ELSE
              CALL DCRFT(INIT,W,INC2X,G,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IACR),NAUX1,AUX2,NAUX2)
            ENDIF
            INIT=1
            INC2X=INCG
            INC2Y=INCW
            N=IMAX
            M=KMAX
            ISIGN=+1
            SCALE=1./IMAX
            IF(DIGITS(1.).LT.DIGITS(1._8)) THEN
              CALL SRCFT(INIT,G,INC2X,W,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IARC),NAUX1,AUX2,NAUX2,AUX3,NAUX3)
            ELSE
              CALL DRCFT(INIT,G,INC2X,W,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IARC),NAUX1,AUX2,NAUX2)
            ENDIF

C  FOURIER TO PHYSICAL TRANSFORM.
          CASE(1:)
            INIT=0
            INC2X=INCW
            INC2Y=INCG
            N=IMAX
            M=KMAX
            ISIGN=-1
            SCALE=1.
            IF(DIGITS(1.).LT.DIGITS(1._8)) THEN
              CALL SCRFT(INIT,W,INC2X,G,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IACR),NAUX1,AUX2,NAUX2,AUX3,NAUX3)
            ELSE
              CALL DCRFT(INIT,W,INC2X,G,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IACR),NAUX1,AUX2,NAUX2)
            ENDIF

C  PHYSICAL TO FOURIER TRANSFORM.
          CASE(:-1)
            INIT=0
            INC2X=INCG
            INC2Y=INCW
            N=IMAX
            M=KMAX
            ISIGN=+1
            SCALE=1./IMAX
            IF(DIGITS(1.).LT.DIGITS(1._8)) THEN
              CALL SRCFT(INIT,G,INC2X,W,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IARC),NAUX1,AUX2,NAUX2,AUX3,NAUX3)
            ELSE
              CALL DRCFT(INIT,G,INC2X,W,INC2Y,N,M,ISIGN,SCALE,
     &                   AFFT(IARC),NAUX1,AUX2,NAUX2)
            ENDIF
        END SELECT
      END SUBROUTINE
