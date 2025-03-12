C> @file
C> @brief A concatenation of several subroutines from the
C> <a target="_blank" href="https://netlib.org/fftpack/">FFTPACK</a> collection.
C> FFTPACK is a package of Fortran subprograms for the fast Fourier
C> transform of periodic and other symmetric sequences. It includes
C> complex, real, sine, cosine, and quarter-wave transforms.
C>
C> See
C> <a target="_blank" href="https://netlib.org/fftpack/doc">FFTPACK documentation</a>,
C> <a target="_blank" href="https://help.graphica.com.au/irix-6.5.30/man/3S/dzfft">IRIX man pages</a>, and
C> <a target="_blank" href="https://www.ibm.com/docs/en/essl/6.2?topic=subroutines-scrft-dcrft-complex-real-fourier-transform">IBM ESSL documentation</a>
C> for further historical context for various subroutines included in this file.
C>      
C> Reference:
C> - P.N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations
C> (G. Rodrigue, ed.), Academic Press, 1982, pp. 51--83.
C> 
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO

C> @brief Computes a set of `m` real discrete `n`-point Fourier transforms of complex conjugate even data
C>
C> The 'd' stands for double precision, but this is historical and may be ignored,
C> and this subroutine is identical to `scrft()`.
C>
C> @param init  Initialization flag; if nonzero, initializes the transformation
C> @param x     Input array of real values, dimensions (`2*ldx`, `m`)
C> @param ldx   Leading dimension of the input array `x`
C> @param y     Output array storing the transformed values, dimensions (`ldy`, `m`)
C> @param ldy   Leading dimension of the output array `y`
C> @param n     Number of data points in each transform (first dimension of `y`)
C> @param m     Number of transforms to be computed (second dimension of `y`)
C> @param isign Sign indicator for the transform direction (not explicitly used)
C> @param scale Scaling factor applied to the transformed values after computation
C> @param table Work array of size 44002, precomputed during initialization
C> @param n1    Auxiliary parameter
C> @param wrk   Work array used for intermediate storage during computation
C> @param n2    Auxiliary parameter
C> @param z     Optional dummy parameter for compatibility
C> @param nz    Optional dummy parameter for compatibility
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE dcrft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)
 
      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,i,j
      real x(2*ldx,*),y(ldy,*),scale,table(44002),wrk(*)
      real, optional :: z
      integer, optional :: nz
 
      IF (init.ne.0) THEN
        CALL rffti(n,table)
      ELSE
!OCL NOVREC
        DO j=1,m
          y(1,j)=x(1,j)
          DO i=2,n
            y(i,j)=x(i+1,j)
          ENDDO
          CALL rfftb(n,y(1,j),table)
          DO i=1,n
            y(i,j)=scale*y(i,j)
          ENDDO
        ENDDO
      ENDIF
 
      RETURN
      END

C> @brief Compute a set of m real discrete n-point Fourier transforms of complex conjugate even data
C>
C> The 's' stands for single precision, but this is historical and may be ignored,
C> and this subroutine is identical to `dcrft()`.
C>
C> @param init  Initialization flag; if nonzero, initializes the transformation
C> @param x     Input array of real values, dimensions (`2*ldx`, `m`)
C> @param ldx   Leading dimension of the input array `x`
C> @param y     Output array storing the transformed values, dimensions (`ldy`, `m`)
C> @param ldy   Leading dimension of the output array `y`
C> @param n     Number of data points in each transform (first dimension of `y`)
C> @param m     Number of transforms to be computed (second dimension of `y`)
C> @param isign Sign indicator for the transform direction (not explicitly used)
C> @param scale Scaling factor applied to the transformed values after computation
C> @param table Work array of size 44002, precomputed during initialization
C> @param n1    Auxiliary parameter
C> @param wrk   Work array used for intermediate storage during computation
C> @param n2    Auxiliary parameter
C> @param z     Optional dummy parameter for compatibility
C> @param nz    Optional dummy parameter for compatibility
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE scrft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)
 
      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,i,j
      real x(2*ldx,*),y(ldy,*),scale,table(44002),wrk(*)
      real, optional :: z
      integer, optional :: nz
 
      IF (init.ne.0) THEN
        CALL rffti(n,table)
      ELSE
!OCL NOVREC
        DO j=1,m
          y(1,j)=x(1,j)
          DO i=2,n
            y(i,j)=x(i+1,j)
          ENDDO
          CALL rfftb(n,y(1,j),table)
          DO i=1,n
            y(i,j)=scale*y(i,j)
          ENDDO
        ENDDO
      ENDIF
 
      RETURN
      END

C> @brief Compute a complex discrete Fourier transform of real data
C>
C> This subroutine performs a Fourier transform on real input data.
C> It supports initialization of the trigonometric coefficient 
C> table and forward or inverse transformations. The `isign` parameter determines 
C> the operation mode.
C>
C> @param isign  Operation mode: 0 for initialization, 1 for inverse transform
C> @param n      Number of data points in the transform
C> @param scale  Scaling factor applied to the transformed values
C> @param x      Input array of real values, length at least `n+1`
C> @param y      Output array storing the transformed values, length at least `n`
C> @param table  Work array used for storing trigonometric coefficients
C> @param work   Work array used for intermediate computations
C> @param isys   Dummy parameter
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE csfft(isign,n,scale,x,y,table,work,isys)
 
      implicit none
      integer isign,n,isys,i
      real scale,x(*),y(*),table(*),work(*)
 
      IF (isign.eq.0) THEN
        CALL rffti(n,table)
      ENDIF
      IF (isign.eq.1) THEN
        y(1)=x(1)
        DO i=2,n
          y(i)=x(i+1)
        ENDDO
        CALL rfftb(n,y,table)
        DO i=1,n
          y(i)=scale*y(i)
        ENDDO
      ENDIF
 
      RETURN
      END

C> @brief Compute a set of m complex discrete n-point Fourier transforms of real data
C>
C> This subroutine performs multiple (`m`) discrete Fourier transforms of length `n`
C> on real input data. The transforms are computed using the real-to-complex Fast Fourier 
C> Transform (FFT) approach. The subroutine requires two invocations: once with `init` set
C> to a non-zero value to initialize work storage, the second with `init=0` to perform the DFT.
C>
C> The 'd' stands for double precision, but this is historical and may be ignored,
C> and this subroutine is identical to `srcft()`.
C>
C> @param init  Initialization flag; if nonzero, initializes the transformation
C> @param x     Input array of real values, dimensions (`ldx`, `m`)
C> @param ldx   Leading dimension of the input array `x`
C> @param y     Output array storing the transformed values, dimensions (`2*ldy`, `m`)
C> @param ldy   Leading dimension of the output array `y`
C> @param n     Number of data points in each transform (first dimension of `y`)
C> @param m     Number of independent transforms to be computed (second dimension of `y`)
C> @param isign Sign indicator for the transform direction (not explicitly used)
C> @param scale Scaling factor applied to the transformed values
C> @param table Work array of size 44002
C> @param n1    Auxiliary parameter
C> @param wrk   Work array used for intermediate storage during computation
C> @param n2    Auxiliary parameter
C> @param z     Optional dummy parameter for compatibility
C> @param nz    Optional dummy parameter for compatibility
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE drcft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)
 
      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,i,j
      real x(ldx,*),y(2*ldy,*),scale,table(44002),wrk(*)
      real, optional :: z
      integer, optional :: nz
 
      IF (init.ne.0) THEN
        CALL rffti(n,table)
      ELSE
        DO j=1,m
          DO i=1,n
            y(i,j)=x(i,j)
          ENDDO
          CALL rfftf(n,y(1,j),table)
          DO i=1,n
            y(i,j)=scale*y(i,j)
          ENDDO
          DO i=n,2,-1
            y(i+1,j)=y(i,j)
          ENDDO
          y(2,j)=0.
C 01/17/2013 vvvvvvvvvvvvv E.Mirvis added ver 2.0.1 by S.Moorthi request. No +|- demo.
          y(n+2,j) = 0.
        ENDDO
      ENDIF
 
      RETURN
      END

C> @brief Compute a set of m complex discrete n-point Fourier transforms of real data
C>
C> This subroutine performs multiple (`m`) discrete Fourier transforms of length `n`
C> on real input data. The transforms are computed using the real-to-complex Fast Fourier 
C> Transform (FFT) approach. The subroutine requires two invocations: once with `init` set
C> to a non-zero value to initialize work storage, the second with `init=0` to perform the DFT.
C>
C> The 's' stands for single precision, but this is historical and may be ignored,
C> and this subroutine is identical to `drcft()`.
C>
C> @param init  Initialization flag; if nonzero, initializes the transformation
C> @param x     Input array of real values, dimensions (`ldx`, `m`)
C> @param ldx   Leading dimension of the input array `x`
C> @param y     Output array storing the transformed values, dimensions (`2*ldy`, `m`)
C> @param ldy   Leading dimension of the output array `y`
C> @param n     Number of data points in each transform (first dimension of `y`)
C> @param m     Number of independent transforms to be computed (second dimension of `y`)
C> @param isign Sign indicator for the transform direction (not explicitly used)
C> @param scale Scaling factor applied to the transformed values
C> @param table Work array of size 44002
C> @param n1    Auxiliary parameter
C> @param wrk   Work array used for intermediate storage during computation
C> @param n2    Auxiliary parameter
C> @param z     Optional dummy parameter for compatibility
C> @param nz    Optional dummy parameter for compatibility
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE srcft(init,x,ldx,y,ldy,n,m,isign,scale,
     &                 table,n1,wrk,n2,z,nz)
 
      implicit none
      integer init,ldx,ldy,n,m,isign,n1,n2,i,j
      real x(ldx,*),y(2*ldy,*),scale,table(44002),wrk(*)
      real, optional :: z
      integer, optional :: nz
 
      IF (init.ne.0) THEN
        CALL rffti(n,table)
      ELSE
        DO j=1,m
          DO i=1,n
            y(i,j)=x(i,j)
          ENDDO
          CALL rfftf(n,y(1,j),table)
          DO i=1,n
            y(i,j)=scale*y(i,j)
          ENDDO
          DO i=n,2,-1
            y(i+1,j)=y(i,j)
          ENDDO
          y(2,j)=0.
          y(n+2,j) = 0.
C 01/17/2013 ^^^^^^^^^^E.Mirvis added ver 2.0.1 by S.Moorthi request. No +|- demo.
        ENDDO
      ENDIF
 
      RETURN
      END

C> @brief Compute a real-to-complex discrete Fourier transform
C>
C> This subroutine performs a Fourier transform on real input data. 
C> It supports initialization of the trigonometric coefficient table and forward transformation.
C>
C> @param isign  Operation mode: 0 for initialization, -1 for forward transform
C> @param n      Number of data points in the transform
C> @param scale  Scaling factor applied to the transformed values
C> @param x      Input array of real values, length at least `n`
C> @param y      Output array storing the transformed values, length at least `n+1`
C> @param table  Work array used for storing trigonometric coefficients
C> @param work   Work array used for intermediate computations
C> @param isys   Dummy parameter
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE scfft(isign,n,scale,x,y,table,work,isys)
 
      implicit none
      integer isign,n,isys,i
      real scale,x(*),y(*),table(*),work(*)
 
      IF (isign.eq.0) THEN
        CALL rffti(n,table)
      ENDIF
      IF (isign.eq.-1) THEN
        DO i=1,n
          y(i)=x(i)
        ENDDO
        CALL rfftf(n,y,table)
        DO i=1,n
          y(i)=scale*y(i)
        ENDDO
        DO i=n,2,-1
          y(i+1)=y(i)
        ENDDO
        y(2)=0.
      ENDIF
 
      RETURN
      END

C> @brief Compute a forward transform of a real periodic sequence
C>
C> The wsave array must be initialized by calling subroutine rffti(n,wsave) and a
C> different wsave array must be used for each different value of n. This
C> initialization does not have to be repeated so long as n remains unchanged thus
C> subsequent transforms can be obtained faster than the first.
C>
C> @note This transform is unnormalized since a call of `rfftf()` followed by a
C> call of `rfftb()` will multiply the input sequence by n.
C>
C> @param n     The length of the array `r` to be transformed
C> @param r     A real array of length `n` which contains the sequence to be transformed
C> @param wsave A work array which must be dimensioned at least `2*n+15`; must not be destroyed between calls of `rfftf()` or `rfftb()`
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RFFTF (N,R,WSAVE)
      DIMENSION       R(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
      
C> @brief Compute the real periodic sequence from its Fourier coefs
C>      
C> @param n     The length of the array `r` to be transformed
C> @param r     A real array of length `n` which contains the sequence to be transformed
C> @param wsave A work array which must be dimensioned at least `2*n+15`; must not be destroyed between calls of `rfftf()` or `rfftb()`
C>
C> @note This transform is unnormalized since a call of `rfftf()` followed by a
C> call of `rfftb()` will multiply the input sequence by n.
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RFFTB (N,R,WSAVE)
      DIMENSION       R(*)       ,WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END

C> @brief Initialize the array wsave which is used in both `rfftf()` and `rfftb()`.
C>      
C> @param n     The length of the sequence to be transformed
C> @param wsave A work array which must be dimensioned at least `2*n+15`; must not be destroyed between calls of `rfftf()` or `rfftb()`
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RFFTI (N,WSAVE)
      DIMENSION       WSAVE(*)
      IF (N .EQ. 1) RETURN
      CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END

C> @brief Compute the inverse fast Fourier transform (IFFT) using a mixed-radix algorithm
C>
C> This subroutine is a low-level component of the real-to-complex IFFT computation. 
C> It applies a series of radix-based backward FFTs
C> (`radb2()`, `radb3()`, `radb4()`, `radb5()`, and `radbg()`)
C> to transform real input data back from the frequency domain to the time domain.
C>
C> @param n     Number of data points in the transform
C> @param c     Input/output array of real values; on input, it contains transformed data, 
C>              and on output, it stores the inverse-transformed real data
C> @param ch    Work array used for intermediate calculations
C> @param wa    Sine and cosine trigonometric table used for the transformation
C> @param ifac  Integer factorization array used to determine the radix decomposition
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RFFTB1 (N,C,CH,WA,IFAC)
      REAL       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = INT(IFAC(2))
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = INT(IFAC(K1+2))
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL RADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END

C> @brief Compute the real forward fast Fourier transform (FFT) using a mixed-radix algorithm
C>
C> This subroutine performs the FFT by applying a sequence of radix-based forward 
C> transformations (`radf2()`, `radf3()`, `radf4()`, `radf5()`, and `radfg()`).
C> It is a low-level routine used internally for computing the real-input FFT.
C>
C> @param n     Number of data points in the transform
C> @param c     Input/output array of real values; on input, it contains time-domain data,
C>              and on output, it stores the transformed frequency-domain data
C> @param ch    Work array used for intermediate calculations
C> @param wa    Sine and cosine trigonometric table used for the transformation
C> @param ifac  Integer factorization array used to determine the radix decomposition
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RFFTF1 (N,C,CH,WA,IFAC)
      REAL       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
      NF = INT(IFAC(2))
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = INT(IFAC(KH+3))
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL RADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END

C> @brief Compute the initialization factors for the real-input FFT
C>
C> This subroutine prepares the trigonometric table and factorization 
C> for the real-input Fast Fourier Transform (FFT). It determines the 
C> prime factorization of N and computes the necessary twiddle factors 
C> for efficient computation.
C>
C> @param n     The length of the input data array
C> @param wa    Output array containing computed trigonometric factors
C> @param ifac  Output integer array storing the factorization of `n` for FFT processing
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RFFTI1 (N,WA,IFAC)
      REAL       WA(*)      ,IFAC(*)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF ((J-4).LE.0) THEN
        GO TO 102
      ELSE
        GO TO 103
      ENDIF
  102 NTRY = INT(NTRYH(J))
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR.EQ.0) THEN
        GO TO 105
      ELSE
        GO TO 101
      ENDIF
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 6.28318530717959
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
!OCL NOVREC
      DO 110 K1=1,NFM1
         IP = INT(IFAC(K1+2))
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0
!OCL SCALAR
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

C> @brief Perform the backward radix-2 FFT step
C>
C> This subroutine computes one stage of the backward Fast Fourier Transform 
C> (FFT) using a radix-2 decomposition. It is used as part of the larger 
C> FFT process and applies twiddle factors to transform the input data.
C>
C> @param ido  The leading dimension, representing the number of data points in a transform section
C> @param l1   The number of radix-2 stages processed in this step
C> @param cc   Input complex array containing intermediate FFT data
C> @param ch   Output complex array storing transformed data
C> @param wa1  Twiddle factors array for weighting the FFT computation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADB2 (IDO,L1,CC,CH,WA1)
      DIMENSION       CC(IDO,2,L1)           ,CH(IDO,L1,2)           ,
     1                WA1(*)
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF (IDO.LT.2) THEN
        GO TO 107
      ELSEIF (IDO.EQ.2) THEN
        GO TO 105
      ELSE
        GO TO 102
      ENDIF
  102 IDP2 = IDO+2
!OCL NOVREC
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END

C> @brief Perform the backward radix-3 FFT step
C>
C> This subroutine computes one stage of the backward Fast Fourier Transform 
C> (FFT) using a radix-3 decomposition. It applies twiddle factors and combines
C> intermediate results to reconstruct transformed data.
C>
C> @param ido  The leading dimension, representing the number of data points in a transform section
C> @param l1   The number of radix-3 stages processed in this step
C> @param cc   Input complex array containing intermediate FFT data
C> @param ch   Output complex array storing transformed data
C> @param wa1  Twiddle factors array for the first rotation
C> @param wa2  Twiddle factors array for the second rotation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADB3 (IDO,L1,CC,CH,WA1,WA2)
      DIMENSION       CC(IDO,3,L1)           ,CH(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5,.866025403784439/
      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
!OCL NOVREC
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

C> @brief Perform the backward radix-4 FFT step
C>
C> This subroutine computes one stage of the backward Fast Fourier Transform 
C> (FFT) using a radix-4 decomposition. It applies the FFT formulae by processing 
C> the input data and performing necessary rotations and twiddle factor multiplications.
C> The result is stored in the output array `ch`.
C>
C> @param ido  The leading dimension, representing the number of data points in a transform section
C> @param l1   The number of radix-4 stages processed in this step
C> @param cc   Input complex array containing intermediate FFT data
C> @param ch   Output complex array storing transformed data
C> @param wa1  Twiddle factors array for the first rotation
C> @param wa2  Twiddle factors array for the second rotation
C> @param wa3  Twiddle factors array for the third rotation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      DIMENSION       CC(IDO,4,L1)           ,CH(IDO,L1,4)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
      DATA SQRT2 /1.414213562373095/
      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF (IDO.LT.2) THEN
        GO TO 107
      ELSEIF (IDO.EQ.2) THEN
        GO TO 105
      ELSE
        GO TO 102
      ENDIF
  102 IDP2 = IDO+2
!OCL NOVREC
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END

C> @brief Perform the backward radix-5 FFT step
C>
C> This subroutine computes one stage of the backward Fast Fourier Transform 
C> (FFT) using a radix-5 decomposition. It processes input data in sections, 
C> applies necessary rotations and twiddle factors, and stores the transformed 
C> results in the output array `ch`. The operation involves calculating 
C> intermediate results for each data point and applying appropriate scaling and 
C> twiddle factor multiplication.
C>
C> @param ido  The leading dimension, representing the number of data points in a transform section
C> @param l1   The number of radix-5 stages processed in this step
C> @param cc   Input complex array containing intermediate FFT data
C> @param ch   Output complex array storing transformed data
C> @param wa1  Twiddle factors array for the first rotation
C> @param wa2  Twiddle factors array for the second rotation
C> @param wa3  Twiddle factors array for the third rotation
C> @param wa4  Twiddle factors array for the fourth rotation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      DIMENSION       CC(IDO,5,L1)           ,CH(IDO,L1,5)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

C> @brief Computes the backward FFT stage using a generalized FFT algorithm.
C>
C> This subroutine performs a backward Fast Fourier Transform (FFT) step with 
C> a generalized radix-`ip` approach. It computes the necessary intermediate 
C> steps, applying the inverse FFT transformation to the complex input data 
C> stored in `cc`, and stores the result in the output arrays `ch` and `c2`.
C> Twiddle factors, stored in `wa`, are used for scaling and rotations at 
C> each stage of the computation.
C>
C> @param ido  The leading dimension for data points in a transform section
C> @param ip   The number of subgroups in the FFT (radix)
C> @param l1   The number of radix stages to process
C> @param idl1 The leading dimension for another input array
C> @param cc   Input complex array containing intermediate FFT data
C> @param c1   Array containing the first set of intermediate FFT results
C> @param c2   Array to store additional intermediate FFT results
C> @param ch   Output complex array to store the final transformed data
C> @param ch2  Output complex array storing intermediate results for the 
C>             second stage of transformation
C> @param wa   Array containing the twiddle factors used in the FFT
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
      DATA TPI/6.28318530717959/
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
!OCL NOVREC
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
!OCL NOVREC
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
!OCL NOVREC
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
!OCL NOVREC
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
!OCL NOVREC
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
!OCL NOVREC
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
!OCL NOVREC
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
!OCL NOVREC
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END

C> @brief Perform a radix-2 FFT step for even-odd decomposition of data
C>
C> This subroutine computes a single stage of the radix-2 Fast Fourier Transform 
C> (FFT) by separating the even and odd indexed values from the input data 
C> array `cc`, and stores the results in the output array `ch`. It processes 
C> intermediate data with the twiddle factors from `wa1` for scaling and rotation.
C>
C> @param ido  The size of the FFT sub-segment (number of data points to process)
C> @param l1   The number of stages to process in the FFT
C> @param cc   Input complex array of size `ido x l1 x 2`, containing the data 
C>             for the FFT transformation
C> @param ch   Output complex array of size `ido x 2 x l1`, where the 
C>             transformed data will be stored
C> @param wa1  Array of twiddle factors used to scale and rotate the data 
C>             during the FFT calculation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)
      DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           ,
     1                WA1(*)
      DO 101 K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
      IF (IDO.LT.2) THEN
        GO TO 107
      ELSEIF (IDO.EQ.2) THEN
        GO TO 105
      ELSE
        GO TO 102
      ENDIF
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END

C> @brief Performs a radix-3 FFT step for decomposing the data into three components.
C>
C> This subroutine computes a single stage of the radix-3 Fast Fourier Transform 
C> (FFT) by separating the input data `cc` into three components (one for each 
C> of the three terms in the FFT). The output data is stored in the array `ch`, 
C> with twiddle factors applied from the arrays `wa1` and `wa2`. This subroutine 
C> operates for various sizes of FFT sub-segments, controlled by the parameter 
C> `ido`.
C>
C> The subroutine handles different configurations based on the value of `ido`, 
C> which affects how the data is processed.
C>
C> @param ido  The size of the FFT sub-segment (number of data points to process)
C> @param l1   The number of stages to process in the FFT
C> @param cc   Input complex array of size `ido x l1 x 3`, containing the data 
C>              for the FFT transformation
C> @param ch   Output complex array of size `ido x 3 x l1`, where the transformed 
C>              data will be stored
C> @param wa1  Array of twiddle factors for the first part of the FFT calculation
C> @param wa2  Array of twiddle factors for the second part of the FFT calculation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADF3 (IDO,L1,CC,CH,WA1,WA2)
      DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           ,
     1                WA1(*)     ,WA2(*)
      DATA TAUR,TAUI /-.5,.866025403784439/
      DO 101 K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

C> @brief Perform a radix-4 FFT step for decomposing the data into four components
C>
C> This subroutine computes a single stage of the radix-4 Fast Fourier Transform 
C> (FFT) by separating the input data `cc` into four components (one for each 
C> of the four terms in the FFT). The output data is stored in the array `ch`, 
C> with twiddle factors applied from the arrays `wa1`, `wa2`, and `wa3`. This 
C> subroutine operates for various sizes of FFT sub-segments, controlled by the 
C> parameter `ido`.
C>
C> @param ido  The size of the FFT sub-segment (number of data points to process)
C> @param l1   The number of stages to process in the FFT
C> @param cc   Input complex array of size `ido x l1 x 4`, containing the data 
C>             for the FFT transformation
C> @param ch   Output complex array of size `ido x 4 x l1`, where the transformed 
C>             data will be stored
C> @param wa1  Array of twiddle factors for the first part of the FFT calculation
C> @param wa2  Array of twiddle factors for the second part of the FFT calculation
C> @param wa3  Array of twiddle factors for the third part of the FFT calculation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)
      DATA HSQT2 /.7071067811865475/
      DO 101 K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(IDO,4,K) = TR2-TR1
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  101 CONTINUE
      IF (IDO.LT.2) THEN
        GO TO 107
      ELSEIF (IDO.EQ.2) THEN
        GO TO 105
      ELSE
        GO TO 102
      ENDIF
  102 IDP2 = IDO+2
!OCL NOVREC
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1+CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(IDO,K,3)
         CH(1,4,K) = TI1+CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END

C> @brief Perform a radix-5 FFT step for decomposing the data into five components
C>
C> This subroutine computes a single stage of the radix-5 Fast Fourier Transform 
C> (FFT) by separating the input data `cc` into five components (one for each 
C> of the five terms in the FFT). The output data is stored in the array `ch`, 
C> with twiddle factors applied from the arrays `wa1`, `wa2`, `wa3`, and `wa4`. 
C> This subroutine operates for various sizes of FFT sub-segments, controlled by 
C> the parameter `ido`.
C>
C> @param ido  The size of the FFT sub-segment (number of data points to process)
C> @param l1   The number of stages to process in the FFT
C> @param cc   Input complex array of size `ido x l1 x 5`, containing the data 
C>              for the FFT transformation
C> @param ch   Output complex array of size `ido x 5 x l1`, where the transformed 
C>              data will be stored
C> @param wa1  Array of twiddle factors for the first part of the FFT calculation
C> @param wa2  Array of twiddle factors for the second part of the FFT calculation
C> @param wa3  Array of twiddle factors for the third part of the FFT calculation
C> @param wa4  Array of twiddle factors for the fourth part of the FFT calculation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           ,
     1                WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
      DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154,
     1-.809016994374947,.587785252292473/
      DO 101 K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END

C> @brief Compute a general radix-based FFT step
C>
C> This subroutine implements a mixed-radix Fast Fourier Transform (FFT) step,
C> decomposing the input data for an arbitrary radix `ip`. It reorganizes the
C> input data `cc` and applies twiddle factors from `wa` to generate the output
C> in `ch` and `ch2`.
C>
C> The transformation is influenced by `ido`, which controls the number of
C> data points per FFT sub-segment, and `l1`, which determines the number of
C> processing stages.
C>
C> @param ido  The number of data points per sub-segment
C> @param ip   The radix of the FFT step
C> @param l1   The number of stages to process
C> @param idl1 The transformed data length parameter, dependent on `l1`
C> @param cc   Input complex array of size `ido x ip x l1`, containing the data to transform
C> @param c1   Temporary working array for intermediate calculations
C> @param c2   Temporary working array for intermediate calculations
C> @param ch   Output complex array of size `ido x l1 x ip`, where transformed data is stored
C> @param ch2  Another working output array used for intermediate transformations
C> @param wa   Twiddle factor array used to apply phase shifts in the FFT calculation
C>
C> @author Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO
      SUBROUTINE RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          ,
     1                C1(IDO,L1,IP)          ,C2(IDL1,IP),
     2                CH2(IDL1,IP)           ,WA(*)
      DATA TPI/6.28318530717959/
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,K,J) = C1(1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(I,1,K) = CH(I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(I,1,K) = CH(I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
