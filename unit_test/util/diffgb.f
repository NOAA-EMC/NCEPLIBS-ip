C-----------------------------------------------------------------------
      PROGRAM DIFFGB
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C
C MAIN PROGRAM:  DIFFGB      COMPARE GRIB FILES
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: The command diffgb compares two GRIB files or alternatively
C   (-d option) prints statistics from a single GRIB file.  Unless
C   otherwise directed (-x option) the GRIB index files are also used
C   to speed the reading.  Either scalars or vectors (-v option) may be
C   compared.  There are five statistics types (-t option) for scalars
C   and three types for vectors.  These types are described later.  The
C   respective GRIB fields may be on virtually any grid; they may be
C   interpolated to a verification grid for computing statistics.
C   Unless the verification grid is directly specified (-g option)
C   the grid is the same as for the first GRIB field in GRIB file 1.
C   The interpolation type defaults to bilinear but may be specified
C   directly (-i option).  For certain grids, the statistics may be done
C   for a specific range of zonal Fourier modes (-z option).  The
C   statistics may be limited to specific fields (-k option).  They may
C   also be limited to a specified subgrid of the verification grid or
C   to a subrange of the input fields (-B and -b, -A, and -K options).
C   Unless otherwise directed (-w option), statistics are area-weighted.
C   Header information may be omitted (-h option) to ease further
C   processing.  Some statistics types have additional required or
C   optional options.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C COMMAND LINE OPTIONS:
C   -a "[<>] thresholds"
C      Threshold(s) for computing threat scores (difftype=3)
C      The thresholds may be led by a > or < to indicate
C      whether the thresholds represent minimum or maximum
C      values, respectively (defaults to >0.).
C
C   -A "<> mapthreshold"
C      Inequality and threshold used in determining
C      where on the map statistics will be computed.
C      Statistics are computed only where the given 
C      map field is on the correct side of the threshold.
C      The mapthresh defaults to '>-1.e30'; in this case,
C      only the map field's bitmap will limit the domain.
C
C   -b mapindex   
C      Optional index file used to get the map field.
C
C   -B mapgrib    
C      GRIB file used to get the map field.  The map field
C      is read from the GRIB file and compared to the
C      map threshold to determine for which region on the map
C      the statistics will be computed.  mapgrib can be the
C      name of an actual GRIB file (in which case the index
C      file may be specified with the -b option) or it can
C      be either '-1', '-2', or '-3'.  If mapgrib is '-1',
C      then GRIB file 1 (first positional argument) is used.
C      If mapgrib is '-2', then GRIB file 2 is used.  If
C      mapgrib is '-3', then the climatology is used.
C      The -K option specifies which field to read from
C      the mapgrib GRIB file.  If -K is not specified, then
C      the first field is taken if mapgrib is an actual file;
C      otherwise one of the fields currently being compared
C      is taken as the map field.
C
C   -c climindex   
C      Optional index file for the climatology GRIB file.
C
C   -C climgrib    
C      GRIB file used to get the climatology when difftype=5
C      or 'anomaly'.  At present the climatology must be
C      expressed in monthly averages.
C
C   -d
C      Turns off comparison and activates statistics
C      for a single GRIB file.
C
C   -g "grid [kgds]"
C      Verification grid identification.  If grid=-1
C      (the default), then the grid is taken from the first
C      GRIB field in GRIB file 1.  If grid=-2,
C      then the grid is taken from the first GRIB field
C      in GRIB file 2.  IF grid=-3, then the grid is taken
C      from the first GRIB field in the climatology.
C      If 0<grid<255, then grid designates an NCEP grid.
C      If grid=255, then the grid must be specified by the
C      full set of kgds parameters determining a GRIB GDS
C      (grid description section) in the W3FI63 format.
C
C   -h    
C      Turns off printing of header and summary information
C      to ease further processing.
C
C   -i "ip [ipopts]"
C      Interpolation options.  The default is bilinear
C      interpolation (ip=0).  Other interpolation options
C      are bicubic (ip=1), neighbor (ip=2), budget (ip=3),
C      and spectral (ip=4).  See the documentation for iplib
C      for further details about interpolation.
C
C   -k "kpds"
C      Full set of kpds parameters determing a GRIB PDS
C      (product definition section) in the W3FI63 format
C      determining the field(s) from GRIB file 1
C      for which statistics are computed.  Note that
C      kpds(5) is the parameter indicator (PDS octet 9).
C      A wildcard is specified by -1 (the defaults).
C      If the -k is not specified, then diffgb will attempt
C      to compare every field in GRIB file 1.
C
C      Note that this option does not apply to GRIB file 2.
C      After each field from GRIB file 1 has been found,
C      GRIB file 2 is searched from the beginning for a field
C      with the same parameter table version number,
C      parameter identifier, level type and level value
C      (PDS octets 4 and 9-12) in order to make a comparison.
C
C   -K "mapkpds"
C      Full set of kpds parameters determing a GRIB PDS
C      (product definition section) in the W3FI63 format
C      determining the map field to be used to determine
C      where on the map statistics will be computed.  
C      A wildcard is specified by -1 (the defaults).
C
C   -t difftype
C      Statistics type.  Possible values are 1, 2, 3, 4, 5, 
C      min, max, avg, sig, thr, cor, rms, or ano.
C
C      If difftype=1 or 'min' or 'max' (default),
C      diffgb outputs each field's minimum and maximum
C      as well as the number of differences and the
C      maximum difference (along with its array index).
C      For vectors, diffgb outputs the minimum and maximum
C      of each field's vector magnitude as well as the number
C      of differences and the maximum magnitude of the vector
C      difference (along with its array index).
C
C      If difftype=2 or 'avg' or 'sig',
C      diffgb outputs each field's average and standard
C      deviation as well as the mean difference and standard
C      deviation of the difference.  For vectors, diffgb
C      outputs the average and standard deviation of each
C      field's vector magnitude as well as the mean
C      difference in magnitude and the RMS vector difference.
C
C      If difftype=3 or 'thr',
C      for each threshold given by the -a option,
C      diffgb outputs each field's percentage coverage
C      and hit rate with respect to the other field
C      as well as the simple and equitable threat scores.
C      For vectors, diffgb outputs equivalent statistics
C      with respect to the field's vector magnitude.
C
C      If difftype=4 or 'cor' or 'rms',
C      diffgb outputs each field's average and standard
C      deviation as well as the correlation (between each
C      field's deviation from its own respective average)
C      and the RMS of the difference.
C
C      If difftype=5 or 'ano',
C      a climatology must be specified with the -C option
C      with which diffgb gets each field's climate anomalies
C      and outputs each field's anomaly average and RMS
C      as well as the anomaly correlation and the RMS of the
C      difference.  This option is not valid for vectors.
C
C   -v
C      Turns on vector rather than scalar comparison.  The
C      vector statistics are similar but separate from
C      the scalar statistics.  No correlation statistics
C      can be computed for vectors.  The horizontal vector
C      components are assumed to have consecutive parameter
C      identifications.  If the -k option specifies the
C      parameter identification kpds(5), then parameters
C      kpds(5) and kpds(5)+1 are taken as the respective
C      vector components; otherwise the wind components
C      (parameters 33 and 34) are taken.
C
C   -w
C      Turns off area-weighting.  The statistics are computed
C      weighting each point in the verification grid equally.
C
C   -x 
C      Turns off the use of index files.  The index records
C      are then extracted from the GRIB file(s), which
C      will increase the time taken by diffgb.
C
C   -z "[kz1] kz2"
C      Zonal Fourier mode range for which the statistics
C      are computed.  If kz1 is not specified, it defaults
C      to 0.  The zonal Fourier filter is performed on the
C      verification grid.  Therefore this option can only
C      be specified for an appropriate verification grid,
C      i.e. a cylindrical grid with full cyclic zonal range.
C
C
C INPUT FILES:
C   UNIT   11    GRIB FILE 1
C   UNIT   12    GRIB FILE 2
C   UNIT   13    CLIMATOLOGY GRIB FILE
C   UNIT   14    MAP GRIB FILE
C   UNIT   31    GRIB INDEX FILE 1
C   UNIT   32    GRIB INDEX FILE 2
C   UNIT   33    CLIMATOLOGY GRIB INDEX FILE
C   UNIT   34    MAP GRIB INDEX FILE
C
C SUBPROGRAMS CALLED:
C   FILENV
C   IARGC
C   GETARG
C   ERRMSG
C   EUSAGE
C   EXIT
C   FPARSEI
C   FPARSER
C   BAOPENR
C   GETGBMH
C   LENGDSF
C   DFGB
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CHARACTER*256 CARG,CG1,CX1,CG2,CX2,CGC,CXC,CGB,CXB
      CHARACTER*3 CTHREE
      INTEGER KGDSI(200),IPOPT(20),KARG(100)
      INTEGER JPDS1(200),JPDSB(200)
      REAL A3(100),RARG(100)
      PARAMETER(MBUF=256*1024)
      CHARACTER CBUF1(MBUF),CBUF2(MBUF),CBUFC(MBUF),CBUFB(MBUF)
      INTEGER JPDS(200),JGDS(200)
      INTEGER KPDS1(200),KGDS1(200),KPDS2(200),KGDS2(200)
      INTEGER KPDSC(200),KGDSC(200),KPDSB(200),KGDSB(200),KGDSIE(200)
      INTEGER KGDS1E(200),KGDS2E(200),KGDSCE(200),KGDSBE(200)
      CHARACTER*400 GDS
      DATA IGI/-1/,KGDSIE/19*0,255,180*0/
      DATA IP/0/,IPOPT/20*-1/
      DATA JPDS1/200*-1/,JPDSB/200*-1/
      DATA LD/1/,LH/1/,LT/1/,LV/1/,LW/1/,LX/1/,KZ1/-1/,KZ2/-2/
      DATA L3/1/,N3/1/,A3/0.,99*1.E30/
      DATA CGC/' '/,CXC/' '/,IGC/0/
      DATA LAB/1/,AB/-1.E30/
      DATA CGB/' '/,CXB/' '/
      DATA MI/1/,M1/1/,M2/1/,MC/1/,MB/1/
      DATA MIE/1/,M1E/1/,M2E/1/,MBE/1/,MCE/1/
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PARSE COMMAND LINE OPTIONS
      CALL FILENV
      NARG=IARGC()
      IARG=1
      LSTOPT=0
      DOWHILE(IARG.LE.NARG.AND.LSTOPT.EQ.0)
        CALL GETARG(IARG,CARG)
        LARG=LEN_TRIM(CARG)
        IARG=IARG+1
        IF(CARG(1:1).NE.'-') THEN
          LSTOPT=1
          IARG=IARG-1
        ELSEIF(LARG.EQ.1) THEN
          CALL ERRMSG('diffgb: invalid option -')
          CALL EUSAGE
          CALL EXIT(1)
        ELSE
          L=2
          DOWHILE(L.LE.LARG)
            IF(CARG(L:L).EQ.'-') THEN
              LSTOPT=1
            ELSEIF(CARG(L:L).EQ.'A') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              IF(CARG(L+1:L+1).EQ.'>') THEN
                LAB=1
                L=L+1
              ELSEIF(CARG(L+1:L+1).EQ.'<') THEN
                LAB=-1
                L=L+1
              ELSE
                CALL ERRMSG('diffgb: invalid threshold '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              CALL FPARSER(CARG(L+1:LARG),1,AB)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'a') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              IF(CARG(L+1:L+1).EQ.'>') THEN
                L3=1
                L=L+1
              ELSEIF(CARG(L+1:L+1).EQ.'<') THEN
                L3=-1
                L=L+1
              ENDIF
              CALL FPARSER(CARG(L+1:LARG),100,A3)
              N3=1
              DO I3=2,100
                IF(A3(I3).NE.1.E30) N3=I3
              ENDDO
              L=LARG
            ELSEIF(CARG(L:L).EQ.'B') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              LCGB=LARG-L
              CGB=CARG(L+1:LARG)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'b') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              LCXB=LARG-L
              CXB=CARG(L+1:LARG)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'C') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              LCGC=LARG-L
              CGC=CARG(L+1:LARG)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'c') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              LCXC=LARG-L
              CXC=CARG(L+1:LARG)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'d') THEN
              LD=0
            ELSEIF(CARG(L:L).EQ.'g') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              KARG(1)=IGI
              KARG(2:100)=KGDSIE(1:99)
              CALL FPARSEI(CARG(L+1:LARG),100,KARG)
              IGI=KARG(1)
              IF(IGI.GT.0.AND.IGI.LT.255) THEN
!                CALL MAKGDS(IGI,KGDSIE,GDS,LGDS,IRET)
                  stop 2
                IF(IRET.NE.0) IGI=-1
              ELSEIF(IGI.EQ.255) THEN
                KGDSIE(1:99)=KARG(2:100)
              ENDIF
              IF(IGI.LT.-3.OR.IGI.EQ.0.OR.IGI.GT.255) THEN
                CALL ERRMSG('diffgb: invalid verification grid '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              MI=LENGDSF(KGDSIE,KGDSI)
              IF(MI.LE.0) THEN
                CALL ERRMSG('diffgb: unsupported verification grid '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'h') THEN
              LH=0
            ELSEIF(CARG(L:L).EQ.'i') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              KARG(1)=IP
              KARG(2:21)=IPOPT
              CALL FPARSEI(CARG(L+1:LARG),21,KARG)
              IP=KARG(1)
              IPOPT=KARG(2:21)
              L=LARG
            ELSEIF(CARG(L:L).EQ.'K') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              CALL FPARSEI(CARG(L+1:LARG),100,JPDSB)
              IF(JPDSB(5).EQ.0) THEN
                CALL ERRMSG('diffgb: invalid PDS parms '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'k') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              CALL FPARSEI(CARG(L+1:LARG),100,JPDS1)
              IF(JPDS1(5).EQ.0) THEN
                CALL ERRMSG('diffgb: invalid PDS parms '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'t') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              IF(LARG-L.GE.3) THEN
                CTHREE=CARG(L+1:L+3)
                IF(CTHREE.EQ.'min'.OR.CTHREE.EQ.'max'.OR.
     &             CTHREE.EQ.'MIN'.OR.CTHREE.EQ.'MAX') THEN
                  LT=1
                ELSEIF(CTHREE.EQ.'avg'.OR.CTHREE.EQ.'sig'.OR.
     &                 CTHREE.EQ.'AVG'.OR.CTHREE.EQ.'SIG') THEN
                  LT=2
                ELSEIF(CTHREE.EQ.'thr'.OR.
     &                 CTHREE.EQ.'THR') THEN
                  LT=3
                ELSEIF(CTHREE.EQ.'cor'.OR.CTHREE.EQ.'rms'.OR.
     &                 CTHREE.EQ.'COR'.OR.CTHREE.EQ.'RMS') THEN
                  LT=4
                ELSEIF(CTHREE.EQ.'ano'.OR.
     &                 CTHREE.EQ.'ANO') THEN
                  LT=5
                ELSE
                  LT=0
                ENDIF
              ELSE
                CALL FPARSEI(CARG(L+1:LARG),1,LT)
              ENDIF
              IF(LT.LT.1.OR.LT.GT.9) THEN
                CALL ERRMSG('diffgb: invalid difference type '//
     &                      CARG(L+1:LARG))
                CALL EUSAGE
                CALL EXIT(1)
              ENDIF
              L=LARG
            ELSEIF(CARG(L:L).EQ.'v') THEN
              LV=2
            ELSEIF(CARG(L:L).EQ.'w') THEN
              LW=0
            ELSEIF(CARG(L:L).EQ.'x') THEN
              LX=0
            ELSEIF(CARG(L:L).EQ.'z') THEN
              IF(L.EQ.LARG) THEN
                L=0
                CALL GETARG(IARG,CARG)
                LARG=LEN_TRIM(CARG)
                IARG=IARG+1
              ENDIF
              KARG(1)=KZ1
              KARG(2)=KZ2
              CALL FPARSEI(CARG(L+1:LARG),2,KARG)
              IF(KARG(2).GE.0) THEN
                KZ1=MIN(KARG(1),KARG(2))
                KZ2=MAX(KARG(1),KARG(2))
              ELSE
                KZ1=0
                KZ2=KARG(1)
              ENDIF
              L=LARG
            ELSE
              CALL ERRMSG('diffgb: invalid option '//CARG(L:L))
              CALL EUSAGE
              CALL EXIT(1)
            ENDIF
            L=L+1
          ENDDO
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PARSE COMMAND LINE POSITIONAL ARGUMENTS
      NXARG=(LD+1)*(LX+1)
      IF(NARG-IARG+1.NE.NXARG) THEN
        CALL ERRMSG('diffgb: incorrect number of arguments')
        CALL EUSAGE
        CALL EXIT(NXARG)
      ENDIF
      LG1=11
      CALL GETARG(IARG,CG1)
      LCG1=LEN_TRIM(CG1)
      IARG=IARG+1
      CALL BAOPENR(LG1,CG1(1:LCG1),IRETBA)
      IF(IRETBA.NE.0) THEN
        CALL ERRMSG('diffgb:  error accessing file '//CG1(1:LCG1))
        CALL EXIT(8)
      ENDIF
      IF(LX.GT.0) THEN
        LX1=31
        CALL GETARG(IARG,CX1)
        LCX1=LEN_TRIM(CX1)
        IARG=IARG+1
        CALL BAOPENR(LX1,CX1(1:LCX1),IRETBA)
        IF(IRETBA.NE.0) THEN
          CALL ERRMSG('diffgb:  error accessing file '//CX1(1:LCX1))
          CALL EXIT(8)
        ENDIF
      ELSE
        LX1=0
      ENDIF
      IF(LD.GT.0) THEN
        LG2=12
        CALL GETARG(IARG,CG2)
        LCG2=LEN_TRIM(CG2)
        IARG=IARG+1
        CALL BAOPENR(LG2,CG2(1:LCG2),IRETBA)
        IF(IRETBA.NE.0) THEN
          CALL ERRMSG('diffgb:  error accessing file '//CG2(1:LCG2))
          CALL EXIT(8)
        ENDIF
        IF(LX.NE.0) THEN
          LX2=32
          CALL GETARG(IARG,CX2)
          LCX2=LEN_TRIM(CX2)
          IARG=IARG+1
          CALL BAOPENR(LX2,CX2(1:LCX2),IRETBA)
          IF(IRETBA.NE.0) THEN
            CALL ERRMSG('diffgb:  error accessing file '//CX2(1:LCX2))
            CALL EXIT(8)
          ENDIF
        ELSE
          LX2=0
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  OPEN GRIB FILE 1
      JPDS=-1
      CALL GETGBMH(LG1,LX1,-1,JPDS,JGDS,
     &             MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &             KG1,M1E,KR1,KPDS1,KGDS1E,IRET1)
      IF(IRET1.NE.0) THEN
        CALL ERRMSG('diffgb: error reading file 1')
        CALL EXIT(IRET1)
      ENDIF
      ID1=MOD(KPDS1(8),100)*1000000+KPDS1(9)*10000+
     &    KPDS1(10)*100+KPDS1(11)
      IE1=KPDS1(14)
      IF1=KPDS1(14)
      IF(KPDS1(16).GE.2.AND.KPDS1(16).LE.5) IF1=KPDS1(15)
      IG1=KPDS1(3)
      M1=LENGDSF(KGDS1E,KGDS1)
      IF(M1.LE.0) THEN
        CALL ERRMSG('diffgb: error transforming irregular grid')
        CALL EXIT(4)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  OPEN GRIB FILE 2
      IF(LD.GT.0) THEN
        CALL GETGBMH(LG2,LX2,-1,JPDS,JGDS,
     &               MBUF,CBUF2,NLEN2,NNUM2,MNUM2,
     &               KG2,M2E,KR2,KPDS2,KGDS2E,IRET2)
        IF(IRET2.NE.0) THEN
          CALL ERRMSG('diffgb: error reading file 2')
          CALL EXIT(IRET2)
        ENDIF
        ID2=MOD(KPDS2(8),100)*1000000+KPDS2(9)*10000+
     &      KPDS2(10)*100+KPDS2(11)
        IE2=KPDS2(14)
        IF2=KPDS2(14)
        IF(KPDS2(16).GE.2.AND.KPDS2(16).LE.5) IF2=KPDS2(15)
        IG2=KPDS2(3)
        M2=LENGDSF(KGDS2E,KGDS2)
        IF(M2.LE.0) THEN
          CALL ERRMSG('diffgb: error transforming irregular grid')
          CALL EXIT(4)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  OPEN CLIMATOLOGY GRIB FILE
      IF(LT.EQ.5) THEN
        LGC=13
        IF(CGC(1:1).EQ.' ') THEN
          CALL ERRMSG('diffgb: unspecified climatology')
          CALL EXIT(1)
        ENDIF
        CALL BAOPENR(LGC,CGC(1:LCGC),IRETBA)
        IF(IRETBA.NE.0) THEN
          CALL ERRMSG('diffgb:  error accessing file '//CGC(1:LCGC))
          CALL EXIT(8)
        ENDIF
        IF(CXC(1:1).NE.' ') THEN
          LXC=33
          CALL BAOPENR(LXC,CXC(1:LCXC),IRETBA)
          IF(IRETBA.NE.0) THEN
            CALL ERRMSG('diffgb:  error accessing file '//CXC(1:LCXC))
            CALL EXIT(8)
          ENDIF
        ELSE
          LXC=0
        ENDIF
        JPDS=-1
        CALL GETGBMH(LGC,LXC,-1,JPDS,JGDS,
     &               MBUF,CBUFC,NLENC,NNUMC,MNUMC,
     &               KGC,MCE,KRC,KPDSC,KGDSCE,IRETC)
        IF(IRETC.NE.0) THEN
          CALL ERRMSG('diffgb: error reading climatology')
          CALL EXIT(IRETC)
        ENDIF
        IGC=KPDSC(3)
        MC=LENGDSF(KGDSCE,KGDSC)
        IF(MC.LE.0) THEN
          CALL ERRMSG('diffgb: error transforming irregular grid')
          CALL EXIT(4)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  OPEN MAP FILE
      IF(CGB.NE.' ') THEN
        IF(CGB(1:2).EQ.'-1') THEN
          IF(JPDSB(5).EQ.-1) THEN
            JB=1
          ELSE
            JB=4
            LGB=LG1
            LXB=LX1
            MBE=M1E
            MB=M1
            IGB=IG1
            CBUFB=CBUF1
            NLENB=NLEN1
            NNUMB=NNUM1
            MNUMB=MNUM1
            IGB=IG1
          ENDIF
        ELSEIF(CGB(1:2).EQ.'-2'.AND.LD.GT.0) THEN
          IF(JPDSB(5).EQ.-1) THEN
            JB=2
          ELSE
            JB=4
            LGB=LG2
            LXB=LX2
            MBE=M2E
            MB=M2
            IGB=IG2
            CBUFB=CBUF2
            NLENB=NLEN2
            NNUMB=NNUM2
            MNUMB=MNUM2
          ENDIF
        ELSEIF(CGB(1:2).EQ.'-C'.AND.LT.EQ.5) THEN
          IF(JPDSB(5).EQ.-1) THEN
            JB=3
          ELSE
            JB=4
            LGB=LGC
            LXB=LXC
            MBE=MCE
            MB=MC
            IGB=IGC
            CBUFB=CBUFC
            NLENB=NLENC
            NNUMB=NNUMC
            MNUMB=MNUMC
          ENDIF
        ELSE
          JB=4
          LGB=14
          CALL BAOPENR(LGB,CGB(1:LCGB),IRETBA)
          IF(IRETBA.NE.0) THEN
            CALL ERRMSG('diffgb:  error accessing file '//CGB(1:LCGB))
            CALL EXIT(8)
          ENDIF
          IF(CXB(1:1).NE.' ') THEN
            LXB=34
            CALL BAOPENR(LXB,CXB(1:LCXB),IRETBA)
            IF(IRETBA.NE.0) THEN
              CALL ERRMSG('diffgb:  error accessing file '//CXB(1:LCXB))
              CALL EXIT(8)
            ENDIF
          ELSE
            LXB=0
          ENDIF
          JPDS=-1
          CALL GETGBMH(LGB,LXB,-1,JPDS,JGDS,
     &                 MBUF,CBUFB,NLENB,NNUMB,MNUMB,
     &                 KGB,MBE,KRB,KPDSB,KGDSBE,IRETB)
          IF(IRETB.NE.0) THEN
            CALL ERRMSG('diffgb: error reading bitmap')
            CALL EXIT(IRETB)
          ENDIF
          IGB=KPDSB(3)
          MB=LENGDSF(KGDSBE,KGDSB)
          IF(MB.LE.0) THEN
            CALL ERRMSG('diffgb: error transforming irregular grid')
            CALL EXIT(4)
          ENDIF
        ENDIF
      ELSE
        JB=0
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  DETERMINE VERIFICATION GRID
      IF(IGI.EQ.-1) THEN
        IGI=IG1
        KGDSI=KGDS1
        MI=M1
      ELSEIF(IGI.EQ.-2.AND.LD.GT.0) THEN
        IGI=IG2
        KGDSI=KGDS2
        MI=M2
      ELSEIF(IGI.EQ.-3.AND.LT.EQ.5) THEN
        IGI=IGC
        KGDSI=KGDSC
        MI=MC
      ELSEIF(IGI.EQ.-4.AND.JB.EQ.4) THEN
        IGI=IGB
        KGDSI=KGDSB
        MI=MB
      ENDIF
      IF(IGI.LE.0) THEN
        CALL ERRMSG('diffgb: improper verification grid')
        CALL EXIT(4)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INITIALIZE FOURIER FILTER
      IF(KZ1.LE.KZ2) THEN
        CALL ZFF0(KZ1,KZ2,KGDSI,IRZ)
        IF(IRZ.NE.0) THEN
          CALL ERRMSG('diffgb: improper zonal fourier filter')
          CALL EXIT(IRZ)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  SET DIMENSIONS AND GO
      MF=MAX(M1,M2,MC,MB)
      MG=MAX(M1,M2)
      IF(LV.EQ.1) MG=1
      MGI=MI
      IF(LV.EQ.1)  MGI=1
      M2I=MI
      IF(LD.EQ.0)  M2I=1
      MCI=MI
      IF(IGC.EQ.0)  MCI=1
      MBI=MI
      IF(JB.EQ.0)  MBI=1
      IF(LV.EQ.2.AND.JPDS1(5).EQ.-1) JPDS1(5)=33
      CALL DFGB(LG1,LX1,M1E,M1,IG1,CBUF1,NLEN1,NNUM1,MNUM1,
     &          LG2,LX2,M2E,M2,IG2,CBUF2,NLEN2,NNUM2,MNUM2,
     &          ID1,IE1,IF1,ID2,IE2,IF2,
     &          MBUF,MF,MG,MI,MGI,M2I,MCI,MBI,
     &          IGI,KGDSI,LD,LH,LT,LV,LW,
     &          IP,IPOPT,JPDS1,KZ1,KZ2,
     &          JPDSB,JB,LAB,AB,
     &          LGB,LXB,MBE,MB,IGB,CBUFB,NLENB,NNUMB,MNUMB,
     &          L3,N3,A3,
     &          LGC,LXC,MCE,MC,IGC,CBUFC,NLENC,NNUMC,MNUMC)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE EUSAGE
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    EUSAGE      PRINT PROPER USAGE TO STDERR
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT PROPER USAGE TO STDERR.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL EUSAGE
C
C SUBPROGRAMS CALLED:
C   ERRMSG
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CALL ERRMSG('Usage: diffgb'//
     & ' [-h] [-v] [-w] [-g "grid [kgds]"] [-i "ip [ipopts]"]')
      CALL ERRMSG('             '//
     & ' [-k "kpds"] [-t difftype] [-z "[kz1] kz2"]')
      CALL ERRMSG('             '//
     & ' [-B mapgrib [-b mapindex] [-A "<> mapthreshold"]'//
     & ' [-K "mapkpds"]]')
      CALL ERRMSG('       then either:')
      CALL ERRMSG('             '//
     & ' grib1 index1 grib2 index2')
      CALL ERRMSG('            or:')
      CALL ERRMSG('             '//
     & ' -d grib1 index1')
      CALL ERRMSG('            or:')
      CALL ERRMSG('             '//
     & ' -x grib1 grib2')
      CALL ERRMSG('            or:')
      CALL ERRMSG('             '//
     & ' -dx grib1')
      CALL ERRMSG('       also if difftype=3:')
      CALL ERRMSG('             '//
     & ' [-a "[<>] thresholds"]')
      CALL ERRMSG('            if difftype=5:')
      CALL ERRMSG('             '//
     & ' -C climgrib [-c climindex]')
      END
C-----------------------------------------------------------------------
      SUBROUTINE DFGB(LG1,LX1,M1E,M1,IG1,CBUF1,NLEN1,NNUM1,MNUM1,
     &                LG2,LX2,M2E,M2,IG2,CBUF2,NLEN2,NNUM2,MNUM2,
     &                ID1,IE1,IF1,ID2,IE2,IF2,
     &                MBUF,MF,MG,MI,MGI,M2I,MCI,MBI,
     &                IGI,KGDSI,LD,LH,LT,LV,LW,
     &                IP,IPOPT,JPDS1,KZ1,KZ2,
     &                JPDSB,JB,LAB,AB,
     &                LGB,LXB,MBE,MB,IGB,CBUFB,NLENB,NNUMB,MNUMB,
     &                L3,N3,A3,
     &                LGC,LXC,MCE,MC,IGC,CBUFC,NLENC,NNUMC,MNUMC)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DFGB        COMPARE GRIB FILES
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: COMPARE GRIB FILES.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DFGB(LG1,LX1,M1E,M1,IG1,CBUF1,NLEN1,NNUM1,MNUM1,
C    &                LG2,LX2,M2E,M2,IG2,CBUF2,NLEN2,NNUM2,MNUM2,
C    &                ID1,IF1,ID2,IF2,
C    &                MBUF,MF,MG,MI,MGI,M2I,MCI,MBI,
C    &                IGI,KGDSI,LD,LH,LT,LV,LW,
C    &                IP,IPOPT,JPDS1,KZ1,KZ2,
C    &                JPDSB,JB,LAB,AB,
C    &                LGB,LXB,MBE,MB,IGB,CBUFB,NLENB,NNUMB,MNUMB,
C    &                L3,N3,A3,
C    &                LGC,LXC,MCE,MC,IGC,CBUFC,NLENC,NNUMC,MNUMC)
C   INPUT ARGUMENTS:
C     LG1          INTEGER UNIT NUMBER FOR GRIB FILE 1
C     LX1          INTEGER UNIT NUMBER FOR GRIB INDEX FILE 1
C     M1E          INTEGER DIMENSION OF IRREGULAR GRIB FIELD 1
C     M1           INTEGER DIMENSION OF GRIB FIELD 1
C     IG1          INTEGER GRID IDENTIFICATION 1
C     CBUF1        CHARACTER (MBUF) INDEX BUFFER 1
C     NLEN1        INTEGER RECORD LENGTH OF INDEX BUFFER 1
C     NNUM1        INTEGER NUMBER OF RECORDS IN INDEX BUFFER 1
C     NLEN1        INTEGER LENGTH OF EACH INDEX RECORD 1
C     NNUM1        INTEGER NUMBER OF INDEX RECORDS 1
C     MNUM1        INTEGER NUMBER OF INDEX RECORDS 1 SKIPPED
C     LG2          INTEGER UNIT NUMBER FOR GRIB FILE 2
C     LX2          INTEGER UNIT NUMBER FOR GRIB INDEX FILE 2
C     M2E          INTEGER DIMENSION OF IRREGULAR GRIB FIELD 2
C     M2           INTEGER DIMENSION OF GRIB FIELD 2
C     IG2          INTEGER GRID IDENTIFICATION 2
C     CBUF2        CHARACTER (MBUF) INDEX BUFFER 2
C     NLEN2        INTEGER RECORD LENGTH OF INDEX BUFFER 2
C     NNUM2        INTEGER NUMBER OF RECORDS IN INDEX BUFFER 2
C     NLEN2        INTEGER LENGTH OF EACH INDEX RECORD 2
C     NNUM2        INTEGER NUMBER OF INDEX RECORDS 2
C     MNUM2        INTEGER NUMBER OF INDEX RECORDS 2 SKIPPED
C     ID1          INTEGER INTITIAL TIME 1
C     IE1          INTEGER FORECAST HOUR 1
C     IF1          INTEGER FORECAST HOUR 1
C     ID2          INTEGER INTITIAL TIME 2
C     IE2          INTEGER FORECAST HOUR 2
C     IF2          INTEGER FORECAST HOUR 2
C     MBUF         INTEGER DIMENSION OF INDEX BUFFERS
C     MF           INTEGER DIMENSION OF FIELD
C     MG           INTEGER DIMENSION OF SECOND FIELD
C     MI           INTEGER DIMENSION OF VERIFICATION GRID
C     MGI          INTEGER DIMENSION OF VERIFICATION GRID (SECOND FIELD)
C     M2I          INTEGER DIMENSION OF VERIFICATION GRID (FILE 2)
C     MCI          INTEGER DIMENSION OF VERIFICATION GRID (CLIMATOLOGY)
C     MBI          INTEGER DIMENSION OF VERIFICATION GRID (MAP)
C     IGI          INTEGER VERIFICATION GRID IDENTIFICATION
C     KGDSI        INTEGER (200) VERIFICATION GRID PARAMETERS
C     LD           INTEGER FLAG FOR D OPTION
C     LT           INTEGER FLAG FOR T OPTION
C     LV           INTEGER FLAG FOR V OPTION
C     LW           INTEGER FLAG FOR W OPTION
C     IP           INTEGER INTERPOLATION TYPE
C     IPOPT        INTEGER (20) INTERPOLATION OPTIONS
C     JPDS1        INTEGER (100) KPDS SEARCH OPTIONS
C     KZ1          INTEGER FOURIER MODE RANGE 1
C     KZ2          INTEGER FOURIER MODE RANGE 2
C     JPDSB        INTEGER (100) KPDS SEARCH OPTIONS (MAP)
C     JB           INTEGER FLAG FOR MAP OPTIION
C     LAB          INTEGER FLAG FOR MAP THRESHOLD INEQUALITY
C     AB           REAL MAP THRESHOLD
C     LGB          INTEGER UNIT NUMBER FOR GRIB FILE MAP
C     LXB          INTEGER UNIT NUMBER FOR GRIB INDEX FILE MAP
C     MBE          INTEGER DIMENSION OF IRREGULAR GRIB FIELD MAP
C     MB           INTEGER DIMENSION OF GRIB FIELD MAP
C     IGB          INTEGER GRID IDENTIFICATION MAP
C     CBUFB        CHARACTER (MBUF) INDEX BUFFER MAP
C     NLENB        INTEGER RECORD LENGTH OF INDEX BUFFER MAP
C     NNUMB        INTEGER NUMBER OF RECORDS IN INDEX BUFFER MAP
C     NLENB        INTEGER LENGTH OF EACH INDEX RECORD MAP
C     NNUMB        INTEGER NUMBER OF INDEX RECORDS MAP
C     MNUMB        INTEGER NUMBER OF INDEX RECORDS MAP SKIPPED
C     L3           INTEGER FLAG FOR THRESHOLD INEQUALITY
C     N3           INTEGER NUMBER OF THRESHOLDS
C     A3           REAL (N3) THRESHOLDS
C     LGC          INTEGER UNIT NUMBER FOR GRIB FILE CLIM
C     LXC          INTEGER UNIT NUMBER FOR GRIB INDEX FILE CLIM
C     MCE          INTEGER DIMENSION OF IRREGULAR GRIB FIELD CLIM
C     MC           INTEGER DIMENSION OF GRIB FIELD CLIM
C     IGC          INTEGER GRID IDENTIFICATION CLIM
C     CBUFC        CHARACTER (MBUF) INDEX BUFFER CLIM
C     NLENC        INTEGER RECORD LENGTH OF INDEX BUFFER CLIM
C     NNUMC        INTEGER NUMBER OF RECORDS IN INDEX BUFFER CLIM
C     NLENC        INTEGER LENGTH OF EACH INDEX RECORD CLIM
C     NNUMC        INTEGER NUMBER OF INDEX RECORDS CLIM
C     MNUMC        INTEGER NUMBER OF INDEX RECORDS CLIM SKIPPED
C
C SUBPROGRAMS CALLED:
C   GETGRIB
C   IPOLATES
C   IPOLATEV
C   ERRMSG
C   EXIT
C   ZFF1
C   GETCLIM
C   GDSAWT
C   DGBHS1
C   DGBHS2
C   DGBHS3
C   DGBHS4
C   DGBHS5
C   DGBHV1
C   DGBHV2
C   DGBHV3
C   DGBFS1
C   DGBFS2
C   DGBFS3
C   DGBFS4
C   DGBFS5
C   DGBFV1
C   DGBFV2
C   DGBFV3
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CHARACTER CBUF1(MBUF),CBUF2(MBUF),CBUFC(MBUF),CBUFB(MBUF)
      INTEGER KGDSI(200)
      INTEGER IPOPT(20)
      INTEGER JPDS1(200),JPDSB(200)
      REAL A3(N3)
      INTEGER JPDS(200),JGDS(200)
      INTEGER KPDS1(200),KGDS1(200),KPDS2(200),KGDS2(200)
      INTEGER KPDSC(200),KGDSC(200),KPDSB(200),KGDSB(200)
      LOGICAL*1 LR(MF),L1I(MI),L2I(M2I),LCI(MCI),LBI(MBI)
      REAL FR(MF),F1I(MI),F2I(M2I),FCI(MCI),FBI(MBI)
      REAL GR(MG),G1I(MGI),G2I(MGI)
      REAL RLAT(MI),RLON(MI),DUM1(MI),DUM2(MI),AWT(MI)
      CHARACTER*1 CONE1,CONE2
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET MAP FIELD
      IF(JB.EQ.4) THEN
        KRB=0
        JGDS=-1
        CALL GETGRIB(LGB,LXB,MBE,MB,KRB,JPDSB,JGDS,
     &               MBUF,CBUFB,NLENB,NNUMB,MNUMB,
     &               KB,KRB,KPDSB,KGDSB,LR,FR,IRETB)
        IF(IRETB.EQ.0) THEN
          IF(IP.EQ.4) THEN
            INTB=1
          ELSE
            INTB=0
            DO I=1,11
              INTB=MAX(INTB,ABS(KGDSI(I)-KGDSB(I)))
            ENDDO
            INTB=MIN(INTB,1)
          ENDIF
          IBB=MOD(KPDSB(4)/64,2)
          IF(INTB.EQ.0) THEN
            KBI=KB
            IBBI=IBB
            DO I=1,KB
              LBI(I)=LR(I)
              FBI(I)=FR(I)
            ENDDO
          ELSE
           stop 3
!           CALL IPOLATES(IP,IPOPT,KGDSB,KGDSI,MB,MI,1,IBB,LR,FR,
!    &                    KBI,RLAT,RLON,IBBI,LBI,FBI,IRETB)
          ENDIF
        ENDIF
        IF(IRETB.NE.0) THEN
          CALL ERRMSG('diffgb: error preparing bitmap')
          CALL EXIT(IRETB)
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PRINT OUT HEADERS
      IF(LH.GT.0) THEN
        IV1=MVDATE(ID1,IF1)
        IF(MNUM1.GE.0) THEN
          CONE1=' '
          NR1=MNUM1+NNUM1
        ELSE
          CONE1='>'
          NR1=(-1-MNUM1)+NNUM1
        ENDIF
        IF(LD.EQ.0) THEN
          IF(IE1.EQ.IF1) THEN
            PRINT '("GRIB FILE  : GRID",I4,8X,
     &              "IDATE ",I8.8,4X,"FHOUR ",I4.2,4X,
     &              "VDATE ",I8.8,8X,
     &              A1,I6," GRIB MESSAGES")',
     &       IG1,ID1,IF1,IV1,CONE1,NR1
          ELSE
            PRINT '("GRIB FILE  : GRID",I4,8X,
     &              "IDATE ",I8.8,4X,"F",I4.2,":",I4.2,4X,
     &              "VDATE ",I8.8,8X,
     &              A1,I6," GRIB MESSAGES")',
     &       IG1,ID1,IE1,IF1,IV1,CONE1,NR1
          ENDIF
        ELSE
          IF(MNUM2.GE.0) THEN
            CONE2=' '
            NR2=MNUM2+NNUM2
          ELSE
            CONE2='>'
            NR2=(-1-MNUM2)+NNUM2
          ENDIF
          IV2=MVDATE(ID2,IF2)
          IF(IE1.EQ.IF1) THEN
            PRINT '("GRIB FILE 1: GRID",I4,8X,
     &              "IDATE ",I8.8,4X,"FHOUR ",I4.2,4X,
     &              "VDATE ",I8.8,8X,
     &              A1,I6," GRIB MESSAGES")',
     &     IG1,ID1,IF1,IV1,CONE1,NR1
          ELSE
            PRINT '("GRIB FILE 1: GRID",I4,8X,
     &              "IDATE ",I8.8,4X,"F",I4.2,":",I4.2,4X,
     &              "VDATE ",I8.8,8X,
     &              A1,I6," GRIB MESSAGES")',
     &       IG1,ID1,IE1,IF1,IV1,CONE1,NR1
          ENDIF
          IF(IE2.EQ.IF2) THEN
            PRINT '("GRIB FILE 2: GRID",I4,8X,
     &              "IDATE ",I8.8,4X,"FHOUR ",I4.2,4X,
     &              "VDATE ",I8.8,8X,
     &              A1,I6," GRIB MESSAGES")',
     &     IG2,ID2,IF2,IV2,CONE2,NR2
          ELSE
            PRINT '("GRIB FILE 2: GRID",I4,8X,
     &              "IDATE ",I8.8,4X,"F",I4.2,":",I4.2,4X,
     &              "VDATE ",I8.8,8X,
     &              A1,I6," GRIB MESSAGES")',
     &       IG2,ID2,IE2,IF2,IV2,CONE2,NR2
          ENDIF
        ENDIF
        IF(IGC.GT.0) THEN
          PRINT '("CLIMATOLOGY  GRID",I4)',IGC
        ENDIF
        JX=KGDSI(2)
        JY=KGDSI(3)
        IF(KGDSI(1).EQ.201.OR.KGDSI(1).EQ.202) THEN
          JX=KGDSI(7)
          JY=KGDSI(8)
        ENDIF
        IF(JB.EQ.0) THEN
          PRINT '("VERIFICATION GRID",I4,4X,
     &            " ( ",I4," X ",I4," ) ",4X,I7," POINTS")',
     &     IGI,JX,JY,MI
        ELSE
          PRINT '("VERIFICATION GRID",I4,4X,
     &            " ( ",I4," X ",I4," ) ",4X,I7," POINTS",4X,
     &            "(BITMAPPED)")',
     &     IGI,JX,JY,MI
        ENDIF
        IF(LV.EQ.1) THEN
          IF(LT.EQ.1) THEN
            CALL DGBHS1(LD)
          ELSEIF(LT.EQ.2) THEN
            CALL DGBHS2(LD,LW)
          ELSEIF(LT.EQ.3) THEN
            CALL DGBHS3(LD,LW,L3)
          ELSEIF(LT.EQ.4) THEN
            CALL DGBHS4(LD,LW)
          ELSEIF(LT.EQ.5) THEN
            CALL DGBHS5(LD,LW)
          ENDIF
        ELSE
          IF(LT.EQ.1) THEN
            CALL DGBHV1(LD)
          ELSEIF(LT.EQ.2) THEN
            CALL DGBHV2(LD,LW)
          ELSEIF(LT.EQ.3) THEN
            CALL DGBHV3(LD,LW,L3)
          ELSEIF(LT.EQ.4) THEN
            CALL DGBHV4(LD,LW)
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE AREA WEIGHTS
      IF(LT.GT.1.AND.LW.EQ.1) THEN
!       CALL GDSAWT(KGDSI,MI,KI,DUM1,DUM2,RLAT,RLON,AWT)
        stop 4
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET FIRST FIELD FROM FILE 1
      NTR=0
      NTC=0
      NTD=0
      KR1=0
      JGDS=-1
      CALL GETGRIB(LG1,LX1,M1E,M1,KR1,JPDS1,JGDS,
     &             MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &             K1,KR1,KPDS1,KGDS1,LR,FR,IRET1)
      IF(IRET1.EQ.0.AND.LV.EQ.2) THEN
        KRV=0
        JPDS=-1
        JPDS(1:21)=KPDS1(1:21)
        JPDS(5)=KPDS1(5)+1
        JGDS=KGDS1
        CALL GETGRIB(LG1,LX1,M1E,M1,KRV,JPDS,JGDS,
     &               MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &               K1,KRV,KPDS1,KGDS1,LR,GR,IRET1)
        KPDS1(5)=JPDS(5)-1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  LOOP UNTIL DONE
      DOWHILE(IRET1.EQ.0)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE FIELD 1
        NTR=NTR+1
        IRET=0
        IF(IP.EQ.4) THEN
          INT1=1
        ELSE
          INT1=0
          DO I=1,11
            INT1=MAX(INT1,ABS(KGDSI(I)-KGDS1(I)))
          ENDDO
          INT1=MIN(INT1,1)
        ENDIF
        IB1=MOD(KPDS1(4)/64,2)
        IF(INT1.EQ.0) THEN
          K1I=K1
          IB1I=IB1
          DO I=1,K1
            L1I(I)=LR(I)
            F1I(I)=FR(I)
            IF(LV.EQ.2) G1I(I)=GR(I)
          ENDDO
        ELSE
          IF(LV.NE.2) THEN
!            CALL IPOLATES(IP,IPOPT,KGDS1,KGDSI,M1,MI,1,IB1,LR,FR,
!     &                    K1I,RLAT,RLON,IB1I,L1I,F1I,IRET)
           stop 5
          ELSE
!            CALL IPOLATEV(IP,IPOPT,KGDS1,KGDSI,M1,MI,1,IB1,LR,FR,GR,
!     &                    K1I,RLAT,RLON,DUM1,DUM2,IB1I,L1I,F1I,G1I,IRET)
           stop 6
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FILTER FIELD 1
        IF(IRET.EQ.0.AND.KZ1.LE.KZ2) THEN
          IRET=IB1I
          IF(IRET.EQ.0) THEN
            CALL ZFF1(F1I)
            IF(LV.EQ.2) CALL ZFF1(G1I)
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET CLIMATOLOGY
        IF(IRET.EQ.0.AND.LV.NE.2.AND.IGC.GT.0) THEN
          KRC=0
          CALL GETCLIM(LGC,LXC,MCE,MC,KRC,MBUF,CBUFC,NLENC,NNUMC,MNUMC,
     &                 KPDS1,KPDSC,KGDSC,KC,LR,FR,IRET)
          IF(IRET.EQ.0) THEN
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE CLIMATOLOGY
            IF(IP.EQ.4) THEN
              INTC=1
            ELSE
              INTC=0
              DO I=1,11
                INTC=MAX(INTC,ABS(KGDSI(I)-KGDSC(I)))
              ENDDO
              INTC=MIN(INTC,1)
            ENDIF
            IBC=MOD(KPDSC(4)/64,2)
            IF(INTC.EQ.0) THEN
              KCI=KC
              IBCI=IBC
              DO I=1,KC
                LCI(I)=LR(I)
                FCI(I)=FR(I)
              ENDDO
            ELSE
!              CALL IPOLATES(IP,IPOPT,KGDSC,KGDSI,MC,MI,1,IBC,LR,FR,
!     &                      KCI,RLAT,RLON,IBCI,LCI,FCI,IRET)
              stop 7
            ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FILTER CLIMATOLOGY
            IF(IRET.EQ.0.AND.KZ1.LE.KZ2) THEN
              IRET=IBCI
              IF(IRET.EQ.0) THEN
                CALL ZFF1(FCI)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET FIELD FROM FILE 2
        IF(IRET.EQ.0.AND.LD.GT.0) THEN
          KR2=0
          JPDS=-1
          JPDS(5)=KPDS1(5)
          JPDS(6)=KPDS1(6)
          JPDS(7)=KPDS1(7)
cggg  next three lines for degrib of same date as file number 1
          jpds(8)=kpds1(8)
          jpds(9)=kpds1(9)
          jpds(10)=kpds1(10)
          jpds(13)=kpds1(13)
          jpds(14)=kpds1(14)
          jpds(15)=kpds1(15)
          jpds(16)=kpds1(16)
          jpds(19)=kpds1(19)
          IF(KPDS1(19).GT.128) JPDS(19)=KPDS1(19)
          JGDS=-1
          CALL GETGRIB(LG2,LX2,M2E,M2,KR2,JPDS,JGDS,
     &                 MBUF,CBUF2,NLEN2,NNUM2,MNUM2,
     &                 K2,KR2,KPDS2,KGDS2,LR,FR,IRET)
          IF(IRET.EQ.0.AND.LV.EQ.2) THEN
            KRV=0
            JPDS=-1
            JPDS(1:21)=KPDS2(1:21)
            JPDS(5)=KPDS2(5)+1
            JGDS=KGDS2
            CALL GETGRIB(LG2,LX2,M2E,M2,KRV,JPDS,JGDS,
     &                   MBUF,CBUF2,NLEN2,NNUM2,MNUM2,
     &                   K2,KRV,KPDS2,KGDS2,LR,GR,IRET)
            KPDS2(5)=JPDS(5)-1
          ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INTERPOLATE FIELD 2
          IF(IRET.EQ.0) THEN
            IF(IP.EQ.4) THEN
              INT2=1
            ELSE
              INT2=0
              DO I=1,11
                INT2=MAX(INT2,ABS(KGDSI(I)-KGDS2(I)))
              ENDDO
              INT2=MIN(INT2,1)
            ENDIF
            IB2=MOD(KPDS2(4)/64,2)
            IF(INT2.EQ.0) THEN
              K2I=K2
              IB2I=IB2
              DO I=1,K2
                L2I(I)=LR(I)
                F2I(I)=FR(I)
                IF(LV.EQ.2) G2I(I)=GR(I)
              ENDDO
            ELSE
              IF(LV.NE.2) THEN
!                CALL IPOLATES(IP,IPOPT,KGDS2,KGDSI,M2,MI,1,IB2,LR,FR,
!     &                        K2I,RLAT,RLON,IB2I,L2I,F2I,IRET)
           stop 8
              ELSE
!                CALL IPOLATEV(IP,IPOPT,KGDS2,KGDSI,M2,MI,1,IB2,LR,FR,GR,
!     &                        K2I,RLAT,RLON,DUM1,DUM2,
!     &                        IB2I,L2I,F2I,G2I,IRET)
           stop 9
              ENDIF
            ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  FILTER FIELD 2
            IF(IRET.EQ.0.AND.KZ1.LE.KZ2) THEN
              IRET=IB2I
              IF(IRET.EQ.0) THEN
                CALL ZFF1(F2I)
                IF(LV.EQ.2) CALL ZFF1(G1I)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  INVOKE MAP MASK
        IF(IRET.EQ.0) THEN
          NTC=NTC+1
          KI=K1I
          ID=0
          IF(JB.EQ.1) THEN
            DO I=1,KI
              IF(L1I(I)) THEN
                IF((LAB.EQ.1.AND.F1I(I).LE.AB).OR.
     &             (LAB.EQ.-1.AND.F1I(I).GE.AB)) THEN
                  L1I(I)=.FALSE.
                  IF(LD.GT.0) L2I(I)=.FALSE.
                  IF(IGC.GT.0) LCI(I)=.FALSE.
                ENDIF
              ELSE
                L1I(I)=.FALSE.
                IF(LD.GT.0) L2I(I)=.FALSE.
                IF(IGC.GT.0) LCI(I)=.FALSE.
              ENDIF
            ENDDO
          ELSEIF(JB.EQ.2) THEN
            DO I=1,KI
              IF(L2I(I)) THEN
                IF((LAB.EQ.1.AND.F2I(I).LE.AB).OR.
     &             (LAB.EQ.-1.AND.F2I(I).GE.AB)) THEN
                  L1I(I)=.FALSE.
                  IF(LD.GT.0) L2I(I)=.FALSE.
                  IF(IGC.GT.0) LCI(I)=.FALSE.
                ENDIF
              ELSE
                L1I(I)=.FALSE.
                IF(LD.GT.0) L2I(I)=.FALSE.
                IF(IGC.GT.0) LCI(I)=.FALSE.
              ENDIF
            ENDDO
          ELSEIF(JB.EQ.3) THEN
            DO I=1,KI
              IF(LCI(I)) THEN
                IF((LAB.EQ.1.AND.FCI(I).LE.AB).OR.
     &             (LAB.EQ.-1.AND.FCI(I).GE.AB)) THEN
                  L1I(I)=.FALSE.
                  IF(LD.GT.0) L2I(I)=.FALSE.
                  IF(IGC.GT.0) LCI(I)=.FALSE.
                ENDIF
              ELSE
                L1I(I)=.FALSE.
                IF(LD.GT.0) L2I(I)=.FALSE.
                IF(IGC.GT.0) LCI(I)=.FALSE.
              ENDIF
            ENDDO
          ELSEIF(JB.EQ.4) THEN
            DO I=1,KI
              IF(LBI(I)) THEN
                IF((LAB.EQ.1.AND.FBI(I).LE.AB).OR.
     &             (LAB.EQ.-1.AND.FBI(I).GE.AB)) THEN
                  L1I(I)=.FALSE.
                  IF(LD.GT.0) L2I(I)=.FALSE.
                  IF(IGC.GT.0) LCI(I)=.FALSE.
                ENDIF
              ELSE
                L1I(I)=.FALSE.
                IF(LD.GT.0) L2I(I)=.FALSE.
                IF(IGC.GT.0) LCI(I)=.FALSE.
              ENDIF
            ENDDO
          ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  COMPUTE AND PRINT STATISTICS
          IF(LV.EQ.1) THEN
            IF(LT.EQ.1) THEN
              CALL DGBFS1(LD,KI,
     &                    KR1,KPDS1,L1I,F1I,
     &                    KR2,KPDS2,L2I,F2I,ID)
            ELSEIF(LT.EQ.2) THEN
              CALL DGBFS2(LD,KI,LW,AWT,
     &                    KR1,KPDS1,L1I,F1I,
     &                    KR2,KPDS2,L2I,F2I,ID)
            ELSEIF(LT.EQ.3) THEN
              CALL DGBFS3(LD,KI,LW,AWT,L3,N3,A3,
     &                    KR1,KPDS1,L1I,F1I,
     &                    KR2,KPDS2,L2I,F2I,ID)
            ELSEIF(LT.EQ.4) THEN
              CALL DGBFS4(LD,KI,LW,AWT,
     &                    KR1,KPDS1,L1I,F1I,
     &                    KR2,KPDS2,L2I,F2I,ID)
            ELSEIF(LT.EQ.5) THEN
              CALL DGBFS5(LD,KI,LW,AWT,LCI,FCI,
     &                    KR1,KPDS1,L1I,F1I,
     &                    KR2,KPDS2,L2I,F2I,ID)
            ENDIF
          ELSE
            IF(LT.EQ.1) THEN
              CALL DGBFV1(LD,KI,
     &                    KR1,KPDS1,L1I,F1I,G1I,
     &                    KR2,KPDS2,L2I,F2I,G2I,ID)
            ELSEIF(LT.EQ.2) THEN
              CALL DGBFV2(LD,KI,LW,AWT,
     &                    KR1,KPDS1,L1I,F1I,G1I,
     &                    KR2,KPDS2,L2I,F2I,G2I,ID)
            ELSEIF(LT.EQ.3) THEN
              CALL DGBFV3(LD,KI,LW,AWT,L3,N3,A3,
     &                    KR1,KPDS1,L1I,F1I,G1I,
     &                    KR2,KPDS2,L2I,F2I,G2I,ID)
            ELSEIF(LT.EQ.4) THEN
              CALL DGBFV4(LD,KI,LW,AWT,
     &                    KR1,KPDS1,L1I,F1I,G1I,
     &                    KR2,KPDS2,L2I,F2I,G2I,ID)
            ENDIF
          ENDIF
          NTD=NTD+ID
        ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET FIELD FROM FILE 1
        JGDS=-1
        CALL GETGRIB(LG1,LX1,M1E,M1,KR1,JPDS1,JGDS,
     &               MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &               K1,KR1,KPDS1,KGDS1,LR,FR,IRET1)
        IF(IRET1.EQ.0.AND.LV.EQ.2) THEN
          KRV=0
          JPDS=-1
          JPDS(1:21)=KPDS1(1:21)
          JPDS(5)=KPDS1(5)+1
          JGDS=KGDS1
          CALL GETGRIB(LG1,LX1,M1E,M1,KRV,JPDS,JGDS,
     &                 MBUF,CBUF1,NLEN1,NNUM1,MNUM1,
     &                 K1,KRV,KPDS1,KGDS1,LR,GR,IRET1)
          KPDS1(5)=JPDS(5)-1
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  PRINT SUMMARY
      IF(LH.GT.0) THEN
        IF(LD.GT.0) THEN
          PRINT '("SUMMARY: ",I8," RECORDS READ,",4X,
     &            I8," RECORDS COMPARED,",
     &            I8," RECORDS DIFFER.")',NTR,NTC,NTD
        ELSE
          PRINT '("SUMMARY: ",I8," RECORDS READ.")',NTR
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHS1(LD)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHS1      PRINT HEADERS FOR SCALAR TYPE 1
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR SCALAR TYPE 1.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHS1(LD)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.GT.0) THEN
cggg
c       PRINT '("    PARM LTYPE LEVEL",
c    &          T25,2(" REC#  VPTS#",9X,"MIN",9X,"MAX"),
c    &          T97,"  VPTS#  NDIF#  IDIF#",6X,"MAXDIF")'
        PRINT '("    PARM LTYPE LEVEL",
     &          T25,2(" REC#  VPTS#",9X,"MIN",9X,"MAX"),
     &          T97,"  VPTS#  NDIF#   IDIF#",6X,"MAXDIF")'
      ELSE
        PRINT '("    PARM LTYPE LEVEL",
     &          T25," REC#  VPTS#",9X,"MIN",9X,"MAX")'
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFS1(LD,KI,
     &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFS1      PRINT STATISTICS FOR SCALAR TYPE 1
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR SCALAR TYPE 1.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFS1(LD,KI,
C    &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) DATA VALUES FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) DATA VALUES FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),F2I(KI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV2=0
      NV3=0
      FMIN1=1.E30
      FMAX1=-1.E30
      FMIN2=1.E30
      FMAX2=-1.E30
      ND3=0
      IM3=0
      DM3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,KI
        IF(L1I(I)) THEN
          NV1=NV1+1
          FMIN1=MIN(FMIN1,F1I(I))
          FMAX1=MAX(FMAX1,F1I(I))
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FMIN1,FMAX1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          IF(L2I(I)) THEN
            NV2=NV2+1
            FMIN2=MIN(FMIN2,F2I(I))
            FMAX2=MAX(FMAX2,F2I(I))
          ENDIF
          IF(L1I(I).AND.L2I(I)) THEN
            NV3=NV3+1
            IF(F1I(I).NE.F2I(I)) THEN
              ND3=ND3+1
              IF(ABS(F1I(I)-F2I(I)).GT.DM3) THEN
                DM3=ABS(F1I(I)-F2I(I))
                IM3=I
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(IM3.GT.0) DM3=F1I(IM3)-F2I(IM3)
cggg
c        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),3I7,1PE12.4)',
c    &   KPDS1(5),KPDS1(6),KPDS1(7),
c    &   KR1,NV1,FMIN1,FMAX1,
c    &   KR2,NV2,FMIN2,FMAX2,
c    &   NV3,ND3,IM3,DM3
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),2I7,1X,I7,1PE12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FMIN1,FMAX1,
     &   KR2,NV2,FMIN2,FMAX2,
     &   NV3,ND3,IM3,DM3
        IF(ND3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHS2(LD,LW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHS2      PRINT HEADERS FOR SCALAR TYPE 2
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR SCALAR TYPE 2.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHS2(LD,LW)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",7X,"GMEAN",6X,"GSIGMA"),
     &            T97,"  VPTS#",9X,"GMEAN",6X,"GSIGMA")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",7X,"GMEAN",6X,"GSIGMA")'
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",7X,"WMEAN",6X,"WSIGMA"),
     &            T97,"  VPTS#",9X,"WMEAN",6X,"WSIGMA")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",7X,"WMEAN",6X,"WSIGMA")'
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFS2(LD,KI,LW,AWT,
     &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFS2      PRINT STATISTICS FOR SCALAR TYPE 2
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR SCALAR TYPE 2.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFS2(LD,KI,LW,AWT,
C    &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) DATA VALUES FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) DATA VALUES FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),F2I(KI)
      DOUBLE PRECISION WV1,WV2,WV3
      DOUBLE PRECISION FAVG1,FAVG2,FAVG3,FRMS1,FRMS2,FRMS3
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO=0.)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV2=0
      NV3=0
      WV1=0
      WV2=0
      WV3=0
      FAVG1=0.
      FRMS1=0.
      FAVG2=0.
      FRMS2=0.
      FAVG3=0.
      FRMS3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,KI
        W=1.
        IF(LW.EQ.1) W=AWT(I)
        IF(L1I(I)) THEN
          NV1=NV1+1
          WV1=WV1+W
          FAVG1=FAVG1+W*F1I(I)
          FRMS1=FRMS1+W*F1I(I)**2
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        IF(WV1.GT.0) THEN
          FAVG1=FAVG1/WV1
          FRMS1=SQRT(MAX(FRMS1/WV1-FAVG1**2,DZERO))
        ENDIF
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L2I(I)) THEN
            NV2=NV2+1
            WV2=WV2+W
            FAVG2=FAVG2+W*F2I(I)
            FRMS2=FRMS2+W*F2I(I)**2
          ENDIF
          IF(L1I(I).AND.L2I(I)) THEN
            NV3=NV3+1
            WV3=WV3+W
            FAVG3=FAVG3+W*(F1I(I)-F2I(I))
            FRMS3=FRMS3+W*(F1I(I)-F2I(I))**2
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV1.GT.0) THEN
          FAVG1=FAVG1/WV1
          FRMS1=SQRT(MAX(FRMS1/WV1-FAVG1**2,DZERO))
        ENDIF
        IF(WV2.GT.0) THEN
          FAVG2=FAVG2/WV2
          FRMS2=SQRT(MAX(FRMS2/WV2-FAVG2**2,DZERO))
        ENDIF
        IF(WV3.GT.0) THEN
          FAVG3=FAVG3/WV3
          FRMS3=SQRT(MAX(FRMS3/WV3-FAVG3**2,DZERO))
        ENDIF
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),I7,2X,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1,
     &   KR2,NV2,FAVG2,FRMS2,
     &   NV3,FAVG3,FRMS3
        IF(FRMS3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHS3(LD,LW,L3)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHS3      PRINT HEADERS FOR SCALAR TYPE 3
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR SCALAR TYPE 3.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHS3(LD,LW,L3)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C     L3           INTEGER FLAG FOR THRESHOLD INEQUALITY
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41,2(" REC#  VPTS#",3X,"GTOT%",3X,"GVER%"),
     &              T97,"  VPTS#",4X,"GTHREAT%",2X,"GEQTHRT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41,2(" REC#  VPTS#",3X,"GTOT%",3X,"GVER%"),
     &              T97,"  VPTS#",4X,"GTHREAT%",2X,"GEQTHRT%")'
          ENDIF
        ELSE
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41," REC#  VPTS#",3X,"GTOT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41," REC#  VPTS#",3X,"GTOT%")'
          ENDIF
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41,2(" REC#  VPTS#",3X,"WTOT%",3X,"WVER%"),
     &              T97,"  VPTS#",4X,"WTHREAT%",2X,"WEQTHRT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41,2(" REC#  VPTS#",3X,"WTOT%",3X,"WVER%"),
     &              T97,"  VPTS#",4X,"WTHREAT%",2X,"WEQTHRT%")'
          ENDIF
        ELSE
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41," REC#  VPTS#",3X,"WTOT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41," REC#  VPTS#",3X,"WTOT%")'
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFS3(LD,KI,LW,AWT,L3,N3,A3,
     &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFS3      PRINT STATISTICS FOR SCALAR TYPE 3
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR SCALAR TYPE 3.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFS3(LD,KI,LW,AWT,L3,N3,A3,
C    &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     L3           INTEGER FLAG FOR THRESHOLD INEQUALITY
C     N3           INTEGER NUMBER OF THRESHOLDS
C     A3           REAL (N3) THRESHOLDS
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) DATA VALUES FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) DATA VALUES FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      REAL A3(N3)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),F2I(KI)
      DOUBLE PRECISION WV0,WV1,WV2,WV3
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      DO I3=1,N3
        AT=A3(I3)
        NV1=0
        NV2=0
        NV3=0
        WV0=0
        WV1=0
        WV2=0
        WV3=0
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(LD.EQ.0) THEN
          DO I=1,KI
            W=1.
            IF(LW.EQ.1) W=AWT(I)
            IF(L1I(I)) THEN
              NV1=NV1+1
              WV0=WV0+W
              IF(L3.GT.0) THEN
                IF(F1I(I).GT.AT) WV1=WV1+W
              ELSE
                IF(F1I(I).LT.AT) WV1=WV1+W
              ENDIF
            ENDIF
          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          FTOT1=0
          IF(WV0.GT.0) THEN
            FTOT1=WV1/WV0
          ENDIF
          PRINT '(2X,3I6,4X,1PE14.4,2X,I5,I7,2PF8.2)',
     &     KPDS1(5),KPDS1(6),KPDS1(7),A3(I3),
     &     KR1,NV1,FTOT1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ELSE
          DO I=1,KI
            W=1.
            IF(LW.EQ.1) W=AWT(I)
            IF(L1I(I).AND.L2I(I)) THEN
              NV1=NV1+1
              NV2=NV2+1
              NV3=NV3+1
              WV0=WV0+W
              IF(L3.GT.0) THEN
                IF(F1I(I).GT.AT) WV1=WV1+W
                IF(F2I(I).GT.AT) WV2=WV2+W
                IF(F1I(I).GT.AT.AND.F2I(I).GT.AT) WV3=WV3+W
              ELSE
                IF(F1I(I).LT.AT) WV1=WV1+W
                IF(F2I(I).LT.AT) WV2=WV2+W
                IF(F1I(I).LT.AT.AND.F2I(I).LT.AT) WV3=WV3+W
              ENDIF
            ENDIF
          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          FTOT1=0
          FVER1=0
          FTOT2=0
          FVER2=0
          FTHREAT=0
          FEQTHRT=0
          IF(WV0.GT.0) THEN
            FTOT1=WV1/WV0
            FTOT2=WV2/WV0
            FTOT3=WV3/WV0
            IF(FTOT1.GT.0) FVER1=FTOT3/FTOT1
            IF(FTOT2.GT.0) FVER2=FTOT3/FTOT2
            IF(FTOT1+FTOT2-FTOT3.GT.0) THEN
              FTHREAT=FTOT3/(FTOT1+FTOT2-FTOT3)
              IF(FTHREAT.LT.0.999999) ID=1
            ENDIF
            FEXP=FTOT1*FTOT2
            IF(FTOT1+FTOT2-FTOT3-FEXP.GT.0) THEN
              FEQTHRT=(FTOT3-FEXP)/(FTOT1+FTOT2-FTOT3-FEXP)
            ENDIF
          ENDIF
          PRINT '(2X,3I6,4X,1PE14.4,2X,2(I5,I7,2P2F8.2),
     &            I7,2X,2P2F10.2)',
     &     KPDS1(5),KPDS1(6),KPDS1(7),A3(I3),
     &     KR1,NV1,FTOT1,FVER1,
     &     KR2,NV2,FTOT2,FVER2,
     &     NV3,FTHREAT,FEQTHRT
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHS4(LD,LW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHS4      PRINT HEADERS FOR SCALAR TYPE 4
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR SCALAR TYPE 4.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHS4(LD,LW)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",7X,"GMEAN",6X,"GSIGMA"),
     &            T97,"  VPTS#",6X,"GCORR%",7X,"GRMSDIF")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",7X,"GMEAN",6X,"GSIGMA")'
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",7X,"WMEAN",6X,"WSIGMA"),
     &            T97,"  VPTS#",6X,"WCORR%",7X,"WRMSDIF")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",7X,"WMEAN",6X,"WSIGMA")'
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFS4(LD,KI,LW,AWT,
     &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFS4      PRINT STATISTICS FOR SCALAR TYPE 4
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR SCALAR TYPE 4.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFS4(LD,KI,LW,AWT,
C    &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) DATA VALUES FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) DATA VALUES FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),F2I(KI)
      DOUBLE PRECISION WV1,WV3
      DOUBLE PRECISION FAVG1,FAVG2,FRMS1,FRMS2,FRMS3,FCOR3
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO=0.)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV3=0
      WV1=0
      WV3=0
      FAVG1=0.
      FRMS1=0.
      FAVG2=0.
      FRMS2=0.
      FCOR3=0.
      FRMS3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L1I(I)) THEN
            NV1=NV1+1
            WV1=WV1+W
            FAVG1=FAVG1+W*F1I(I)
            FRMS1=FRMS1+W*F1I(I)**2
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV1.GT.0) THEN
          FAVG1=FAVG1/WV1
          FRMS1=SQRT(MAX(FRMS1/WV1-FAVG1**2,DZERO))
        ENDIF
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L1I(I).AND.L2I(I)) THEN
            NV3=NV3+1
            WV3=WV3+W
            FAVG1=FAVG1+W*F1I(I)
            FRMS1=FRMS1+W*F1I(I)**2
            FAVG2=FAVG2+W*F2I(I)
            FRMS2=FRMS2+W*F2I(I)**2
            FCOR3=FCOR3+W*F1I(I)*F2I(I)
            FRMS3=FRMS3+W*(F1I(I)-F2I(I))**2
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV3.GT.0) THEN
          FAVG1=FAVG1/WV3
          FRMS1=SQRT(MAX(FRMS1/WV3-FAVG1**2,DZERO))
          FAVG2=FAVG2/WV3
          FRMS2=SQRT(MAX(FRMS2/WV3-FAVG2**2,DZERO))
          IF(FRMS1*FRMS2.GT.0) THEN
            FCOR3=(FCOR3/WV3-FAVG1*FAVG2)/(FRMS1*FRMS2)
          ELSE
            FCOR3=0.
          ENDIF
          FRMS3=SQRT(FRMS3/WV3)
        ENDIF
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),
     &          I7,2X,2PF10.2,1PE14.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV3,FAVG1,FRMS1,
     &   KR2,NV3,FAVG2,FRMS2,
     &   NV3,FCOR3,FRMS3
        IF(FRMS3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHS5(LD,LW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHS5      PRINT HEADERS FOR SCALAR TYPE 5
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR SCALAR TYPE 5.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHS5(LD,LW)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",6X,"GAMEAN",7X,"GARMS"),
     &            T97,"  VPTS#",5X,"GACORR%",7X,"GRMSDIF")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",6X,"GAMEAN",7X,"GARMS")'
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",6X,"WAMEAN",7X,"WARMS"),
     &            T97,"  VPTS#",5X,"WACORR%",7X,"WRMSDIF")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",6X,"WAMEAN",7X,"WARMS")'
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFS5(LD,KI,LW,AWT,LCI,FCI,
     &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFS5      PRINT STATISTICS FOR SCALAR TYPE 5
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR SCALAR TYPE 5.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFS5(LD,KI,LW,AWT,
C    &                  KR1,KPDS1,L1I,F1I,KR2,KPDS2,L2I,F2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) DATA VALUES FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) DATA VALUES FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 LCI(KI),L1I(KI),L2I(KI)
      REAL FCI(KI),F1I(KI),F2I(KI)
      DOUBLE PRECISION WV1,WV3
      DOUBLE PRECISION FAVG1,FAVG2,FRMS1,FRMS2,FRMS3,FCOR3
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO=0.)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV3=0
      WV1=0
      WV3=0
      FAVG1=0.
      FRMS1=0.
      FAVG2=0.
      FRMS2=0.
      FCOR3=0.
      FRMS3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L1I(I).AND.LCI(I)) THEN
            NV1=NV1+1
            WV1=WV1+W
            FAVG1=FAVG1+W*(F1I(I)-FCI(I))
            FRMS1=FRMS1+W*(F1I(I)-FCI(I))**2
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV1.GT.0) THEN
          FAVG1=FAVG1/WV1
          FRMS1=SQRT(FRMS1/WV1)
        ENDIF
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L1I(I).AND.L2I(I)) THEN
            NV3=NV3+1
            WV3=WV3+W
            FAVG1=FAVG1+W*(F1I(I)-FCI(I))
            FRMS1=FRMS1+W*(F1I(I)-FCI(I))**2
            FAVG2=FAVG2+W*(F2I(I)-FCI(I))
            FRMS2=FRMS2+W*(F2I(I)-FCI(I))**2
            FCOR3=FCOR3+W*(F1I(I)-FCI(I))*(F2I(I)-FCI(I))
            FRMS3=FRMS3+W*(F1I(I)-F2I(I))**2
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV3.GT.0) THEN
          FAVG1=FAVG1/WV3
          FRMS1=SQRT(FRMS1/WV3)
          FAVG2=FAVG2/WV3
          FRMS2=SQRT(FRMS2/WV3)
          IF(FRMS1*FRMS2.GT.0) THEN
            FCOR3=FCOR3/WV3/(FRMS1*FRMS2)
          ELSE
            FCOR3=0.
          ENDIF
          FRMS3=SQRT(FRMS3/WV3)
        ENDIF
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),
     &          I7,2X,2PF10.2,1PE14.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV3,FAVG1,FRMS1,
     &   KR2,NV3,FAVG2,FRMS2,
     &   NV3,FCOR3,FRMS3
        IF(FRMS3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHV1(LD)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHV1      PRINT HEADERS FOR VECTOR TYPE 1
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR VECTOR TYPE 1.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHV1(LD)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.GT.0) THEN
        PRINT '("    PARM LTYPE LEVEL",
     &          T25,2(" REC#  VPTS#",8X,"SMIN",8X,"SMAX"),
     &          T97,"  VPTS#  NDIF#  IDIF#",5X,"MAXVDIF")'
      ELSE
        PRINT '("    PARM LTYPE LEVEL",
     &          T25," REC#  VPTS#",8X,"SMIN",8X,"SMAX")'
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFV1(LD,KI,
     &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFV1      PRINT STATISTICS FOR VECTOR TYPE 1
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR VECTOR TYPE 1.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFV1(LD,KI,
C    &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) U-COMPONENT FOR FIELD 1
C     G1I          REAL (KI) V-COMPONENT FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) U-COMPONENT FOR FIELD 2
C     G2I          REAL (KI) V-COMPONENT FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),G1I(KI),F2I(KI),G2I(KI)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV2=0
      NV3=0
      FMIN1=1.E30
      FMAX1=0.
      FMIN2=1.E30
      FMAX2=0.
      ND3=0
      IM3=1
      DM3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,KI
        IF(L1I(I)) THEN
          NV1=NV1+1
          SPD1=F1I(I)**2+G1I(I)**2
          FMIN1=MIN(FMIN1,SPD1)
          FMAX1=MAX(FMAX1,SPD1)
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        FMIN1=SQRT(FMIN1)
        FMAX1=SQRT(FMAX1)
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FMIN1,FMAX1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          IF(L2I(I)) THEN
            NV2=NV2+1
            SPD2=F2I(I)**2+G2I(I)**2
            FMIN2=MIN(FMIN2,SPD2)
            FMAX2=MAX(FMAX2,SPD2)
          ENDIF
          IF(L1I(I).AND.L2I(I)) THEN
            NV3=NV3+1
            IF(F1I(I).NE.F2I(I).AND.G1I(I).NE.G2I(I)) THEN
              ND3=ND3+1
              SPD3=(F1I(I)-F2I(I))**2+(G1I(I)-G2I(I))**2
              IF(SPD3.GT.DM3) THEN
                DM3=SPD3
                IM3=I
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        FMIN1=SQRT(FMIN1)
        FMAX1=SQRT(FMAX1)
        FMIN2=SQRT(FMIN2)
        FMAX2=SQRT(FMAX2)
        DM3=SQRT(DM3)
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),
     &          3I7,1PE12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FMIN1,FMAX1,
     &   KR2,NV2,FMIN2,FMAX2,
     &   NV3,ND3,IM3,DM3
        IF(ND3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHV2(LD,LW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHV2      PRINT HEADERS FOR VECTOR TYPE 2
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR VECTOR TYPE 2.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHV2(LD,LW)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",6X,"GSMEAN",5X,"GSSIGMA"),
     &            T97,"  VPTS#",8X,"GSMEAN",7X,"GVRMS")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",6X,"GSMEAN",5X,"GSSIGMA")'
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",6X,"WSMEAN",5X,"WSSIGMA"),
     &            T97,"  VPTS#",8X,"WSMEAN",7X,"WVRMS")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",6X,"WSMEAN",5X,"WSSIGMA")'
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFV2(LD,KI,LW,AWT,
     &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFV2      PRINT STATISTICS FOR VECTOR TYPE 2
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR VECTOR TYPE 2.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFV2(LD,KI,LW,AWT,
C    &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) U-COMPONENT FOR FIELD 1
C     G1I          REAL (KI) V-COMPONENT FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) U-COMPONENT FOR FIELD 2
C     G2I          REAL (KI) V-COMPONENT FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),G1I(KI),F2I(KI),G2I(KI)
      DOUBLE PRECISION WV1,WV2,WV3
      DOUBLE PRECISION FAVG1,FAVG2,FAVG3,FRMS1,FRMS2,FRMS3
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO=0.)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV2=0
      NV3=0
      WV1=0
      WV2=0
      WV3=0
      FAVG1=0.
      FRMS1=0.
      FAVG2=0.
      FRMS2=0.
      FAVG3=0.
      FRMS3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO I=1,KI
        W=1.
        IF(LW.EQ.1) W=AWT(I)
        IF(L1I(I)) THEN
          NV1=NV1+1
          WV1=WV1+W
          SPD1=F1I(I)**2+G1I(I)**2
          FAVG1=FAVG1+W*SQRT(SPD1)
          FRMS1=FRMS1+W*SPD1
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        IF(WV1.GT.0) THEN
          FAVG1=FAVG1/WV1
          FRMS1=SQRT(MAX(FRMS1/WV1-FAVG1**2,DZERO))
        ENDIF
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L2I(I)) THEN
            NV2=NV2+1
            WV2=WV2+W
            SPD2=F2I(I)**2+G2I(I)**2
            FAVG2=FAVG2+W*SQRT(SPD2)
            FRMS2=FRMS2+W*SPD2
          ENDIF
          IF(L1I(I).AND.L2I(I)) THEN
            NV3=NV3+1
            WV3=WV3+W
            SPD1=F1I(I)**2+G1I(I)**2
            SPD2=F2I(I)**2+G2I(I)**2
            SPD3=(F1I(I)-F2I(I))**2+(G1I(I)-G2I(I))**2
            FAVG3=FAVG3+W*(SQRT(SPD1)-SQRT(SPD2))
            FRMS3=FRMS3+W*SPD3
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV1.GT.0) THEN
          FAVG1=FAVG1/WV1
          FRMS1=SQRT(MAX(FRMS1/WV1-FAVG1**2,DZERO))
        ENDIF
        IF(WV2.GT.0) THEN
          FAVG2=FAVG2/WV2
          FRMS2=SQRT(MAX(FRMS2/WV2-FAVG2**2,DZERO))
        ENDIF
        IF(WV3.GT.0) THEN
          FAVG3=FAVG3/WV3
          FRMS3=SQRT(FRMS3/WV3)
        ENDIF
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),
     &          I7,2X,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1,
     &   KR2,NV2,FAVG2,FRMS2,
     &   NV3,FAVG3,FRMS3
        IF(FRMS3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHV3(LD,LW,L3)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHV3      PRINT HEADERS FOR VECTOR TYPE 3
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR VECTOR TYPE 3.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHV3(LD,LW,L3)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C     L3           INTEGER FLAG FOR THRESHOLD INEQUALITY
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41,2(" REC#  VPTS#",2X,"GSTOT%",2X,"GSVER%"),
     &              T97,"  VPTS#",3X,"GSTHREAT%",1X,"GSEQTHRT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41,2(" REC#  VPTS#",2X,"GSTOT%",2X,"GSVER%"),
     &              T97,"  VPTS#",3X,"GSTHREAT%",1X,"GSEQTHRT%")'
          ENDIF
        ELSE
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41," REC#  VPTS#",2X,"GSTOT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41," REC#  VPTS#",2X,"GSTOT%")'
          ENDIF
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41,2(" REC#  VPTS#",2X,"WSTOT%",2X,"WSVER%"),
     &              T97,"  VPTS#",3X,"WSTHREAT%",1X,"WSEQTHRT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41,2(" REC#  VPTS#",2X,"WSTOT%",2X,"WSVER%"),
     &              T97,"  VPTS#",3X,"WSTHREAT%",1X,"WSEQTHRT%")'
          ENDIF
        ELSE
          IF(L3.GT.0) THEN
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"> THRESHOLD",
     &              T41," REC#  VPTS#",2X,"WSTOT%")'
          ELSE
            PRINT '("    PARM LTYPE LEVEL",T25,3X,"< THRESHOLD",
     &              T41," REC#  VPTS#",2X,"WSTOT%")'
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFV3(LD,KI,LW,AWT,L3,N3,A3,
     &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFV3      PRINT STATISTICS FOR VECTOR TYPE 3
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR VECTOR TYPE 3.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFV3(LD,KI,LW,AWT,L3,N3,A3,
C    &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     L3           INTEGER FLAG FOR THRESHOLD INEQUALITY
C     N3           INTEGER NUMBER OF THRESHOLDS
C     A3           REAL (N3) THRESHOLDS
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) U-COMPONENT FOR FIELD 1
C     G1I          REAL (KI) V-COMPONENT FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) U-COMPONENT FOR FIELD 2
C     G2I          REAL (KI) V-COMPONENT FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      REAL A3(N3)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),G1I(KI),F2I(KI),G2I(KI)
      DOUBLE PRECISION WV0,WV1,WV2,WV3
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      DO I3=1,N3
        AT=A3(I3)**2
        NV1=0
        NV2=0
        NV3=0
        WV0=0
        WV1=0
        WV2=0
        WV3=0
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(LD.EQ.0) THEN
          DO I=1,KI
            W=1.
            IF(LW.EQ.1) W=AWT(I)
            IF(L1I(I)) THEN
              NV1=NV1+1
              WV0=WV0+W
              SPD1=F1I(I)**2+G1I(I)**2
              IF(L3.GT.0) THEN
                IF(SPD1.GT.AT) WV1=WV1+W
              ELSE
                IF(SPD1.LT.AT) WV1=WV1+W
              ENDIF
            ENDIF
          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          FTOT1=0
          IF(WV0.GT.0) THEN
            FTOT1=WV1/WV0
          ENDIF
          PRINT '(2X,3I6,4X,1PE14.4,2X,
     &            I5,I7,2PF8.2)',
     &     KPDS1(5),KPDS1(6),KPDS1(7),A3(I3),
     &     KR1,NV1,FTOT1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ELSE
          DO I=1,KI
            W=1.
            IF(LW.EQ.1) W=AWT(I)
            IF(L1I(I).AND.L2I(I)) THEN
              NV1=NV1+1
              NV2=NV2+1
              NV3=NV3+1
              WV0=WV0+W
              SPD1=F1I(I)**2+G1I(I)**2
              SPD2=F2I(I)**2+G2I(I)**2
              IF(L3.GT.0) THEN
                IF(SPD1.GT.AT) WV1=WV1+W
                IF(SPD2.GT.AT) WV2=WV2+W
                IF(SPD1.GT.AT.AND.SPD2.GT.AT) WV3=WV3+W
              ELSE
                IF(SPD1.LT.AT) WV1=WV1+W
                IF(SPD2.LT.AT) WV2=WV2+W
                IF(SPD1.LT.AT.AND.SPD2.LT.AT) WV3=WV3+W
              ENDIF
            ENDIF
          ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          FTOT1=0
          FVER1=0
          FTOT2=0
          FVER2=0
          FTHREAT=0
          FEQTHRT=0
          IF(WV0.GT.0) THEN
            FTOT1=WV1/WV0
            FTOT2=WV2/WV0
            FTOT3=WV3/WV0
            IF(FTOT1.GT.0) FVER1=FTOT3/FTOT1
            IF(FTOT2.GT.0) FVER2=FTOT3/FTOT2
            IF(FTOT1+FTOT2-FTOT3.GT.0) THEN
              FTHREAT=FTOT3/(FTOT1+FTOT2-FTOT3)
              IF(FTHREAT.LT.0.999999) ID=1
            ENDIF
            FEXP=FTOT1*FTOT2
            IF(FTOT1+FTOT2-FTOT3-FEXP.GT.0) THEN
              FEQTHRT=(FTOT3-FEXP)/(FTOT1+FTOT2-FTOT3-FEXP)
            ENDIF
          ENDIF
          PRINT '(2X,3I6,4X,1PE14.4,2X,
     &            2(I5,I7,2P2F8.2),I7,2X,2P2F10.2)',
     &     KPDS1(5),KPDS1(6),KPDS1(7),A3(I3),
     &     KR1,NV1,FTOT1,FVER1,
     &     KR2,NV2,FTOT2,FVER2,
     &     NV3,FTHREAT,FEQTHRT
        ENDIF
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBHV4(LD,LW)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBHV4      PRINT HEADERS FOR VECTOR TYPE 4
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT HEADERS FOR VECTOR TYPE 4.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBHV4(LD,LW)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     LW           INTEGER FLAG FOR W OPTION
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LW.EQ.0) THEN
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",6X,"GSMEAN",5X,"GSSIGMA"),
     &            T97,"  VPTS#",5X,"GVCORR%",9X,"GVRMS")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",6X,"GSMEAN",5X,"GSSIGMA")'
        ENDIF
      ELSE
        IF(LD.GT.0) THEN
          PRINT '("    PARM LTYPE LEVEL",
     &            T25,2(" REC#  VPTS#",6X,"WSMEAN",5X,"WSSIGMA"),
     &            T97,"  VPTS#",5X,"WVCORR%",9X,"WVRMS")'
        ELSE
          PRINT '("    PARM LTYPE LEVEL",
     &            T25," REC#  VPTS#",6X,"WSMEAN",5X,"WSSIGMA")'
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE DGBFV4(LD,KI,LW,AWT,
     &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    DGBFV4      PRINT STATISTICS FOR VECTOR TYPE 4
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: PRINT STATISTICS FOR VECTOR TYPE 4.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL DGBFV4(LD,KI,LW,AWT,
C    &                  KR1,KPDS1,L1I,F1I,G1I,KR2,KPDS2,L2I,F2I,G2I,ID)
C   INPUT ARGUMENTS:
C     LD           INTEGER FLAG FOR D OPTION
C     KI           INTEGER NUMBER OF POINTS IN GRID
C     LW           INTEGER FLAG FOR W OPTION
C     AWT          REAL (KI) AREA WEIGHTS (M**2)
C     KR1          INTEGER RECORD NUMBER FOR FIELD 1
C     KPDS1        INTEGER (200) PDS PARM FOR FIELD 1
C     L1I          LOGICAL*1 (KI) BITMAP FOR FIELD 1
C     F1I          REAL (KI) U-COMPONENT FOR FIELD 1
C     G1I          REAL (KI) V-COMPONENT FOR FIELD 1
C     KR2          INTEGER RECORD NUMBER FOR FIELD 2
C     KPDS2        INTEGER (200) PDS PARM FOR FIELD 2
C     L2I          LOGICAL*1 (KI) BITMAP FOR FIELD 2
C     F2I          REAL (KI) U-COMPONENT FOR FIELD 2
C     G2I          REAL (KI) V-COMPONENT FOR FIELD 2
C   OUTPUT ARGUMENTS:
C     ID           INTEGER DIFF FLAG (0 IF IDENTICAL, 1 IF DIFFERENT)
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL AWT(KI)
      INTEGER KPDS1(200),KPDS2(200)
      LOGICAL*1 L1I(KI),L2I(KI)
      REAL F1I(KI),G1I(KI),F2I(KI),G2I(KI)
      DOUBLE PRECISION WV1,WV2,WV3
      DOUBLE PRECISION FAVG1,GAVG1,FAVG2,GAVG2,FCOR3,FRMS1,FRMS2,FRMS3
      DOUBLE PRECISION DZERO
      PARAMETER(DZERO=0.)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ID=0
      NV1=0
      NV2=0
      NV3=0
      WV1=0
      WV2=0
      WV3=0
      FAVG1=0.
      GAVG1=0.
      FRMS1=0.
      FAVG2=0.
      GAVG2=0.
      FRMS2=0.
      FCOR3=0.
      FRMS3=0.
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LD.EQ.0) THEN
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L1I(I)) THEN
            NV1=NV1+1
            WV1=WV1+W
            FAVG1=FAVG1+W*F1I(I)
            GAVG1=GAVG1+W*G1I(I)
            FRMS1=FRMS1+W*(F1I(I)**2+G1I(I)**2)
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV1.GT.0) THEN
          FAVG1=SQRT(FAVG1**2+GAVG1**2)/WV1
          FRMS1=SQRT(MAX(FRMS1/WV1-FAVG1**2,DZERO))
        ENDIF
        PRINT '(2X,3I6,4X,I5,I7,1P2E12.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        DO I=1,KI
          W=1.
          IF(LW.EQ.1) W=AWT(I)
          IF(L1I(I).AND.L2I(I)) THEN
            NV1=NV1+1
            NV2=NV2+1
            NV3=NV3+1
            WV3=WV3+W
            FAVG1=FAVG1+W*F1I(I)
            GAVG1=GAVG1+W*G1I(I)
            FRMS1=FRMS1+W*(F1I(I)**2+G1I(I)**2)
            FAVG2=FAVG2+W*F2I(I)
            GAVG2=GAVG2+W*G2I(I)
            FRMS2=FRMS2+W*(F2I(I)**2+G2I(I)**2)
            FCOR3=FCOR3+W*(F1I(I)*F2I(I)+G1I(I)*G2I(I))
            FRMS3=FRMS3+W*((F1I(I)-F2I(I))**2+(G1I(I)-G2I(I))**2)
          ENDIF
        ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF(WV3.GT.0) THEN
          FAVG1=FAVG1/WV3
          GAVG1=GAVG1/WV3
          FRMS1=SQRT(MAX(FRMS1/WV3-(FAVG1**2+GAVG1**2),DZERO))
          FAVG2=FAVG2/WV3
          GAVG2=GAVG2/WV3
          FRMS2=SQRT(MAX(FRMS2/WV3-(FAVG2**2+GAVG2**2),DZERO))
          IF(FRMS1*FRMS2.GT.0) THEN
            FCOR3=(FCOR3/WV3-(FAVG1*FAVG2+GAVG1*GAVG2))/(FRMS1*FRMS2)
          ELSE
            FCOR3=0.
          ENDIF
          FAVG1=SQRT(FAVG1**2+GAVG1**2)
          FAVG2=SQRT(FAVG2**2+GAVG2**2)
          FRMS3=SQRT(FRMS3/WV3)
        ENDIF
        PRINT '(2X,3I6,4X,2(I5,I7,1P2E12.4),
     &          I7,2X,2PF10.2,1PE14.4)',
     &   KPDS1(5),KPDS1(6),KPDS1(7),
     &   KR1,NV1,FAVG1,FRMS1,
     &   KR2,NV2,FAVG2,FRMS2,
     &   NV3,FCOR3,FRMS3
        IF(FRMS3.GT.0) ID=1
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      FUNCTION LENGDSF(KGDS,KGDSF)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    LENGDSF     RETURN THE LENGTH OF A FILLED GRID
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: GIVEN A GRID DESCRIPTION SECTION (IN W3FI63 FORMAT),
C   RETURN THE GRID DESCRIPTION SECTION AND SIZE OF ITS REGULARIZED
C   COUNTERPART.  THAT IS, IF THE INPUT GRID IS REGULAR, THEN ITSELF
C   IS RETURNED ALONG WITH ITS GRID SIZE; HOWEVER IF THE INPUT GRID IS
C   ONLY QUASI-REGULAR (SUCH AS THE WAFS GRIDS), THEN ITS FILLED REGULAR
C   VERSION IS RETURNED ALONG WITH ITS FILLED GRID SIZE.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL LENGDSF(KGDS,KGDSF)
C   INPUT ARGUMENTS:
C     KGDS         INTEGER (200) GDS PARAMETERS IN W3FI63 FORMAT
C   OUTPUT ARGUMENTS:
C     KGDSF        INTEGER (200) REGULAR GDS PARAMETERS IN W3FI63 FORMAT
C     LENGDSF      INTEGER SIZE OF REGULARIZED GRID
C
C SUBPROGRAMS CALLED:
C   IPXWAFS
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER KGDS(200),KGDSF(200)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(KGDS(1).EQ.201) THEN
        KGDSF=KGDS
        LENGDSF=KGDS(7)*KGDS(8)-KGDS(8)/2
      ELSEIF(KGDS(1).EQ.202) THEN
        KGDSF=KGDS
        LENGDSF=KGDS(7)*KGDS(8)
      ELSEIF(KGDS(19).EQ.0.AND.KGDS(20).NE.255) THEN
!        CALL IPXWAFS(1,1,1,0,KGDS,DUM,KGDSF,DUMF,IRET)
         stop 10
        IF(IRET.EQ.0) THEN
          LENGDSF=KGDSF(2)*KGDSF(3)
        ELSE
          LENGDSF=0
        ENDIF
      ELSE
        KGDSF=KGDS
        LENGDSF=KGDS(2)*KGDS(3)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE GETGRIB(LG,LX,ME,M,JR,JPDS,JGDS,
     &                   MBUF,CBUF,NLEN,NNUM,MNUM,
     &                   K,KR,KPDS,KGDS,L,F,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GETGRIB     READ AND UNPACK GRIB MESSAGE
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: READ AND UNPACK GRIB MESSAGE.
C   SEE DOCBLOCK FOR GETGBM FOR FULLER DESCRIPTION.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL GETGRIB(LG,LX,ME,M,JR,JPDS,JGDS,
C    &                   MBUF,CBUF,NLEN,NNUM,MNUM,
C    &                   K,KR,KPDS,KGDS,L,F,IRET)
C   INPUT ARGUMENTS:
C     LG           INTEGER UNIT OF THE UNBLOCKED GRIB DATA FILE
C     LX           INTEGER UNIT OF THE UNBLOCKED GRIB INDEX FILE
C     ME           INTEGER ACTUAL GRID SIZE FOR IRREGULAR GRID
C     M            INTEGER MAXIMUM NUMBER OF DATA POINTS TO UNPACK
C     JR           INTEGER NUMBER OF MESSAGES TO SKIP
C     JPDS         INTEGER (200) PDS PARAMETERS FOR WHICH TO SEARCH
C     JGDS         INTEGER (200) GDS PARAMETERS FOR WHICH TO SEARCH
C     MBUF         INTEGER LENGTH OF INDEX BUFFER IN BYTES
C     CBUF         CHARACTER*1 (MBUF) INDEX BUFFER
C     NLEN         INTEGER LENGTH OF EACH INDEX RECORD IN BYTES
C     NNUM         INTEGER NUMBER OF INDEX RECORDS
C     MNUM         INTEGER NUMBER OF INDEX RECORDS SKIPPED
C   OUTPUT ARGUMENTS:
C     CBUF         CHARACTER*1 (MBUF) INDEX BUFFER
C     NLEN         INTEGER LENGTH OF EACH INDEX RECORD IN BYTES
C     NNUM         INTEGER NUMBER OF INDEX RECORDS
C     MNUM         INTEGER NUMBER OF INDEX RECORDS SKIPPED
C     K            INTEGER NUMBER OF DATA POINTS UNPACKED
C     KR           INTEGER MESSAGE NUMBER UNPACKED
C     KPDS         INTEGER (200) UNPACKED PDS PARAMETERS
C     KGDS         INTEGER (200) UNPACKED GDS PARAMETERS
C     L            LOGICAL*1 (M) UNPACKED BITMAP IF PRESENT
C     F            REAL (M) UNPACKED DATA
C     IRET         INTEGER RETURN CODE
C                  
C SUBPROGRAMS CALLED:
C   GETGBM 
C   IPXWAFS
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CHARACTER CBUF(MBUF)
      INTEGER KPDS(200),KGDS(200)
      LOGICAL*1 L(M)
      REAL F(M)
      INTEGER KGDSE(200)
      LOGICAL*1 LE(ME)
      REAL FE(ME)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(ME.EQ.M) THEN
        CALL GETGBM(LG,LX,M,JR,JPDS,JGDS,
     &              MBUF,CBUF,NLEN,NNUM,MNUM,
     &              K,KR,KPDS,KGDS,L,F,IRET)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ELSE
        CALL GETGBM(LG,LX,ME,JR,JPDS,JGDS,
     &              MBUF,CBUF,NLEN,NNUM,MNUM,
     &              K,KR,KPDS,KGDSE,LE,FE,IRET)
        IF(IRET.EQ.0) THEN
          IF(MOD(KPDS(4)/64,2).EQ.0.AND.KGDSE(20).NE.255) THEN
            K=LENGDSF(KGDSE,KGDS)
!            CALL IPXWAFS(1,ME,M,1,KGDSE,FE,KGDS,F,IRET)
            stop 11
            IF(IRET.EQ.0) L(1:M)=.TRUE.
          ELSE
            IRET=3
          ENDIF
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE GETCLIM(LG,LX,ME,M,KR,MBUF,CBUF,NLEN,NNUM,MNUM,
     &                   KPDS1,KPDS,KGDS,K,L,F,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    GETCLIM     GET CLIMATOLOGY
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: GET CLIMATOLOGICAL FIELD FROM MONTHLY CLIMATOLOGY.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL GETCLIM(LG,LX,ME,M,KR,MBUF,CBUF,NLEN,NNUM,MNUM,
C    &                   KPDS1,KPDS,KGDS,K,L,F,IRET)
C   INPUT ARGUMENTS:
C     LG           INTEGER UNIT OF THE UNBLOCKED GRIB DATA FILE
C     LX           INTEGER UNIT OF THE UNBLOCKED GRIB INDEX FILE
C     ME           INTEGER ACTUAL GRID SIZE FOR IRREGULAR GRID
C     M            INTEGER MAXIMUM NUMBER OF DATA POINTS TO UNPACK
C     KR           INTEGER NUMBER OF MESSAGES TO SKIP
C     MBUF         INTEGER LENGTH OF INDEX BUFFER IN BYTES
C     CBUF         CHARACTER*1 (MBUF) INDEX BUFFER
C     NLEN         INTEGER LENGTH OF EACH INDEX RECORD IN BYTES
C     NNUM         INTEGER NUMBER OF INDEX RECORDS
C     MNUM         INTEGER NUMBER OF INDEX RECORDS SKIPPED
C     KPDS1        INTEGER (200) UNPACKED PDS PARAMETERS
C   OUTPUT ARGUMENTS:
C     CBUF         CHARACTER*1 (MBUF) INDEX BUFFER
C     NLEN         INTEGER LENGTH OF EACH INDEX RECORD IN BYTES
C     NNUM         INTEGER NUMBER OF INDEX RECORDS
C     MNUM         INTEGER NUMBER OF INDEX RECORDS SKIPPED
C     KPDS         INTEGER (200) UNPACKED PDS PARAMETERS
C     KGDS         INTEGER (200) UNPACKED GDS PARAMETERS
C     K            INTEGER NUMBER OF DATA POINTS UNPACKED
C     L            LOGICAL*1 (M) UNPACKED BITMAP IF PRESENT
C     F            REAL (M) UNPACKED DATA
C     IRET         INTEGER RETURN CODE
C                  
C SUBPROGRAMS CALLED:
C   MVDATE
C   GETGRIB
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      CHARACTER CBUF(MBUF)
      INTEGER KPDS1(200),KPDS(200),KGDS(200)
      LOGICAL*1 L(M),L2(M)
      REAL F(M),F2(M)
      INTEGER JPDS(200),JGDS(200)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IDR=MOD(KPDS1(8),100)*1000000+KPDS1(9)*10000+
     &    KPDS1(10)*100+KPDS1(11)
      IFR=MAX(KPDS1(14),KPDS1(15))
      IVR=MVDATE(IDR,IFR)
      IMON=MOD(IVR/10000,100)
      IDAY=MOD(IVR/100,100)
      JPDS=-1
      JPDS(5)=KPDS1(5)
      JPDS(6)=KPDS1(6)
      JPDS(7)=KPDS1(7)
      JPDS(9)=IMON
      JPDS(16)=51
      KR=0
      CALL GETGRIB(LG,LX,ME,M,KR,JPDS,JGDS,
     &             MBUF,CBUF,NLEN,NNUM,MNUM,
     &             K,KR,KPDS,KGDS,L,F,IRET)
      IF(IRET.EQ.0.AND.IDAY.NE.15) THEN
        KR=0
        JPDS=KPDS
        JGDS=KGDS
        JPDS(9)=MOD(IMON+ISIGN(1,IDAY-15)+11,12)+1
        CALL GETGRIB(LG,LX,ME,M,KR,JPDS,JGDS,
     &               MBUF,CBUF,NLEN,NNUM,MNUM,
     &               K,KR,KPDS,KGDS,L2,F2,IRET)
        IF(IRET.EQ.0) THEN
          DO I=1,K
            IF(L(I).AND.L2(I)) F(I)=F(I)+ABS(IDAY-15)/30.*(F2(I)-F(I))
          ENDDO
        ENDIF
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE ZFF0(NW1,NW2,KGDS,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    ZFF0        SET UP ZONAL FOURIER FILTER
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: SET UP ZONAL FOURIER FILTER.  CALL ZFF1 TO INVOKE FILTER.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL ZFF0(NW1,NW2,KGDS,IRET)
C   INPUT ARGUMENTS:
C     NW1          INTEGER LOWEST ZONAL WAVENUMBER
C     NW2          INTEGER HIGHEST ZONAL WAVENUMBER
C     KGDS         INTEGER (200) UNPACKED GDS PARAMETERS
C   OUTPUT ARGUMENTS:
C     IRET         INTEGER RETURN CODE
C                  
C SUBPROGRAMS CALLED:
C   SPFFTE
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      INTEGER KGDS(200)
      PARAMETER(IGX=1024)
      COMMON/ZFFCOM/ KW1,KW2,IG,IM,JM,NSCAN,AFFT(100000+8*IGX)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IRET=0
      IDRTI=KGDS(1)
      IF(IDRTI.NE.0.AND.IDRTI.NE.1.AND.IDRTI.NE.4) IRET=1
      IF(IRET.EQ.0) THEN
        IM=KGDS(2)
        JM=KGDS(3)
        RLON1=KGDS(5)*1.E-3
        RLON2=KGDS(8)*1.E-3
        ISCAN=MOD(KGDS(11)/128,2)
        NSCAN=MOD(KGDS(11)/32,2)
        IF(ISCAN.EQ.0) THEN
          DLON=(MOD(RLON2-RLON1-1+3600,360.)+1)/(IM-1)
        ELSE
          DLON=-(MOD(RLON1-RLON2-1+3600,360.)+1)/(IM-1)
        ENDIF
        IG=NINT(360/ABS(DLON))
        IF(MOD(IG,2).NE.0.OR.IM.LT.IG) IRET=2
        IF(IG.GT.IGX) IRET=3
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IRET.EQ.0) THEN
        KW1=NW1
        KW2=NW2
       stop 12
c        CALL SPFFTE(IG,IG/2+1,IG,JM,W1,G1,0,AFFT)
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      SUBROUTINE ZFF1(F)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    ZFF1        INVOKE ZONAL FOURIER FILTER
C   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 96-07-19
C
C ABSTRACT: INVOKE ZONAL FOURIER FILTER.  CALL ZFF0 TO SET UP FILTER.
C
C PROGRAM HISTORY LOG:
C   96-07-19  IREDELL
C
C USAGE:    CALL ZFF1(F)
C   INPUT ARGUMENTS:
C     F            REAL (*) FIELD TO FOURIER FILTER
C                  
C SUBPROGRAMS CALLED:
C   SPFFTE
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN
C
C$$$
      REAL F(*)
      PARAMETER(IGX=1024)
      COMMON/ZFFCOM/ KW1,KW2,IG,IM,JM,NSCAN,AFFT(100000+8*IGX)
      REAL G1(IG,JM),W1(IG+2,JM)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NSCAN.EQ.0) THEN
        DO J=1,JM
          DO I=1,IG
            G1(I,J)=F(I+(J-1)*IM)
          ENDDO
        ENDDO
      ELSE
        DO J=1,JM
          DO I=1,IG
            G1(I,J)=F((I-1)*JM+J)
          ENDDO
        ENDDO
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO J=1,JM
        DO I=1,IG
          IF(NSCAN.EQ.0) THEN
            G1(I,J)=F(I+(J-1)*IM)
          ELSE
            G1(I,J)=F((I-1)*JM+J)
          ENDIF
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      stop 13
C     CALL SPFFTE(IG,IG/2+1,IG,JM,W1,G1,-1,AFFT)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      DO J=1,JM
        DO I=1,2*KW1
          W1(I,J)=0.
        ENDDO
        DO I=2*KW2+3,IG+2
          W1(I,J)=0.
        ENDDO
      ENDDO
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      stop 14
c     CALL SPFFTE(IG,IG/2+1,IG,JM,W1,G1,+1,AFFT)
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NSCAN.EQ.0) THEN
        DO J=1,JM
          DO I=1,IM
            IW=MOD(I-1,IG)+1
            F(I+(J-1)*IM)=G1(IW,J)
          ENDDO
        ENDDO
      ELSE
        DO J=1,JM
          DO I=1,IM
            IW=MOD(I-1,IG)+1
            F((I-1)*JM+J)=G1(IW,J)
          ENDDO
        ENDDO
      ENDIF
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
C-----------------------------------------------------------------------
      FUNCTION MVDATE(IYMDH,NH)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    MVDATE      RETURN A VERIFYING DATE IN YYMMDDHH FORMAT
C   PRGMMR: IREDELL          ORG: W/NMC2     DATE: 92-09-28
C
C ABSTRACT: THIS INTEGER FUNCTION COMPUTES THE VERIFYING DATE
C   GIVEN ITS INITIAL DATE AND HOUR AND ITS FORECAST LENGTH IN HOURS.
C   IF ONLY THE LAST TWO DIGITS ARE GIVEN FOR THE INITIAL YEAR,
C   THE YEAR IS ASSUMED TO BE FROM THE 100 YEARS ENDING ON MAXYR=2040.
C   THE VERIFYING DATE IS RETURNED AS AN INTEGER IN YYMMDDHH FORMAT.
C
C PROGRAM HISTORY LOG:
C   92-09-28  IREDELL
C   94-02-07  IREDELL    CHANGED ARGUMENT LIST
C
C USAGE:  ...=MVDATE(IYMDH,NH)
C   INPUT ARGUMENT LIST:
C     IYMDH    - INTEGER INITIAL DATE IN YYMMDDHH FORMAT
C     NH       - INTEGER FORECAST HOUR
C   OUTPUT ARGUMENT LIST:
C     MVDATE   - INTEGER VERIFYING DATE IN YYMMDDHH FORMAT
C
C$$$
      PARAMETER(MAXYR=2040)
      DIMENSION MT(13)
      DATA MT/365,396,59,90,120,151,181,212,243,273,304,334,365/
      MDY(KY)=KY*365+KY/4-KY/100+KY/400
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET VERIFYING REFERENCE DAY AND HOUR
      IY=IYMDH/1000000
      IM=(IYMDH-IY*1000000)/10000
      ID=(IYMDH-IY*1000000-IM*10000)/100
      IH=IYMDH-IY*1000000-IM*10000-ID*100
      JY=IY
      IF(IY.LT.100) JY=MAXYR-MOD(MAXYR-JY,100)
      IF(IM.LT.3) JY=JY-1
      KH=IH+NH
      IF(KH.GE.0) THEN
        KD=KH/24
      ELSE
        KD=(KH+1)/24-1
      ENDIF
      JH=KH-KD*24
      MDAY=MDY(JY)+MT(IM)+ID+KD
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C  GET VERIFYING DATE
      JY=MDAY/365.2425
      IF(MDAY-MDY(JY).LE.MT(3)) JY=JY-1
      MD=MDAY-MDY(JY)
      IF(MD.GT.MT(1)) THEN
        JY=JY+1
        JM=2
      ELSE
        JM=MD/30+1
      ENDIF
      IF(MD-MT(JM).LE.0) JM=JM-1
      JD=MD-MT(JM)
      MVDATE=MOD(JY,100)*1000000+JM*10000+JD*100+JH
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
C-----------------------------------------------------------------------
C  HP WORKSTATION
C  MACHINE DEPENDENT ROUTINES
C-----------------------------------------------------------------------
      SUBROUTINE FILENV
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FPARSEI(CARG,MARG,KARG)
      CHARACTER*(*) CARG
      INTEGER KARG(MARG)
      READ(CARG,*,IOSTAT=IOS) KARG
      END
C-----------------------------------------------------------------------
      SUBROUTINE FPARSER(CARG,MARG,RARG)
      CHARACTER*(*) CARG
      REAL RARG(MARG)
      READ(CARG,*,IOSTAT=IOS) RARG
      END
