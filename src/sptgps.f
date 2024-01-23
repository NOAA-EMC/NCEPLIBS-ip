C> @file
C> @brief Transform spectral scalar to polar stereo.
C>
C> ### Program History Log
C> Date | Programmer | Comments
C> -----|------------|---------
C> 96-02-29 | Iredell | Initial.
C> 1998-12-15 | Iredell | Openmp directives inserted.
C>
C> @author Iredell @date 96-02-29

C> This subprogram performs a spherical transform
C> from spectral coefficients of scalar quantities
C> to scalar fields on a pair of polar stereographic grids.
C>
C> The wave-space can be either triangular or rhomboidal.
C>
C> The wave and grid fields may have general indexing,
C> but each wave field is in sequential 'IBM order',
C> i.e. with zonal wavenumber as the slower index.
C>
C> The two square polar stereographic grids are centered
C> on the respective poles, with the orientation longitude
C> of the southern hemisphere grid 180 degrees opposite
C> that of the northern hemisphere grid.
C>
C> The transform is made efficient
C> by combining points in eight sectors
C> of each polar stereographic grid,   
C> numbered as in the diagram below.
C>
C> The pole and the sector boundaries  
C> are treated specially in the code.  
C>
C> Unfortunately, this approach induces
C> some hairy indexing and code loquacity.
C>
C> <pre>
C>              \ 4 | 5 /
C>               \  |  /
C>              3 \ | / 6
C>                 \|/
C>              ----+----
C>                 /|\
C>              2 / | \ 7
C>               /  |  \
C>              / 1 | 8 \
C> </pre>
C>
C> The transforms are all multiprocessed over sector points.
C>
C> Transform several fields at a time to improve vectorization.
C>
C> Subprogram can be called from a multiprocessing environment.
C>
C> @param IROMB spectral domain shape
C> (0 for triangular, 1 for rhomboidal)
C> @param MAXWV spectral truncation
C> @param KMAX number of fields to transform.
C> @param NPS odd order of the polar stereographic grids.
C> @param KWSKIP skip number between wave fields
C> (defaults to (MAXWV+1)*((IROMB+1)*MAXWV+2) if KWSKIP=0)
C> @param KGSKIP skip number between grid fields
C> (defaults to NPS*NPS if KGSKIP=0)
C> @param NISKIP skip number between grid i-points
C> (defaults to 1 if NISKIP=0)
C> @param NJSKIP skip number between grid j-points
C> (defaults to NPS if NJSKIP=0)
C> @param TRUE latitude at which ps grid is true (usually 60.)
C> @param XMESH grid length at true latitude (m)
C> @param ORIENT longitude at bottom of northern ps grid
C> (southern ps grid will have opposite orientation.)
C> @param WAVE wave fields
C> @param GN northern polar stereographic fields
C> @param GS southern polar stereographic fields
C>
C> @author Iredell @date 96-02-29
      SUBROUTINE SPTGPS(IROMB,MAXWV,KMAX,NPS,
     &                  KWSKIP,KGSKIP,NISKIP,NJSKIP,
     &                  TRUE,XMESH,ORIENT,WAVE,GN,GS)

      REAL WAVE(*),GN(*),GS(*)
      REAL EPS((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EPSTOP(MAXWV+1)
      REAL ENN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL ELONN1((MAXWV+1)*((IROMB+1)*MAXWV+2)/2)
      REAL EON((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),EONTOP(MAXWV+1)
      INTEGER MP(KMAX)
      REAL SLON(MAXWV,8),CLON(MAXWV,8),SROT(0:3),CROT(0:3)
      REAL WTOP(2*(MAXWV+1),KMAX)
      REAL PLN((MAXWV+1)*((IROMB+1)*MAXWV+2)/2),PLNTOP(MAXWV+1)
      REAL F(2*MAXWV+3,2,KMAX)
      DATA SROT/0.,1.,0.,-1./,CROT/1.,0.,-1.,0./
      PARAMETER(RERTH=6.3712E6)
      PARAMETER(PI=3.14159265358979,DPR=180./PI)

C  CALCULATE PRELIMINARY CONSTANTS
      CALL SPWGET(IROMB,MAXWV,EPS,EPSTOP,ENN1,ELONN1,EON,EONTOP)
      MX=(MAXWV+1)*((IROMB+1)*MAXWV+2)/2
      MXTOP=MAXWV+1
      IDIM=2*MAXWV+3
      KW=KWSKIP
      KG=KGSKIP
      NI=NISKIP
      NJ=NJSKIP
      IF(KW.EQ.0) KW=2*MX
      IF(KG.EQ.0) KG=NPS*NPS
      IF(NI.EQ.0) NI=1
      IF(NJ.EQ.0) NJ=NPS
      MP=0
      NPH=(NPS-1)/2
      GQ=((1.+SIN(TRUE/DPR))*RERTH/XMESH)**2
C$OMP PARALLEL DO
      DO K=1,KMAX
        WTOP(1:2*MXTOP,K)=0
      ENDDO

C  CALCULATE POLE POINT
      I1=NPH+1
      J1=NPH+1
      IJ1=(I1-1)*NI+(J1-1)*NJ+1
      SLAT1=1.
      CLAT1=0.
      CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &              PLN,PLNTOP)
      CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,KW,2*MXTOP,KMAX,
     &             CLAT1,PLN,PLNTOP,MP,WAVE,WTOP,F)
CDIR$ IVDEP
      DO K=1,KMAX
        IJK1=IJ1+(K-1)*KG
        GN(IJK1)=F(1,1,K)
        GS(IJK1)=F(1,2,K)
      ENDDO

C  CALCULATE POINTS ALONG THE ROW AND COLUMN OF THE POLE,
C  STARTING AT THE ORIENTATION LONGITUDE AND GOING CLOCKWISE.
C$OMP PARALLEL DO PRIVATE(I1,J2,I2,J3,I3,J4,I4,J5,I5,J6,I6,J7,I7,J8,I8)
C$OMP&       PRIVATE(IJ1,IJ2,IJ3,IJ4,IJ5,IJ6,IJ7,IJ8)
C$OMP&       PRIVATE(IJK1,IJK2,IJK3,IJK4,IJK5,IJK6,IJK7,IJK8)
C$OMP&       PRIVATE(DJ1,DI1,RQ,RADLON,RADLON1,RADLON2,SLAT1,CLAT1)
C$OMP&       PRIVATE(PLN,PLNTOP,F,SLON,CLON,LR,LI)
      DO J1=1,NPH
        I1=NPH+1
        RADLON=ORIENT/DPR
        J3=NPS+1-I1
        I3=J1
        J5=NPS+1-J1
        I5=NPS+1-I1
        J7=I1
        I7=NPS+1-J1
        IJ1=(I1-1)*NI+(J1-1)*NJ+1
        IJ3=(I3-1)*NI+(J3-1)*NJ+1
        IJ5=(I5-1)*NI+(J5-1)*NJ+1
        IJ7=(I7-1)*NI+(J7-1)*NJ+1
        DI1=I1-NPH-1
        DJ1=J1-NPH-1
        RQ=DI1**2+DJ1**2
        SLAT1=(GQ-RQ)/(GQ+RQ)
        CLAT1=SQRT(1.-SLAT1**2)
        CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,KW,2*MXTOP,KMAX,
     &               CLAT1,PLN,PLNTOP,MP,WAVE,WTOP,F)
        DO L=1,MAXWV
          SLON(L,1)=SIN(L*RADLON)
          CLON(L,1)=COS(L*RADLON)
          SLON(L,3)=SLON(L,1)*CROT(MOD(1*L,4))
     &             -CLON(L,1)*SROT(MOD(1*L,4))
          CLON(L,3)=CLON(L,1)*CROT(MOD(1*L,4))
     &             +SLON(L,1)*SROT(MOD(1*L,4))
          SLON(L,5)=SLON(L,1)*CROT(MOD(2*L,4))
     &             -CLON(L,1)*SROT(MOD(2*L,4))
          CLON(L,5)=CLON(L,1)*CROT(MOD(2*L,4))
     &             +SLON(L,1)*SROT(MOD(2*L,4))
          SLON(L,7)=SLON(L,1)*CROT(MOD(3*L,4))
     &             -CLON(L,1)*SROT(MOD(3*L,4))
          CLON(L,7)=CLON(L,1)*CROT(MOD(3*L,4))
     &             +SLON(L,1)*SROT(MOD(3*L,4))
        ENDDO
CDIR$ IVDEP
        DO K=1,KMAX
          IJK1=IJ1+(K-1)*KG
          IJK3=IJ3+(K-1)*KG
          IJK5=IJ5+(K-1)*KG
          IJK7=IJ7+(K-1)*KG
          GN(IJK1)=F(1,1,K)
          GN(IJK3)=F(1,1,K)
          GN(IJK5)=F(1,1,K)
          GN(IJK7)=F(1,1,K)
          GS(IJK1)=F(1,2,K)
          GS(IJK3)=F(1,2,K)
          GS(IJK5)=F(1,2,K)
          GS(IJK7)=F(1,2,K)
        ENDDO
        IF(KMAX.EQ.1) THEN
          DO L=1,MAXWV
            LR=2*L+1
            LI=2*L+2
            GN(IJ1)=GN(IJ1)+2*(F(LR,1,1)*CLON(L,1)
     &                        -F(LI,1,1)*SLON(L,1))
            GN(IJ3)=GN(IJ3)+2*(F(LR,1,1)*CLON(L,3)
     &                        -F(LI,1,1)*SLON(L,3))
            GN(IJ5)=GN(IJ5)+2*(F(LR,1,1)*CLON(L,5)
     &                        -F(LI,1,1)*SLON(L,5))
            GN(IJ7)=GN(IJ7)+2*(F(LR,1,1)*CLON(L,7)
     &                        -F(LI,1,1)*SLON(L,7))
            GS(IJ1)=GS(IJ1)+2*(F(LR,2,1)*CLON(L,5)
     &                        -F(LI,2,1)*SLON(L,5))
            GS(IJ3)=GS(IJ3)+2*(F(LR,2,1)*CLON(L,3)
     &                        -F(LI,2,1)*SLON(L,3))
            GS(IJ5)=GS(IJ5)+2*(F(LR,2,1)*CLON(L,1)
     &                        -F(LI,2,1)*SLON(L,1))
            GS(IJ7)=GS(IJ7)+2*(F(LR,2,1)*CLON(L,7)
     &                        -F(LI,2,1)*SLON(L,7))
          ENDDO
        ELSE
          DO L=1,MAXWV
            LR=2*L+1
            LI=2*L+2
CDIR$ IVDEP
            DO K=1,KMAX
              IJK1=IJ1+(K-1)*KG
              IJK3=IJ3+(K-1)*KG
              IJK5=IJ5+(K-1)*KG
              IJK7=IJ7+(K-1)*KG
              GN(IJK1)=GN(IJK1)+2*(F(LR,1,K)*CLON(L,1)
     &                            -F(LI,1,K)*SLON(L,1))
              GN(IJK3)=GN(IJK3)+2*(F(LR,1,K)*CLON(L,3)
     &                            -F(LI,1,K)*SLON(L,3))
              GN(IJK5)=GN(IJK5)+2*(F(LR,1,K)*CLON(L,5)
     &                            -F(LI,1,K)*SLON(L,5))
              GN(IJK7)=GN(IJK7)+2*(F(LR,1,K)*CLON(L,7)
     &                            -F(LI,1,K)*SLON(L,7))
              GS(IJK1)=GS(IJK1)+2*(F(LR,2,K)*CLON(L,5)
     &                            -F(LI,2,K)*SLON(L,5))
              GS(IJK3)=GS(IJK3)+2*(F(LR,2,K)*CLON(L,3)
     &                            -F(LI,2,K)*SLON(L,3))
              GS(IJK5)=GS(IJK5)+2*(F(LR,2,K)*CLON(L,1)
     &                            -F(LI,2,K)*SLON(L,1))
              GS(IJK7)=GS(IJK7)+2*(F(LR,2,K)*CLON(L,7)
     &                            -F(LI,2,K)*SLON(L,7))
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  CALCULATE POINTS ON THE MAIN DIAGONALS THROUGH THE POLE,
C  STARTING CLOCKWISE OF THE ORIENTATION LONGITUDE AND GOING CLOCKWISE.
C$OMP PARALLEL DO PRIVATE(I1,J2,I2,J3,I3,J4,I4,J5,I5,J6,I6,J7,I7,J8,I8)
C$OMP&       PRIVATE(IJ1,IJ2,IJ3,IJ4,IJ5,IJ6,IJ7,IJ8)
C$OMP&       PRIVATE(IJK1,IJK2,IJK3,IJK4,IJK5,IJK6,IJK7,IJK8)
C$OMP&       PRIVATE(DJ1,DI1,RQ,RADLON,RADLON1,RADLON2,SLAT1,CLAT1)
C$OMP&       PRIVATE(PLN,PLNTOP,F,SLON,CLON,LR,LI)
      DO J1=1,NPH
        I1=J1
        RADLON=(ORIENT-45)/DPR
        J3=NPS+1-I1
        I3=J1
        J5=NPS+1-J1
        I5=NPS+1-I1
        J7=I1
        I7=NPS+1-J1
        IJ1=(I1-1)*NI+(J1-1)*NJ+1
        IJ3=(I3-1)*NI+(J3-1)*NJ+1
        IJ5=(I5-1)*NI+(J5-1)*NJ+1
        IJ7=(I7-1)*NI+(J7-1)*NJ+1
        DI1=I1-NPH-1
        DJ1=J1-NPH-1
        RQ=DI1**2+DJ1**2
        SLAT1=(GQ-RQ)/(GQ+RQ)
        CLAT1=SQRT(1.-SLAT1**2)
        CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                PLN,PLNTOP)
        CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,KW,2*MXTOP,KMAX,
     &               CLAT1,PLN,PLNTOP,MP,WAVE,WTOP,F)
        DO L=1,MAXWV
          SLON(L,1)=SIN(L*RADLON)
          CLON(L,1)=COS(L*RADLON)
          SLON(L,3)=SLON(L,1)*CROT(MOD(1*L,4))
     &             -CLON(L,1)*SROT(MOD(1*L,4))
          CLON(L,3)=CLON(L,1)*CROT(MOD(1*L,4))
     &             +SLON(L,1)*SROT(MOD(1*L,4))
          SLON(L,5)=SLON(L,1)*CROT(MOD(2*L,4))
     &             -CLON(L,1)*SROT(MOD(2*L,4))
          CLON(L,5)=CLON(L,1)*CROT(MOD(2*L,4))
     &             +SLON(L,1)*SROT(MOD(2*L,4))
          SLON(L,7)=SLON(L,1)*CROT(MOD(3*L,4))
     &             -CLON(L,1)*SROT(MOD(3*L,4))
          CLON(L,7)=CLON(L,1)*CROT(MOD(3*L,4))
     &             +SLON(L,1)*SROT(MOD(3*L,4))
        ENDDO
CDIR$ IVDEP
        DO K=1,KMAX
          IJK1=IJ1+(K-1)*KG
          IJK3=IJ3+(K-1)*KG
          IJK5=IJ5+(K-1)*KG
          IJK7=IJ7+(K-1)*KG
          GN(IJK1)=F(1,1,K)
          GN(IJK3)=F(1,1,K)
          GN(IJK5)=F(1,1,K)
          GN(IJK7)=F(1,1,K)
          GS(IJK1)=F(1,2,K)
          GS(IJK3)=F(1,2,K)
          GS(IJK5)=F(1,2,K)
          GS(IJK7)=F(1,2,K)
        ENDDO
        IF(KMAX.EQ.1) THEN
          DO L=1,MAXWV
            LR=2*L+1
            LI=2*L+2
            GN(IJ1)=GN(IJ1)+2*(F(LR,1,1)*CLON(L,1)
     &                        -F(LI,1,1)*SLON(L,1))
            GN(IJ3)=GN(IJ3)+2*(F(LR,1,1)*CLON(L,3)
     &                        -F(LI,1,1)*SLON(L,3))
            GN(IJ5)=GN(IJ5)+2*(F(LR,1,1)*CLON(L,5)
     &                        -F(LI,1,1)*SLON(L,5))
            GN(IJ7)=GN(IJ7)+2*(F(LR,1,1)*CLON(L,7)
     &                        -F(LI,1,1)*SLON(L,7))
            GS(IJ1)=GS(IJ1)+2*(F(LR,2,1)*CLON(L,3)
     &                        -F(LI,2,1)*SLON(L,3))
            GS(IJ3)=GS(IJ3)+2*(F(LR,2,1)*CLON(L,1)
     &                        -F(LI,2,1)*SLON(L,1))
            GS(IJ5)=GS(IJ5)+2*(F(LR,2,1)*CLON(L,7)
     &                        -F(LI,2,1)*SLON(L,7))
            GS(IJ7)=GS(IJ7)+2*(F(LR,2,1)*CLON(L,5)
     &                        -F(LI,2,1)*SLON(L,5))
          ENDDO
        ELSE
          DO L=1,MAXWV
            LR=2*L+1
            LI=2*L+2
CDIR$ IVDEP
            DO K=1,KMAX
              IJK1=IJ1+(K-1)*KG
              IJK3=IJ3+(K-1)*KG
              IJK5=IJ5+(K-1)*KG
              IJK7=IJ7+(K-1)*KG
              GN(IJK1)=GN(IJK1)+2*(F(LR,1,K)*CLON(L,1)
     &                            -F(LI,1,K)*SLON(L,1))
              GN(IJK3)=GN(IJK3)+2*(F(LR,1,K)*CLON(L,3)
     &                            -F(LI,1,K)*SLON(L,3))
              GN(IJK5)=GN(IJK5)+2*(F(LR,1,K)*CLON(L,5)
     &                            -F(LI,1,K)*SLON(L,5))
              GN(IJK7)=GN(IJK7)+2*(F(LR,1,K)*CLON(L,7)
     &                            -F(LI,1,K)*SLON(L,7))
              GS(IJK1)=GS(IJK1)+2*(F(LR,2,K)*CLON(L,3)
     &                            -F(LI,2,K)*SLON(L,3))
              GS(IJK3)=GS(IJK3)+2*(F(LR,2,K)*CLON(L,1)
     &                            -F(LI,2,K)*SLON(L,1))
              GS(IJK5)=GS(IJK5)+2*(F(LR,2,K)*CLON(L,7)
     &                            -F(LI,2,K)*SLON(L,7))
              GS(IJK7)=GS(IJK7)+2*(F(LR,2,K)*CLON(L,5)
     &                            -F(LI,2,K)*SLON(L,5))
            ENDDO
          ENDDO
        ENDIF
      ENDDO

C  CALCULATE THE REMAINDER OF THE POLAR STEREOGRAPHIC DOMAIN,
C  STARTING AT THE SECTOR JUST CLOCKWISE OF THE ORIENTATION LONGITUDE
C  AND GOING CLOCKWISE UNTIL ALL EIGHT SECTORS ARE DONE.
C$OMP PARALLEL DO PRIVATE(I1,J2,I2,J3,I3,J4,I4,J5,I5,J6,I6,J7,I7,J8,I8)
C$OMP&       PRIVATE(IJ1,IJ2,IJ3,IJ4,IJ5,IJ6,IJ7,IJ8)
C$OMP&       PRIVATE(IJK1,IJK2,IJK3,IJK4,IJK5,IJK6,IJK7,IJK8)
C$OMP&       PRIVATE(DJ1,DI1,RQ,RADLON,RADLON1,RADLON2,SLAT1,CLAT1)
C$OMP&       PRIVATE(PLN,PLNTOP,F,SLON,CLON,LR,LI)
      DO J1=1,NPH-1
        DO I1=J1+1,NPH
          J2=I1
          I2=J1
          J3=NPS+1-I1
          I3=J1
          J4=NPS+1-J1
          I4=I1
          J5=NPS+1-J1
          I5=NPS+1-I1
          J6=NPS+1-I1
          I6=NPS+1-J1
          J7=I1
          I7=NPS+1-J1
          J8=J1
          I8=NPS+1-I1
          IJ1=(I1-1)*NI+(J1-1)*NJ+1
          IJ2=(I2-1)*NI+(J2-1)*NJ+1
          IJ3=(I3-1)*NI+(J3-1)*NJ+1
          IJ4=(I4-1)*NI+(J4-1)*NJ+1
          IJ5=(I5-1)*NI+(J5-1)*NJ+1
          IJ6=(I6-1)*NI+(J6-1)*NJ+1
          IJ7=(I7-1)*NI+(J7-1)*NJ+1
          IJ8=(I8-1)*NI+(J8-1)*NJ+1
          DI1=I1-NPH-1
          DJ1=J1-NPH-1
          RQ=DI1**2+DJ1**2
          SLAT1=(GQ-RQ)/(GQ+RQ)
          CLAT1=SQRT(1.-SLAT1**2)
          RADLON1=ORIENT/DPR+ATAN(-DI1/DJ1)
          RADLON2=(ORIENT-45)/DPR*2-RADLON1
          CALL SPLEGEND(IROMB,MAXWV,SLAT1,CLAT1,EPS,EPSTOP,
     &                  PLN,PLNTOP)
          CALL SPSYNTH(IROMB,MAXWV,2*MAXWV,IDIM,KW,2*MXTOP,KMAX,
     &                 CLAT1,PLN,PLNTOP,MP,WAVE,WTOP,F)
          DO L=1,MAXWV
            SLON(L,1)=SIN(L*RADLON1)
            CLON(L,1)=COS(L*RADLON1)
            SLON(L,2)=SIN(L*RADLON2)
            CLON(L,2)=COS(L*RADLON2)
            SLON(L,3)=SLON(L,1)*CROT(MOD(1*L,4))
     &               -CLON(L,1)*SROT(MOD(1*L,4))
            CLON(L,3)=CLON(L,1)*CROT(MOD(1*L,4))
     &               +SLON(L,1)*SROT(MOD(1*L,4))
            SLON(L,4)=SLON(L,2)*CROT(MOD(1*L,4))
     &               -CLON(L,2)*SROT(MOD(1*L,4))
            CLON(L,4)=CLON(L,2)*CROT(MOD(1*L,4))
     &               +SLON(L,2)*SROT(MOD(1*L,4))
            SLON(L,5)=SLON(L,1)*CROT(MOD(2*L,4))
     &               -CLON(L,1)*SROT(MOD(2*L,4))
            CLON(L,5)=CLON(L,1)*CROT(MOD(2*L,4))
     &               +SLON(L,1)*SROT(MOD(2*L,4))
            SLON(L,6)=SLON(L,2)*CROT(MOD(2*L,4))
     &               -CLON(L,2)*SROT(MOD(2*L,4))
            CLON(L,6)=CLON(L,2)*CROT(MOD(2*L,4))
     &               +SLON(L,2)*SROT(MOD(2*L,4))
            SLON(L,7)=SLON(L,1)*CROT(MOD(3*L,4))
     &               -CLON(L,1)*SROT(MOD(3*L,4))
            CLON(L,7)=CLON(L,1)*CROT(MOD(3*L,4))
     &               +SLON(L,1)*SROT(MOD(3*L,4))
            SLON(L,8)=SLON(L,2)*CROT(MOD(3*L,4))
     &               -CLON(L,2)*SROT(MOD(3*L,4))
            CLON(L,8)=CLON(L,2)*CROT(MOD(3*L,4))
     &               +SLON(L,2)*SROT(MOD(3*L,4))
          ENDDO
CDIR$ IVDEP
          DO K=1,KMAX
            IJK1=IJ1+(K-1)*KG
            IJK2=IJ2+(K-1)*KG
            IJK3=IJ3+(K-1)*KG
            IJK4=IJ4+(K-1)*KG
            IJK5=IJ5+(K-1)*KG
            IJK6=IJ6+(K-1)*KG
            IJK7=IJ7+(K-1)*KG
            IJK8=IJ8+(K-1)*KG
            GN(IJK1)=F(1,1,K)
            GN(IJK2)=F(1,1,K)
            GN(IJK3)=F(1,1,K)
            GN(IJK4)=F(1,1,K)
            GN(IJK5)=F(1,1,K)
            GN(IJK6)=F(1,1,K)
            GN(IJK7)=F(1,1,K)
            GN(IJK8)=F(1,1,K)
            GS(IJK1)=F(1,2,K)
            GS(IJK2)=F(1,2,K)
            GS(IJK3)=F(1,2,K)
            GS(IJK4)=F(1,2,K)
            GS(IJK5)=F(1,2,K)
            GS(IJK6)=F(1,2,K)
            GS(IJK7)=F(1,2,K)
            GS(IJK8)=F(1,2,K)
          ENDDO
          IF(KMAX.EQ.1) THEN
            DO L=1,MAXWV
              LR=2*L+1
              LI=2*L+2
              GN(IJ1)=GN(IJ1)+2*(F(LR,1,1)*CLON(L,1)
     &                          -F(LI,1,1)*SLON(L,1))
              GN(IJ2)=GN(IJ2)+2*(F(LR,1,1)*CLON(L,2)
     &                          -F(LI,1,1)*SLON(L,2))
              GN(IJ3)=GN(IJ3)+2*(F(LR,1,1)*CLON(L,3)
     &                          -F(LI,1,1)*SLON(L,3))
              GN(IJ4)=GN(IJ4)+2*(F(LR,1,1)*CLON(L,4)
     &                          -F(LI,1,1)*SLON(L,4))
              GN(IJ5)=GN(IJ5)+2*(F(LR,1,1)*CLON(L,5)
     &                          -F(LI,1,1)*SLON(L,5))
              GN(IJ6)=GN(IJ6)+2*(F(LR,1,1)*CLON(L,6)
     &                          -F(LI,1,1)*SLON(L,6))
              GN(IJ7)=GN(IJ7)+2*(F(LR,1,1)*CLON(L,7)
     &                          -F(LI,1,1)*SLON(L,7))
              GN(IJ8)=GN(IJ8)+2*(F(LR,1,1)*CLON(L,8)
     &                          -F(LI,1,1)*SLON(L,8))
              GS(IJ1)=GS(IJ1)+2*(F(LR,2,1)*CLON(L,4)
     &                          -F(LI,2,1)*SLON(L,4))
              GS(IJ2)=GS(IJ2)+2*(F(LR,2,1)*CLON(L,3)
     &                          -F(LI,2,1)*SLON(L,3))
              GS(IJ3)=GS(IJ3)+2*(F(LR,2,1)*CLON(L,2)
     &                          -F(LI,2,1)*SLON(L,2))
              GS(IJ4)=GS(IJ4)+2*(F(LR,2,1)*CLON(L,1)
     &                          -F(LI,2,1)*SLON(L,1))
              GS(IJ5)=GS(IJ5)+2*(F(LR,2,1)*CLON(L,8)
     &                          -F(LI,2,1)*SLON(L,8))
              GS(IJ6)=GS(IJ6)+2*(F(LR,2,1)*CLON(L,7)
     &                          -F(LI,2,1)*SLON(L,7))
              GS(IJ7)=GS(IJ7)+2*(F(LR,2,1)*CLON(L,6)
     &                          -F(LI,2,1)*SLON(L,6))
              GS(IJ8)=GS(IJ8)+2*(F(LR,2,1)*CLON(L,5)
     &                          -F(LI,2,1)*SLON(L,5))
            ENDDO
          ELSE
            DO L=1,MAXWV
              LR=2*L+1
              LI=2*L+2
CDIR$ IVDEP
              DO K=1,KMAX
                IJK1=IJ1+(K-1)*KG
                IJK2=IJ2+(K-1)*KG
                IJK3=IJ3+(K-1)*KG
                IJK4=IJ4+(K-1)*KG
                IJK5=IJ5+(K-1)*KG
                IJK6=IJ6+(K-1)*KG
                IJK7=IJ7+(K-1)*KG
                IJK8=IJ8+(K-1)*KG
                GN(IJK1)=GN(IJK1)+2*(F(LR,1,K)*CLON(L,1)
     &                              -F(LI,1,K)*SLON(L,1))
                GN(IJK2)=GN(IJK2)+2*(F(LR,1,K)*CLON(L,2)
     &                              -F(LI,1,K)*SLON(L,2))
                GN(IJK3)=GN(IJK3)+2*(F(LR,1,K)*CLON(L,3)
     &                              -F(LI,1,K)*SLON(L,3))
                GN(IJK4)=GN(IJK4)+2*(F(LR,1,K)*CLON(L,4)
     &                              -F(LI,1,K)*SLON(L,4))
                GN(IJK5)=GN(IJK5)+2*(F(LR,1,K)*CLON(L,5)
     &                              -F(LI,1,K)*SLON(L,5))
                GN(IJK6)=GN(IJK6)+2*(F(LR,1,K)*CLON(L,6)
     &                              -F(LI,1,K)*SLON(L,6))
                GN(IJK7)=GN(IJK7)+2*(F(LR,1,K)*CLON(L,7)
     &                              -F(LI,1,K)*SLON(L,7))
                GN(IJK8)=GN(IJK8)+2*(F(LR,1,K)*CLON(L,8)
     &                              -F(LI,1,K)*SLON(L,8))
                GS(IJK1)=GS(IJK1)+2*(F(LR,2,K)*CLON(L,4)
     &                              -F(LI,2,K)*SLON(L,4))
                GS(IJK2)=GS(IJK2)+2*(F(LR,2,K)*CLON(L,3)
     &                              -F(LI,2,K)*SLON(L,3))
                GS(IJK3)=GS(IJK3)+2*(F(LR,2,K)*CLON(L,2)
     &                              -F(LI,2,K)*SLON(L,2))
                GS(IJK4)=GS(IJK4)+2*(F(LR,2,K)*CLON(L,1)
     &                              -F(LI,2,K)*SLON(L,1))
                GS(IJK5)=GS(IJK5)+2*(F(LR,2,K)*CLON(L,8)
     &                              -F(LI,2,K)*SLON(L,8))
                GS(IJK6)=GS(IJK6)+2*(F(LR,2,K)*CLON(L,7)
     &                              -F(LI,2,K)*SLON(L,7))
                GS(IJK7)=GS(IJK7)+2*(F(LR,2,K)*CLON(L,6)
     &                              -F(LI,2,K)*SLON(L,6))
                GS(IJK8)=GS(IJK8)+2*(F(LR,2,K)*CLON(L,5)
     &                              -F(LI,2,K)*SLON(L,5))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO

      END
