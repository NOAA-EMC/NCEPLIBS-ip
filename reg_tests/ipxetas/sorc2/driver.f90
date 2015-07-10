 program driver

 use grib_mod

 implicit none

 character(len=100)                   :: input_file, output_file

 integer(kind=4)                      :: lugi, iunit, iret
 integer(kind=4)                      :: j, jdisc, jpdtn, jgdtn, k
 integer(kind=4)                      :: jids(200), jgdt(200), jpdt(200)
 integer                              :: idir, npts_input, npts_output, km
 integer                              :: igdtnum_output, istat
 integer(kind=4), allocatable, target :: igdtmpl_input(:), igdtmpl_output(:)

 logical                              :: unpack
 logical*1, allocatable, target       :: bitmap_input(:), bitmap_output(:)
 
 real, allocatable, target            :: data_input(:), data_output(:)

 type(gribfield)                      :: gfld_input, gfld_output

 input_file="./fort.9"
 print*,"- OPEN AND READ FILE ", trim(input_file)
 iunit=9
 call baopenr (iunit, input_file, iret)

 if (iret /= 0) then
   print*,'- BAD OPEN OF FILE, IRET IS ', iret
   stop 2
 end if

 j       = 0      ! search at beginning of file
 jdisc   = -1     ! search for any discipline
 jpdtn   = -1     ! search for any product definition template number
 jgdtn   =  1     ! search for grid definition template number 1 - rot lat/lon grid
 jids    = -9999  ! array of values in identification section, set to wildcard
 jgdt    = -9999  ! array of values in grid definition template 3.m
 jpdt    = -9999  ! array of values in product definition template 4.n
 unpack  = .true. ! unpack data
 lugi    = 0      ! no index file

 nullify(gfld_input%idsect)
 nullify(gfld_input%local)
 nullify(gfld_input%list_opt)
 nullify(gfld_input%igdtmpl)
 nullify(gfld_input%ipdtmpl)
 nullify(gfld_input%coord_list)
 nullify(gfld_input%idrtmpl)
 nullify(gfld_input%bmap)
 nullify(gfld_input%fld)

 call getgb2(iunit, lugi, j, jdisc, jids, jpdtn, jpdt, jgdtn, jgdt, &
             unpack, k, gfld_input, iret)

 if (iret /= 0) then
   print*,"- BAD DEGRIB OF DATA. IRET IS ", iret
   stop 4
 else
   print*,"- SUCCESSFULL DEGRIB OF DATA."
 end if

 call baclose (iunit, iret)

 print*,'- SCAN MODE IS: ', gfld_input%igdtmpl(19)
 print*,'- DIMENSION OF GRID: ', gfld_input%igdtmpl(8), gfld_input%igdtmpl(9)

 if (gfld_input%igdtmpl(19)==68.or.gfld_input%igdtmpl(19)==72)then
   idir=0
   npts_output = gfld_input%igdtmpl(9) * (gfld_input%igdtmpl(8)*2-1)
 else
   idir=-1   ! full to hpnts
!  idir=-2   ! full to vpnts
   npts_output = gfld_input%igdtmpl(9) * (gfld_input%igdtmpl(8)+1)/2
 endif

 npts_input = gfld_input%ngrdpts

 allocate(igdtmpl_input(gfld_input%igdtlen))
 igdtmpl_input=gfld_input%igdtmpl

 allocate(igdtmpl_output(gfld_input%igdtlen))
 igdtmpl_output=igdtmpl_input

 allocate(data_input(npts_input))
 data_input=gfld_input%fld

 allocate(bitmap_input(npts_input))
 bitmap_input=.true.

 allocate(data_output(npts_output))
 data_output=-999.

 allocate(bitmap_output(npts_output))
 bitmap_input=.true.

 call ipxetas(idir, gfld_input%igdtnum, gfld_input%igdtlen, igdtmpl_input, & 
              npts_input, bitmap_input, data_input, igdtnum_output,  &
              igdtmpl_output, npts_output, bitmap_output, data_output, istat)

 if (istat /= 0) then
   print*,"- PROBLEM IN IPXETAS. ISTAT IS: ", istat
   stop 5
 endif

 nullify(gfld_output%idsect)
 nullify(gfld_output%local)
 nullify(gfld_output%list_opt)
 nullify(gfld_output%igdtmpl)
 nullify(gfld_output%ipdtmpl)
 nullify(gfld_output%coord_list)
 nullify(gfld_output%idrtmpl)
 nullify(gfld_output%bmap)
 nullify(gfld_output%fld)

 gfld_output=gfld_input
 gfld_output%ngrdpts=npts_output
 gfld_output%igdtmpl=>igdtmpl_output
 gfld_output%idrtmpl=>gfld_input%idrtmpl
 gfld_output%idrtmpl=0
 gfld_output%idrtmpl(3)=gfld_input%idrtmpl(3)
 gfld_output%fld=>data_output

 if(any(.not.bitmap_output)) then
   gfld_output%ibmap=0
   gfld_output%bmap=>bitmap_output
 else
   gfld_output%ibmap=255
 endif

 output_file="./output.grb2"
 print*,"- OPEN AND WRITE FILE ", trim(output_file)
 iunit=10
 call baopenw(iunit, output_file, iret)

 call putgb2(iunit, gfld_output, iret)
 print*,'iret after write ', iret

 call baclose(iunit, iret)

 end program driver

 SUBROUTINE IPXETAS_NEW(IDIR, GDTNUMI, GDTLEN, GDTMPLI, NPTS_INPUT,  &
                        BITMAP_INPUT, DATA_INPUT, GDTNUMO, GDTMPLO, &
                        NPTS_OUTPUT, BITMAP_OUTPUT, DATA_OUTPUT, IRET)

 IMPLICIT NONE

 INTEGER, INTENT(IN   )            :: IDIR
 INTEGER, INTENT(IN   )            :: GDTNUMI, GDTLEN
 INTEGER(KIND=4), INTENT(IN   )    :: GDTMPLI(GDTLEN)
 INTEGER, INTENT(IN   )            :: NPTS_INPUT, NPTS_OUTPUT
 INTEGER, INTENT(  OUT)            :: GDTNUMO
 INTEGER(KIND=4), INTENT(  OUT)    :: GDTMPLO(GDTLEN)
 INTEGER, INTENT(  OUT)            :: IRET

 LOGICAL(KIND=1), INTENT(IN   )    :: BITMAP_INPUT(NPTS_INPUT)
 LOGICAL(KIND=1), INTENT(  OUT)    :: BITMAP_OUTPUT(NPTS_OUTPUT)

 REAL, INTENT(IN  )                :: DATA_INPUT(NPTS_INPUT)
 REAL, INTENT(OUT)                 :: DATA_OUTPUT(NPTS_OUTPUT)

 INTEGER                           :: SCAN_MODE, ISCALE, IP, IPOPT(20)
 INTEGER                           :: IBI(1), IBO(1), J, KM, NO

 REAL                              :: DLONS
 REAL, ALLOCATABLE                 :: OUTPUT_RLAT(:), OUTPUT_RLON(:)

 IRET = 0

 IF (GDTNUMI/=1) THEN
   IRET=1
   RETURN
 ENDIF

 SCAN_MODE=GDTMPLI(19)

 IF((SCAN_MODE==68.OR.SCAN_MODE==72).AND.(IDIR<-2.OR.IDIR>-1))THEN
   GDTNUMO=GDTNUMI
   GDTMPLO=GDTMPLI
   GDTMPLO(19)=64              
   GDTMPLO(8)=GDTMPLO(8)*2-1
   IF((GDTMPLO(8)*GDTMPLO(9))/=NPTS_OUTPUT)THEN
     IRET=3
     RETURN
   ENDIF
   ISCALE=GDTMPLO(10)*GDTMPLO(11)
   IF(ISCALE==0) ISCALE=1E6
   DLONS=FLOAT(GDTMPLO(17))/FLOAT(ISCALE)
   DLONS=DLONS*0.5
   GDTMPLO(17)=NINT(DLONS*FLOAT(ISCALE))
 ELSEIF(SCAN_MODE==64.AND.IDIR==-1)THEN  ! FULL TO H-GRID
   GDTNUMO=GDTNUMI
   GDTMPLO=GDTMPLI
   GDTMPLO(19)=68
   GDTMPLO(8)=(GDTMPLO(8)+1)/2
   IF((GDTMPLO(8)*GDTMPLO(9))/=NPTS_OUTPUT)THEN
     IRET=3
     RETURN
   ENDIF
   ISCALE=GDTMPLO(10)*GDTMPLO(11)
   IF(ISCALE==0) ISCALE=1E6
   DLONS=FLOAT(GDTMPLO(17))/FLOAT(ISCALE)
   DLONS=DLONS*2.0
   GDTMPLO(17)=NINT(DLONS*FLOAT(ISCALE))
 ELSEIF(SCAN_MODE==64.AND.IDIR==-2)THEN  ! FULL TO V-GRID
   GDTNUMO=GDTNUMI
   GDTMPLO=GDTMPLI
   GDTMPLO(19)=72
   GDTMPLO(8)=(GDTMPLO(8)+1)/2
   IF((GDTMPLO(8)*GDTMPLO(9))/=NPTS_OUTPUT)THEN
     IRET=3
     RETURN
   ENDIF
   ISCALE=GDTMPLO(10)*GDTMPLO(11)
   IF(ISCALE==0) ISCALE=1E6
   DLONS=FLOAT(GDTMPLO(17))/FLOAT(ISCALE)
   DLONS=DLONS*2.0
   GDTMPLO(17)=NINT(DLONS*FLOAT(ISCALE))
 ELSE
   IRET=2
   RETURN
 ENDIF

 KM=1
 IP=0
 IPOPT=0
 IBI=1  
 IBO=0

 ALLOCATE(OUTPUT_RLAT(NPTS_OUTPUT))
 ALLOCATE(OUTPUT_RLON(NPTS_OUTPUT))

 CALL IPOLATES(IP, IPOPT, GDTNUMI, GDTMPLI, GDTLEN, &
               GDTNUMO, GDTMPLO, GDTLEN, &
               NPTS_INPUT, NPTS_OUTPUT, KM, IBI, BITMAP_INPUT, DATA_INPUT, &
               NO, OUTPUT_RLAT, OUTPUT_RLON, IBO, BITMAP_OUTPUT, DATA_OUTPUT, IRET)

 DEALLOCATE(OUTPUT_RLAT, OUTPUT_RLON)

 IF(IRET /= 0)THEN
   PRINT*,'- PROBLEM IN IPOLATES: ', IRET
   RETURN
 ENDIF
 
 DO J=1, GDTMPLO(9)
   BITMAP_OUTPUT(J*GDTMPLO(8))=BITMAP_OUTPUT(J*GDTMPLO(8)-1)
   DATA_OUTPUT(J*GDTMPLO(8))=DATA_OUTPUT(J*GDTMPLO(8)-1)
   BITMAP_OUTPUT((J-1)*GDTMPLO(8)+1)=BITMAP_OUTPUT((J-1)*GDTMPLO(8)+2)
   DATA_OUTPUT((J-1)*GDTMPLO(8)+1)=DATA_OUTPUT((J-1)*GDTMPLO(8)+2)
 ENDDO

 RETURN

 END SUBROUTINE IPXETAS_NEW
