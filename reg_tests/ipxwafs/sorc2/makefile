EXEC    =ipxwafs.exe
FC	=ifort
FOPTS	= -O0 -I/nwprod/lib/g2/v2.5.0/incmod/g2_v2.5.0_4
LIB1	=/nwprod/lib/libbacio_v2.0.1_4.a   
LIB2    =/nwprod/lib/g2/v2.5.0/libg2_v2.5.0_4.a
LIB3    =/nwprod/lib/libjasper.a
LIB4    =/nwprod/lib/libpng.a
LIB5    =/nwprod/lib/libz.a
LIB6    =/nwprod/lib/w3nco/v2.0.6/libw3nco_v2.0.6_4.a
LIB7    =/global/save/George.Gayno/iplib_branches/ip_grib2/ip/v2.0.4/libip_v2.0.4_4.a
                                                                                
$(EXEC) : ipxwafs.o
	$(FC) ${FOPTS} -o ${EXEC} ipxwafs.o ${LIB2} ${LIB1} $(LIB6) ${LIB3} $(LIB4) $(LIB5) $(LIB7)
                                                                                
ipxwafs.o: ipxwafs.f90
	$(FC) $(FOPTS) -c ipxwafs.f90

clean:
	rm -f *.o *.mod $(EXEC)
