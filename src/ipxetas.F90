!> @file
!> @brief Expand or contract eta grids using linear interpolation.
!> @author Iredell @date 96-04-10

!> Expand or contract eta grids using linear interpolation.
!>
!> This subprogram transforms between the staggered eta grids as used
!> in the eta model and for native grid transmission and their full
!> expansion as used for general interpolation and graphics. The eta
!> grids are rotated latitude-longitude grids staggered as defined by
!> the arakawa e-grid, that is with mass data points alternating with
!> wind data points.
!>
!> This is similar to:
!> - ipxwafs2() which uses linear interpolation and accounts for
!> - bitmapped data.
!> - ipxwafs2() which uses neighbor interpolation and accounts for
!> - bitmapped data.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | Iredell | Initial
!> 2015-07-14 | Gayno | Make grib 2 compliant. Replace 4-pt interpolation with call to ipolates.
!> 2022-11-09 | Engle | Made ibi and ibo scalars.
!>
!> @param[in] idir integer transform option.
!> - 0 to expand staggered fields to full fields
!> - -1 to contract full mass fields to staggered fields
!> - -2 to contract full wind fields to staggered fields
!> @param[in] igdtnumi integer grid definition template number - input
!> grid. Corresponds to the gfld%igdtnum component of the ncep g2
!> library gridmod data structure. Must be = 1 (for a rotated lat/lon
!> grid.)
!> @param[in] igdtlen integer number of elements of the grid
!> definition template array - same for input and output grids (=22)
!> which are both rotated lat/lon grids.  corresponds to the
!> gfld%igdtlen component of the ncep g2 library gridmod data
!> structure.

!> @param[in] igdtmpli integer (igdtlen) grid definition template
!> array - input grid. Corresponds to the gfld%igdtmpl component of
!> the ncep g2 library gridmod data structure (section 3 info):

!> - 1 shape of earth, octet 15
!> - 2 scale factor of spherical earth radius, octet 16
!> - 3 scaled value of radius of spherical earth, octets 17-20
!> - 4 scale factor of major axis of elliptical earth, octet 21
!> - 5 scaled value of major axis of elliptical earth, octets 22-25
!> - 6 scale factor of minor axis of elliptical earth, octet 26
!> - 7 scaled value of minor axis of elliptical earth, octets 27-30
!> - 8 number of points along a parallel, octs 31-34
!> - 9 number of points along a meridian, octs 35-38
!> - 10 basic angle of initial production domain, octets 39-42
!> - 11 subdivisions of basic angle, octets 43-46
!> - 12 latitude of first grid point, octets 47-50
!> - 13 longitude of first grid point, octets 51-54
!> - 14 resolution and component flags, octet 55
!> - 15 latitude of last grid point, octets 56-59
!> - 16 longitude of last grid point, octets 60-63
!> - 17 i-direction increment, octets 64-67
!> - 18 j-direction increment, octets 68-71
!> - 19 scanning mode, octet 72
!> - 20 latitude of southern pole of projection, octets 73-76
!> - 21 longitude of southern pole of projection, octets 77-80
!> - 22 angle of rotation of projection, octs 81-84
!> @param[in] npts_input integer number points input grid
!> @param[in] bitmap_input logical (npts_input) input grid bitmap
!> @param[in] data_input real (npts_input) input grid data
!> @param[out] igdtnumo integer grid definition template number -
!> output grid. Corresponds to the gfld%igdtnum component of the ncep
!> g2 library gridmod data structure. Same as igdtnumi (=1 for a
!> rotated lat/lon grid).
!> @param[out] igdtmplo integer (igdtlen) grid definition template
!> array - output grid. Corresponds to the gfld%igdtmpl component of
!> the ncep g2 library gridmod data structure. Array definitions same
!> as "igdtmpli".
!> @param[out] npts_output integer number points output grid. the
!> j-dimension of the input and output grids are the same. When going
!> from a staggered to a full grid the i-dimension increases to
!> idim*2-1. When going from full to staggered the i-dimension
!> decreases to (idim+1)/2.
!>  @param[out] bitmap_output logical (npts_outut) output grid bitmap
!>  @param[out] data_output real (npts_output) output grid data
!> @param[out] iret integer return code
!> - 0 successful transformation
!> - non-0 invalid grid specs or problem in ipolates().
!>
! @author Iredell @date 96-04-10
 subroutine ipxetas(idir, igdtnumi, igdtlen, igdtmpli, npts_input, &
                    bitmap_input, data_input, igdtnumo, igdtmplo, &
                    npts_output, bitmap_output, data_output, iret)
     use ipolates_mod
     implicit none
!
     integer, intent(in)    :: idir
     integer, intent(in)    :: igdtnumi, igdtlen
     integer, intent(in)    :: igdtmpli(igdtlen)
     integer, intent(in)    :: npts_input, npts_output
     integer, intent(out)    :: igdtnumo
     integer, intent(out)    :: igdtmplo(igdtlen)
     integer, intent(out)    :: iret

     logical(KIND=1), intent(in)    :: bitmap_input(npts_input)
     logical(KIND=1), intent(out)    :: bitmap_output(npts_output)

     real, intent(in)    :: data_input(npts_input)
     real, intent(out)    :: data_output(npts_output)

     integer                           :: scan_mode, iscale, ip, ipopt(20)
     integer                           :: ibi, ibo, j, km, no

     real                              :: dlons
     real, allocatable                 :: output_rlat(:), output_rlon(:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     iret = 0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! ROUTINE ONLY WORKS FOR ROTATED LAT/LON GRIDS.
     if (igdtnumi .ne. 1) then
         iret = 1
         return
     end if
!
     scan_mode = igdtmpli(19)
     if ((scan_mode .eq. 68 .or. scan_mode .eq. 72) .and. (idir .lt. -2 .or. idir .gt. -1)) then
         igdtnumo = igdtnumi
         igdtmplo = igdtmpli
         igdtmplo(19) = 64
         igdtmplo(8) = igdtmplo(8)*2-1
         if ((igdtmplo(8)*igdtmplo(9)) .ne. npts_output) then
             iret = 3
             return
         end if
         iscale = igdtmplo(10)*igdtmplo(11)
         if (iscale .eq. 0) iscale = 10**6
         dlons = float(igdtmplo(17))/float(iscale)
         dlons = dlons*0.5
         igdtmplo(17) = nint(dlons*float(iscale))
     elseif (scan_mode .eq. 64 .and. idir .eq. -1) then  ! FULL TO H-GRID
         igdtnumo = igdtnumi
         igdtmplo = igdtmpli
         igdtmplo(19) = 68
         igdtmplo(8) = (igdtmplo(8)+1)/2
         if ((igdtmplo(8)*igdtmplo(9)) .ne. npts_output) then
             iret = 3
             return
         end if
         iscale = igdtmplo(10)*igdtmplo(11)
         if (iscale .eq. 0) iscale = 10**6
         dlons = float(igdtmplo(17))/float(iscale)
         dlons = dlons*2.0
         igdtmplo(17) = nint(dlons*float(iscale))
     elseif (scan_mode .eq. 64 .and. idir .eq. -2) then  ! FULL TO V-GRID
         igdtnumo = igdtnumi
         igdtmplo = igdtmpli
         igdtmplo(19) = 72
         igdtmplo(8) = (igdtmplo(8)+1)/2
         if ((igdtmplo(8)*igdtmplo(9)) .ne. npts_output) then
             iret = 3
             return
         end if
         iscale = igdtmplo(10)*igdtmplo(11)
         if (iscale .eq. 0) iscale = 10**6
         dlons = float(igdtmplo(17))/float(iscale)
         dlons = dlons*2.0
         igdtmplo(17) = nint(dlons*float(iscale))
     else
         iret = 2
         return
     end if

     km = 1
     ip = 0
     ipopt = 0
     ibi = 1
     ibo = 0

     allocate (output_rlat(npts_output))
     allocate (output_rlon(npts_output))

     call ipolates(ip, ipopt, igdtnumi, igdtmpli, igdtlen, &
                   igdtnumo, igdtmplo, igdtlen, &
                   npts_input, npts_output, km, ibi, bitmap_input, data_input, &
                   no, output_rlat, output_rlon, ibo, bitmap_output, data_output, iret)

     deallocate (output_rlat, output_rlon)

     if (iret .ne. 0) then
         print *, '- PROBLEM IN IPOLATES: ', iret
         return
     end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! REPLACE ANY UNDEFINED POINTS ALONG THE LEFT AND RIGHT EDGES.
     do j = 1, igdtmplo(9)
         bitmap_output(j*igdtmplo(8)) = bitmap_output(j*igdtmplo(8)-1)
         data_output(j*igdtmplo(8)) = data_output(j*igdtmplo(8)-1)
         bitmap_output((j-1)*igdtmplo(8)+1) = bitmap_output((j-1)*igdtmplo(8)+2)
         data_output((j-1)*igdtmplo(8)+1) = data_output((j-1)*igdtmplo(8)+2)
     end do

     return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end subroutine ipxetas
