!> @file
!> @brief Expand or contract wafs grids.
!> @author Iredell @date 96-04-10

!> Expand or contract wafs grids.
!>
!> This subprogram transforms between the thinned wafs grids used for
!> transmitting to the aviation community and their full expansion as
!> used for general interpolation and graphics.
!>
!> The thinned wafs grids are latitude-longitude grids where the
!> number of points in each row decrease toward the pole.
!>
!> This information is stored in the grib 2 grid definition template
!> (section 3) starting at octet 73. The full grid counterparts have
!> an equal number of points per row. The transform between the full
!> and thinned wafs grid is done by linear interpolation and is not
!> reversible.
!>
!> This routine does not work for bitmapped data.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 96-04-10 | iredell | initial version
!> 2015-july | gayno | convert to grib 2
!>
!> @param[in] idir integer transform option
!> - (+1 to expand thinned fields to full fields)
!> - (-1 to contract full fields to thinned fields)
!> @param[in] numpts_thin integer number of grid points - thinned
!> grid.  must be 3447.
!> @param[in] numpts_full integer number of grid points - full
!> grid. must be 5329 (73^2).
!> @param[in] km integer number of fields to transform
!> @param[in] num_opt integer number of values to describe the thinned
!> grid. must be 73.  dimension of array opt_pts.
!> @param[inout] opt_pts integer (num_opt) number of grid points per
!> row - thinned grid - if idir=+1
!> @param[in] igdtlen integer grid defintion template array length.
!> must be 19 for lat/lon grids. corresponds to the gfld%igdtlen
!> component of the ncep g2 library gridmod data structure.  same for
!> thin and full grids which are both lat/lon.
!> @param[inout] igdtmpl_thin integer (igdtlen) grid definition
!> template array - thinned grid - if idir=+1. corresponds to the
!> gfld%igdtmpl component of the ncep g2 library gridmod data
!> structure (section 3 info):
!> - 1 shape of earth, octet 15
!> - 2 scale factor of spherical earth radius, octet 16
!> - 3 scaled value of radius of spherical earth, octets 17-20
!> - 4 scale factor of major axis of elliptical earth, octet 21
!> - 5 scaled value of major axis of elliptical earth, octets 22-25
!> - 6 scale factor of minor axis of elliptical earth, octet 26
!> - 7 scaled value of minor axis of elliptical earth, octets 27-30
!> - 8 set to missing for thinned grid., octs 31-34
!> - 9 number of points along a meridian, octs 35-38
!> - 10 basic angle of initial production domain, octets 39-42.
!> - 11 subdivisions of basic angle, octets 43-46
!> - 12 latitude of first grid point, octets 47-50
!> - 13 longitude of first grid point, octets 51-54
!> - 14 resolution and component flags, octet 55
!> - 15 latitude of last grid point, octets 56-59
!> - 16 longitude of last grid point, octets 60-63
!> - 17 set to missing for thinned grid, octets 64-67
!> - 18 j-direction increment, octets 68-71
!> - 19 scanning mode, octet 72
!> @param[inout] data_thin real (numpts_thin,km) thinned grid fields if idir=+1

!> @param[inout] igdtmpl_full integer (igdtlen) grid definition
!> template array - full grid - if idir=-1. corresponds to the
!> gfld%igdtmpl component of the ncep g2 library gridmod data
!> structure.  same as igdtmpl_thin except:
!> - 8 number of points along a parallel, octs 31-34
!> - 17 i-direction increment, octets 64-67
!> @param[inout] data_full real (numpts_full,km) full grid fields if idir=-1
!> @param[out] iret integer return code
!> - 0 successful transformation
!> - 1 improper grid specification
!>
!> @author Iredell @date 96-04-10
 subroutine ipxwafs(idir, numpts_thin, numpts_full, km, num_opt, &
                    opt_pts, igdtlen, igdtmpl_thin, data_thin, &
                    igdtmpl_full, data_full, iret)
     implicit none
!
     integer, intent(in) :: num_opt
     integer, intent(inout) :: opt_pts(num_opt)
     integer, intent(in) :: idir, km, numpts_thin, numpts_full
     integer, intent(in) :: igdtlen
     integer, intent(inout) :: igdtmpl_thin(igdtlen)
     integer, intent(inout) :: igdtmpl_full(igdtlen)
     integer, intent(out) :: iret
!
     real, intent(inout) :: data_thin(numpts_thin, km)
     real, intent(inout) :: data_full(numpts_full, km)
!
     integer, parameter     :: missing = -1
!
     integer                              :: scan_mode, i, j, k, idlat, idlon
     integer                              :: ia, ib, im, im1, im2, npwafs(73)
     integer                              :: is1, is2, iscan, iscale
!
     logical                              :: test1, test2
!
     real                                 :: dlon, hi
     real                                 :: rat1, rat2, rlon1, rlon2
     real                                 :: wa, wb, x1, x2
!
     data npwafs/ &
         73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70, &
         70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60, &
         59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43, &
         42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22, &
         20, 19, 17, 16, 14, 12, 11, 9, 8, 6, 5, 3, 2/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  TRANSFORM GDS
     iret = 0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  REG LAT/LON GRIDS HAVE 19 GDT ELEMENTS.
     if (igdtlen .ne. 19 .or. numpts_thin .ne. 3447 .or. numpts_full .ne. 5329) then
         iret = 1
         return
     end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EXPAND THINNED GDS TO FULL GDS
     if (idir .gt. 0) then
         scan_mode = igdtmpl_thin(19)
         iscale = igdtmpl_thin(10)*igdtmpl_thin(11)
         if (iscale .eq. 0) iscale = 10**6
         idlat = nint(1.25*float(iscale))
         test1 = all(opt_pts .eq. npwafs)
         test2 = all(opt_pts .eq. npwafs(73:1:-1))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SOME CHECKS TO ENSURE THIS IS A WAFS GRID
         if (scan_mode .eq. 64 .and. igdtmpl_thin(9) .eq. 73 .and. &
             idlat .eq. igdtmpl_thin(18) .and. (test1 .or. test2)) then
             igdtmpl_full = igdtmpl_thin
             im = 73
             igdtmpl_full(8) = im
             rlon1 = float(igdtmpl_full(13))/float(iscale)
             rlon2 = float(igdtmpl_full(16))/float(iscale)
             iscan = mod(igdtmpl_full(19)/128, 2)
             hi = (-1.)**iscan
             dlon = hi*(mod(hi*(rlon2-rlon1)-1+3600, 360.)+1)/(im-1)
             igdtmpl_full(17) = nint(dlon*float(iscale))
         else
             iret = 1
         end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  CONTRACT FULL GDS TO THINNED GDS
     elseif (idir .lt. 0) then
         scan_mode = igdtmpl_full(19)
         iscale = igdtmpl_full(10)*igdtmpl_full(11)
         if (iscale .eq. 0) iscale = 10**6
         idlat = nint(1.25*float(iscale))
         idlon = idlat
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SOME CHECKS TO ENSURE THIS IS A WAFS GRID
         if (scan_mode .eq. 64 .and. igdtmpl_full(8) .eq. 73 .and. igdtmpl_full(9) .eq. 73 .and. &
             num_opt .eq. 73 .and. idlat .eq. igdtmpl_full(18) .and. idlon .eq. igdtmpl_full(17)) then
             igdtmpl_thin = igdtmpl_full
             igdtmpl_thin(8) = missing
             igdtmpl_thin(17) = missing
             if (igdtmpl_thin(12) .eq. 0) then  ! IS LATITUDE OF ROW 1 THE EQUATOR?
                 opt_pts = npwafs
             else
                 opt_pts = npwafs(73:1:-1)
             end if
         else
             iret = 1
         end if
     end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  TRANSFORM FIELDS
     if (iret .eq. 0) then
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  EXPAND THINNED FIELDS TO FULL FIELDS
         if (idir .eq. 1) then
             do k = 1, km
                 is1 = 0
                 is2 = 0
                 do j = 1, igdtmpl_full(9)
                     im1 = opt_pts(j)
                     im2 = igdtmpl_full(8)
                     rat1 = float(im1-1)/float(im2-1)
                     do i = 1, im2
                         x1 = (i-1)*rat1+1
                         ia = int(x1)
                         ia = min(max(ia, 1), im1-1)
                         ib = ia+1
                         wa = ib-x1
                         wb = x1-ia
                         data_full(is2+i, k) = wa*data_thin(is1+ia, k)+wb*data_thin(is1+ib, k)
                     end do
                     is1 = is1+im1
                     is2 = is2+im2
                 end do
             end do
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  CONTRACT FULL FIELDS TO THINNED FIELDS
         elseif (idir .eq. -1) then
             do k = 1, km
                 is1 = 0
                 is2 = 0
                 do j = 1, igdtmpl_full(9)
                     im1 = opt_pts(j)
                     im2 = igdtmpl_full(8)
                     rat2 = float(im2-1)/float(im1-1)
                     do i = 1, im1
                         x2 = (i-1)*rat2+1
                         ia = int(x2)
                         ia = min(max(ia, 1), im2-1)
                         ib = ia+1
                         wa = ib-x2
                         wb = x2-ia
                         data_thin(is1+i, k) = wa*data_full(is2+ia, k)+wb*data_full(is2+ib, k)
                     end do
                     is1 = is1+im1
                     is2 = is2+im2
                 end do
             end do
         end if
     end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 end subroutine ipxwafs
