! This is a test for the NCEPLIBS-ip library.
!
! Alex Richert, June 2023

#if (LSIZE==D)
#define REALSIZE 8
#elif (LSIZE==4)
#define REALSIZE 4
#endif
program test_ipxwafs

  implicit none

  integer          :: opt_pts(73)
  integer          :: igdtmpl_thin(19)=0
  integer          :: igdtmpl_full(19)
  integer          :: ib_thin(1)=0, ib_full(1)
  integer          :: iret, i, which_func
  integer, parameter :: nthin = 3447, nfull = 5329
  !
  logical(kind=1)  :: bitmap_thin(nthin,1)=.true.
  logical(kind=1)  :: bitmap_full(nfull,1)
  !
  real(KIND=REALSIZE)             :: data_thin(nthin,1)
  real(KIND=REALSIZE)             :: data_thin_contract(nthin,1)
  real(KIND=REALSIZE)             :: data_full(nfull,1)
  real             :: ref_data(10)
  real, parameter  :: abstol=1e-6

  opt_pts = (/ &
       73, 73, 73, 73, 73, 73, 73, 73, 72, 72, 72, 71, 71, 71, 70,&
       70, 69, 69, 68, 67, 67, 66, 65, 65, 64, 63, 62, 61, 60, 60,&
       59, 58, 57, 56, 55, 54, 52, 51, 50, 49, 48, 47, 45, 44, 43,&
       42, 40, 39, 38, 36, 35, 33, 32, 30, 29, 28, 26, 25, 23, 22,&
       20, 19, 17, 16, 14, 12, 11,  9,  8,  6,  5,  3,  2 /)

  do which_func=1,3
      ! Expand thin to full
      igdtmpl_thin(9) = 73
      igdtmpl_thin(18) = 1250000
      igdtmpl_thin(19) = 64

      do i=1,nthin
        data_thin(i,1) = real(i,KIND=REALSIZE)/nthin
      enddo

      if (which_func .eq. 1) then
          call ipxwafs(1, 3447, 5329, 1, 73, opt_pts, &
             19, igdtmpl_thin, data_thin, &
             igdtmpl_full, data_full, iret)
      elseif (which_func .eq. 2) then
          call ipxwafs2(1, 3447, 5329, 1, 73, opt_pts, &
             19, igdtmpl_thin, data_thin, ib_thin, bitmap_thin,  &
             igdtmpl_full, data_full, ib_full, bitmap_full, iret)
      elseif (which_func .eq. 3) then
          call ipxwafs3(1, 3447, 5329, 1, 73, opt_pts, &
             19, igdtmpl_thin, data_thin, ib_thin, bitmap_thin,  &
             igdtmpl_full, data_full, ib_full, bitmap_full, iret)
      endif

      if (which_func .ne. 3) then
          ref_data = (/ 0.68969796602520717, 0.68991151726138666, 0.69012506849756627, &
            0.69033861973374588, 0.69055217096992549, 0.69076572220610510, 0.69097927344228471, &
            0.69119282467846443, 0.69140637591464404, 0.69161992715082354 /)
      else
          ref_data = (/ 0.68958514650420655, 0.68987525384392223, 0.69016536118363792, &
            0.69045546852335360, 0.69045546852335360, 0.69074557586306928, 0.69103568320278508, &
            0.69132579054250076, 0.69132579054250076, 0.69161589788221645 /)
      endif
      if (.not. all(abs(data_full(2600:2609,1)-ref_data)<abstol)) stop 1

      ! Contract full to thin
      igdtmpl_full(8)=73
      igdtmpl_full(9)=73
      igdtmpl_full(17)=1250000
      if (which_func .eq. 1) then
          call ipxwafs(-1, 3447, 5329, 1, 73, opt_pts, &
             19, igdtmpl_thin, data_thin_contract, &
             igdtmpl_full, data_full, iret)
      elseif (which_func .eq. 2) then
          call ipxwafs2(-1, 3447, 5329, 1, 73, opt_pts, &
             19, igdtmpl_thin, data_thin_contract, ib_thin, bitmap_thin,  &
             igdtmpl_full, data_full, ib_full, bitmap_full, iret)
      elseif (which_func .eq. 3) then
          call ipxwafs3(-1, 3447, 5329, 1, 73, opt_pts, &
             19, igdtmpl_thin, data_thin_contract, ib_thin, bitmap_thin,  &
             igdtmpl_full, data_full, ib_full, bitmap_full, iret)
      endif
      if (.not. all(abs(data_thin-data_thin_contract)<abstol)) stop 2
  enddo ! which_func

end program test_ipxwafs
