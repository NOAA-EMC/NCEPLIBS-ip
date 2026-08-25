! Unit tests for spectral_interp_mod.F90
!
! Covers POLATES4 (scalar) and POLATEV4 (vector) for both GRIB1 and
! GRIB2 descriptors, including error branches and all fast-transform
! special cases (global cylindrical, polar stereographic, Mercator).
!
! Alex Richert, 2026
program test_spectral_interp_mod
  use ip_mod
  implicit none

  integer, parameter :: MISSING = huge(1)
  integer, parameter :: IP_SPEC = SPECTRAL_INTERP_ID    ! 4

  ! ipopt: triangular truncation, auto wave-number
  integer :: ipopt(20)

  ! Standard GRIB2 360x180 global equidistant cylindrical input grid.
  ! ISCAN=0, JSCAN=1, NSCAN=0 (scan mode 64). IDRTI=256.
  integer, parameter :: IN2_NUM = 0
  integer, parameter :: IN2_LEN = 19
  integer :: in2_tmpl(IN2_LEN)
  integer, parameter :: MI2 = 360*180   ! 64800 input points

  ! Standard GRIB1 360x180 global equidistant cylindrical input.
  integer :: kgdsi(200)
  integer, parameter :: MI1 = 360*180

  ! GRIB2 output grids
  integer, parameter :: OUT2_LL_LEN = 19
  integer :: out2_ll_tmpl(OUT2_LL_LEN)
  integer, parameter :: MO2_LL = 360*181

  integer, parameter :: OUT2_GAUSS_LEN = 19
  integer :: out2_gauss_tmpl(OUT2_GAUSS_LEN)
  integer, parameter :: MO2_GAUSS = 192*94

  integer, parameter :: OUT2_PS_LEN = 18
  integer :: out2_ps_tmpl(OUT2_PS_LEN)
  integer, parameter :: MO2_PS = 3*3

  integer, parameter :: OUT2_MERC_LEN = 19
  integer :: out2_merc_tmpl(OUT2_MERC_LEN)
  integer, parameter :: MO2_MERC = 116*44

  ! GRIB1 output grids
  integer :: out1_ll_kgds(200)
  integer :: out1_ps_kgds(200)
  integer :: out1_merc_kgds(200)

  integer, parameter :: MO1_LL   = 360*181
  integer, parameter :: MO1_PS   = 3*3
  integer, parameter :: MO1_MERC = 116*44

  ! Initialise arrays
  ipopt = [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

  in2_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
      360, 180, 0, MISSING, -89500000, -180000000, &
      48, 89500000, 179000000, 1000000, 1000000, 64]

  out2_ll_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
      360, 181, 0, MISSING, -90000000, 0, &
      48, 90000000, 359000000, 1000000, 1000000, 0]

  ! T62 Gaussian: N=47 -> JG=94=JMO; first Gaussian lat ~88.542 deg.
  out2_gauss_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
      192, 94, 0, MISSING, 88542000, 0, &
      48, -88542000, 358125000, 1875000, 47, 0]

  ! NH polar stereo 3x3: NPS=3, Dx=Dy=500 km, scan=64.
  ! res/comp flags=56 -> MOD(56/8,2)=1 (needed by polatev4).
  out2_ps_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
      3, 3, 83194000, 315000000, 56, 60000000, 0, 500000000, 500000000, 0, 64]

  out2_merc_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
      116, 44, -48670000, 3104000, 56, 22500000, 61050000, 0, 64, 0, 318830000, 318830000]

  kgdsi = 0
  kgdsi(1)=0;      kgdsi(2)=360;     kgdsi(3)=180
  kgdsi(4)=-89500; kgdsi(5)=-180000; kgdsi(6)=128
  kgdsi(7)=89500;  kgdsi(8)=179000;  kgdsi(9)=1000
  kgdsi(10)=1000;  kgdsi(11)=64

  out1_ll_kgds = 0
  out1_ll_kgds(1)=0;      out1_ll_kgds(2)=360;    out1_ll_kgds(3)=181
  out1_ll_kgds(4)=90000;  out1_ll_kgds(5)=0;      out1_ll_kgds(6)=128
  out1_ll_kgds(7)=-90000; out1_ll_kgds(8)=-1000;  out1_ll_kgds(9)=1000
  out1_ll_kgds(10)=1000;  out1_ll_kgds(11)=0

  ! NH polar stereo 3x3: kgds(6)=56 -> MOD(56/8,2)=1 (needed by polatev4).
  ! RLAT1=83.193, RLON1=315, ORIENT=0 => NINT(XP)=NINT(YP)=IP=2.
  out1_ps_kgds = 0
  out1_ps_kgds(1)=5;      out1_ps_kgds(2)=3;      out1_ps_kgds(3)=3
  out1_ps_kgds(4)=83193;  out1_ps_kgds(5)=315000; out1_ps_kgds(6)=56
  out1_ps_kgds(7)=0;      out1_ps_kgds(8)=500000; out1_ps_kgds(9)=500000
  out1_ps_kgds(10)=0;     out1_ps_kgds(11)=64

  out1_merc_kgds = 0
  out1_merc_kgds(1)=1;      out1_merc_kgds(2)=116;    out1_merc_kgds(3)=44
  out1_merc_kgds(4)=-48670; out1_merc_kgds(5)=3104;   out1_merc_kgds(6)=128
  out1_merc_kgds(7)=61050;  out1_merc_kgds(8)=0;      out1_merc_kgds(9)=22500
  out1_merc_kgds(11)=64;    out1_merc_kgds(13)=318830

  ! POLATES4_grib2
  call test_s2_iret41_bad_type()     ! line 331: bad input type -> IRET=41
  call test_s2_error_return_path()   ! lines 482-486: IBO=1, LO=F, GO=0
  call test_s2_iscan1()              ! line 340: ISCAN=1 negative DLON
  call test_s2_gaussian_jm_mismatch()! lines 357-359: Gaussian JM mismatch
  call test_s2_nscan1()              ! lines 383-384: NSCAN=1 path
  call test_s2_global_cyl_latlon_out() ! lines 390-411 IDRTO=0: SPTRUN
  call test_s2_global_cyl_gauss_out() ! lines 390-411 IDRTO=4: SPTRUN
  call test_s2_polar_stereo_out()    ! lines 417-469: SPTRUNS
  call test_s2_mercator_out()        ! lines 417-469: SPTRUNM

  ! POLATES4_grib1
  call test_s1_global_cyl_out()      ! line 681: SPTRUN
  call test_s1_polar_stereo_out()    ! line 706: SPTRUNS
  call test_s1_mercator_out()        ! line 734: SPTRUNM

  ! POLATEV4_grib2 (also exercises lines 148-155 via interpolate_spectral_vector)
  call test_v2_global_cyl_out()      ! line 1021: SPTRUNV
  call test_v2_polar_stereo_out()    ! line 1051: SPTRUNSV
  call test_v2_mercator_out()        ! line 1086: SPTRUNMV

  ! POLATEV4_grib1
  call test_v1_global_cyl_out()      ! line 1341: SPTRUNV
  call test_v1_polar_stereo_out()    ! line 1368: SPTRUNSV
  call test_v1_mercator_out()        ! line 1401: SPTRUNMV

  print *, "All spectral_interp_mod tests passed."

contains

  logical function all_close(x, n, ref, tol)
    integer, intent(in) :: n
    real,    intent(in) :: x(n), ref, tol
    all_close = all(abs(x - ref) <= tol)
  end function

  ! ====================================================================
  ! POLATES4_grib2
  ! ====================================================================

  subroutine test_s2_iret41_bad_type()
    ! Use polar-stereo (igdtnumi=20) as input -> IDRTI not 0 or 4 -> IRET=41.
    ! Use the standard 360x181 output so GDSWZD reliably sets NO.
    integer, parameter :: bad_len = 18
    integer :: bad_tmpl(bad_len)
    integer, parameter :: km = 1, mo = MO2_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(1,km), go(mo,km), rlat(mo), rlon(mo), li(1,km), lo(mo,km))
    bad_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
        5, 5, 0, 0, 56, 60000000, 0, 500000000, 500000000, 0, 0]
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, 20, bad_tmpl, bad_len, &
        0, out2_ll_tmpl, OUT2_LL_LEN, &
        1, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 41) stop 1
  end subroutine

  subroutine test_s2_error_return_path()
    ! Same bad input; verify lines 482-486 set IBO=1, LO=.FALSE., GO=0.
    integer, parameter :: bad_len = 18
    integer :: bad_tmpl(bad_len)
    integer, parameter :: km = 1, mo = MO2_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(1,km), go(mo,km), rlat(mo), rlon(mo), li(1,km), lo(mo,km))
    bad_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
        5, 5, 0, 0, 56, 60000000, 0, 500000000, 500000000, 0, 0]
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, 20, bad_tmpl, bad_len, &
        0, out2_ll_tmpl, OUT2_LL_LEN, &
        1, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 41) stop 2
    if (ibo(1) /= 1) stop 3
    if (any(lo(1:no,1))) stop 4
    if (any(abs(go(1:no,1)) .gt. 0.00001)) stop 5
  end subroutine

  subroutine test_s2_iscan1()
    ! Line 340: ISCAN=1 causes DLON to be computed as a negative value.
    ! Use a 4-point ISCAN=1 grid with RLON1=RLON2=0 so IG=3 (odd) which
    ! triggers IRET=41 right after line 340 without ever calling SPTRUN.
    ! The SP library does not correctly handle ISKIPI<0 for general inputs.
    integer :: in_iscan1(IN2_LEN)
    integer, parameter :: km = 1, mi = 4*4, mo = MO2_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(mi,km), go(mo,km), rlat(mo), rlon(mo), li(mi,km), lo(mo,km))
    in_iscan1 = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
        4, 4, 0, MISSING, -89500000, 0, &
        48, 89500000, 0, 1000000, 1000000, 192]
    ! scan=192 -> ISCAN=1; RLON1=RLON2=0 -> DLON=-120 -> IG=3 (odd) -> IRET=41
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, IN2_NUM, in_iscan1, IN2_LEN, &
        0, out2_ll_tmpl, OUT2_LL_LEN, &
        mi, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 41) stop 10
  end subroutine

  subroutine test_s2_gaussian_jm_mismatch()
    integer, parameter :: gauss_len = 19
    integer :: gauss_tmpl(gauss_len)
    integer, parameter :: km = 1, mi = 8*5, mo = 1
    integer :: ibi(km), ibo(km), no, iret
    real    :: gi(mi,km), go(mo,km), rlat(mo), rlon(mo)  ! small arrays, no stack issue
    logical*1 :: li(mi,km), lo(mo,km)
    gauss_tmpl = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
        8, 5, 0, MISSING, -89000000, 0, &
        48, 89000000, 315000000, 45000000, 3, 0]
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.; no=1
    call ipolates(IP_SPEC, ipopt, 40, gauss_tmpl, gauss_len, &
        -1, gauss_tmpl, gauss_len, &
        mi, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 41) stop 20
  end subroutine

  subroutine test_s2_nscan1()
    integer :: in_nscan1(IN2_LEN)
    integer, parameter :: km = 1, mo = MO2_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI2,km), go(mo,km), rlat(mo), rlon(mo), li(MI2,km), lo(mo,km))
    in_nscan1 = [6, 255, MISSING, 255, MISSING, 255, MISSING, &
        360, 180, 0, MISSING, -89500000, -180000000, &
        48, 89500000, 179000000, 1000000, 1000000, 96]
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, IN2_NUM, in_nscan1, IN2_LEN, &
        0, out2_ll_tmpl, OUT2_LL_LEN, &
        MI2, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 30
    if (ibo(1) /= 0) stop 31
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 32
  end subroutine

  subroutine test_s2_global_cyl_latlon_out()
    integer, parameter :: km = 1, mo = MO2_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI2,km), go(mo,km), rlat(mo), rlon(mo), li(MI2,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        0, out2_ll_tmpl, OUT2_LL_LEN, &
        MI2, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 40
    if (ibo(1) /= 0) stop 41
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 42
  end subroutine

  subroutine test_s2_global_cyl_gauss_out()
    integer, parameter :: km = 1, mo = MO2_GAUSS
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI2,km), go(mo,km), rlat(mo), rlon(mo), li(MI2,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        40, out2_gauss_tmpl, OUT2_GAUSS_LEN, &
        MI2, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 50
    if (ibo(1) /= 0) stop 51
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 52
  end subroutine

  subroutine test_s2_polar_stereo_out()
    integer, parameter :: km = 1, mo = MO2_PS
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI2,km), go(mo,km), rlat(mo), rlon(mo), li(MI2,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        20, out2_ps_tmpl, OUT2_PS_LEN, &
        MI2, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 60
    if (ibo(1) /= 0) stop 61
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 62
  end subroutine

  subroutine test_s2_mercator_out()
    integer, parameter :: km = 1, mo = MO2_MERC
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI2,km), go(mo,km), rlat(mo), rlon(mo), li(MI2,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        10, out2_merc_tmpl, OUT2_MERC_LEN, &
        MI2, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 70
    if (ibo(1) /= 0) stop 71
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 72
  end subroutine

  ! ====================================================================
  ! POLATES4_grib1
  ! ====================================================================

  subroutine test_s1_global_cyl_out()
    integer, parameter :: km = 1, mo = MO1_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI1,km), go(mo,km), rlat(mo), rlon(mo), li(MI1,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, kgdsi, out1_ll_kgds, &
        MI1, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 80
    if (ibo(1) /= 0) stop 81
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 82
  end subroutine

  subroutine test_s1_polar_stereo_out()
    integer, parameter :: km = 1, mo = MO1_PS
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI1,km), go(mo,km), rlat(mo), rlon(mo), li(MI1,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, kgdsi, out1_ps_kgds, &
        MI1, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 90
    if (ibo(1) /= 0) stop 91
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 92
  end subroutine

  subroutine test_s1_mercator_out()
    integer, parameter :: km = 1, mo = MO1_MERC
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: gi(:,:), go(:,:), rlat(:), rlon(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(gi(MI1,km), go(mo,km), rlat(mo), rlon(mo), li(MI1,km), lo(mo,km))
    ibi=0; gi=1.0; li=.false.; rlat=0.; rlon=0.
    call ipolates(IP_SPEC, ipopt, kgdsi, out1_merc_kgds, &
        MI1, mo, km, ibi, li, gi, no, rlat, rlon, ibo, lo, go, iret)
    if (iret /= 0) stop 100
    if (ibo(1) /= 0) stop 101
    if (.not. all_close(go(:,1), mo, 1.0, 0.1)) stop 102
  end subroutine

  ! ====================================================================
  ! POLATEV4_grib2  (also covers lines 148-155 via interpolate_spectral_vector)
  ! ====================================================================

  subroutine test_v2_global_cyl_out()
    integer, parameter :: km = 1, mo = MO2_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: ui(:,:), vi(:,:), uo(:,:), vo(:,:)
    real,      allocatable :: rlat(:), rlon(:), crot(:), srot(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(ui(MI2,km), vi(MI2,km), uo(mo,km), vo(mo,km))
    allocate(rlat(mo), rlon(mo), crot(mo), srot(mo), li(MI2,km), lo(mo,km))
    ibi=0; ui=1.0; vi=0.0; li=.false.; rlat=0.; rlon=0.
    call ipolatev(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        0, out2_ll_tmpl, OUT2_LL_LEN, &
        MI2, mo, km, ibi, li, ui, vi, &
        no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    if (iret /= 0) stop 110
  end subroutine

  subroutine test_v2_polar_stereo_out()
    integer, parameter :: km = 1, mo = MO2_PS
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: ui(:,:), vi(:,:), uo(:,:), vo(:,:)
    real,      allocatable :: rlat(:), rlon(:), crot(:), srot(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(ui(MI2,km), vi(MI2,km), uo(mo,km), vo(mo,km))
    allocate(rlat(mo), rlon(mo), crot(mo), srot(mo), li(MI2,km), lo(mo,km))
    ibi=0; ui=1.0; vi=0.0; li=.false.; rlat=0.; rlon=0.
    call ipolatev(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        20, out2_ps_tmpl, OUT2_PS_LEN, &
        MI2, mo, km, ibi, li, ui, vi, &
        no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    if (iret /= 0) stop 120
  end subroutine

  subroutine test_v2_mercator_out()
    integer, parameter :: km = 1, mo = MO2_MERC
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: ui(:,:), vi(:,:), uo(:,:), vo(:,:)
    real,      allocatable :: rlat(:), rlon(:), crot(:), srot(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(ui(MI2,km), vi(MI2,km), uo(mo,km), vo(mo,km))
    allocate(rlat(mo), rlon(mo), crot(mo), srot(mo), li(MI2,km), lo(mo,km))
    ibi=0; ui=1.0; vi=0.0; li=.false.; rlat=0.; rlon=0.
    call ipolatev(IP_SPEC, ipopt, IN2_NUM, in2_tmpl, IN2_LEN, &
        10, out2_merc_tmpl, OUT2_MERC_LEN, &
        MI2, mo, km, ibi, li, ui, vi, &
        no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    if (iret /= 0) stop 130
  end subroutine

  ! ====================================================================
  ! POLATEV4_grib1
  ! ====================================================================

  subroutine test_v1_global_cyl_out()
    integer, parameter :: km = 1, mo = MO1_LL
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: ui(:,:), vi(:,:), uo(:,:), vo(:,:)
    real,      allocatable :: rlat(:), rlon(:), crot(:), srot(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(ui(MI1,km), vi(MI1,km), uo(mo,km), vo(mo,km))
    allocate(rlat(mo), rlon(mo), crot(mo), srot(mo), li(MI1,km), lo(mo,km))
    ibi=0; ui=1.0; vi=0.0; li=.false.; rlat=0.; rlon=0.
    call ipolatev(IP_SPEC, ipopt, kgdsi, out1_ll_kgds, &
        MI1, mo, km, ibi, li, ui, vi, &
        no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    if (iret /= 0) stop 140
  end subroutine

  subroutine test_v1_polar_stereo_out()
    integer, parameter :: km = 1, mo = MO1_PS
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: ui(:,:), vi(:,:), uo(:,:), vo(:,:)
    real,      allocatable :: rlat(:), rlon(:), crot(:), srot(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(ui(MI1,km), vi(MI1,km), uo(mo,km), vo(mo,km))
    allocate(rlat(mo), rlon(mo), crot(mo), srot(mo), li(MI1,km), lo(mo,km))
    ibi=0; ui=1.0; vi=0.0; li=.false.; rlat=0.; rlon=0.
    call ipolatev(IP_SPEC, ipopt, kgdsi, out1_ps_kgds, &
        MI1, mo, km, ibi, li, ui, vi, &
        no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    if (iret /= 0) stop 150
  end subroutine

  subroutine test_v1_mercator_out()
    integer, parameter :: km = 1, mo = MO1_MERC
    integer :: ibi(km), ibo(km), no, iret
    real,      allocatable :: ui(:,:), vi(:,:), uo(:,:), vo(:,:)
    real,      allocatable :: rlat(:), rlon(:), crot(:), srot(:)
    logical*1, allocatable :: li(:,:), lo(:,:)
    allocate(ui(MI1,km), vi(MI1,km), uo(mo,km), vo(mo,km))
    allocate(rlat(mo), rlon(mo), crot(mo), srot(mo), li(MI1,km), lo(mo,km))
    ibi=0; ui=1.0; vi=0.0; li=.false.; rlat=0.; rlon=0.
    call ipolatev(IP_SPEC, ipopt, kgdsi, out1_merc_kgds, &
        MI1, mo, km, ibi, li, ui, vi, &
        no, rlat, rlon, crot, srot, ibo, lo, uo, vo, iret)
    if (iret /= 0) stop 160
  end subroutine

end program test_spectral_interp_mod
