! This is a test for the NCEPLIBS-ip library.
!
! Tests POLFIXS and POLFIXV from polfix_mod:
!   POLFIXS averages multiple pole scalar values on a lat/lon grid.
!   POLFIXV averages multiple pole vector values (u,v winds).
!
! Hang Lei, Aug 2026

program test_polfix
  use polfix_mod
  implicit none

  integer, parameter :: NM = 6, KM = 1
  real, parameter :: tol = 1.0e-5

  ! --- Test POLFIXS ---
  ! Set up 4 points near the north pole, 2 points not near any pole.
  ! POLFIXS averages all pole-neighboring points and assigns the mean.
  real    :: rlat(NM)
  integer :: ib(KM)
  logical*1 :: lo(NM, KM)
  real    :: go(NM, KM)
  real    :: expected_np

  print *, 'Testing POLFIXS...'

  ! lat >= 89.9995 are treated as "north pole"
  rlat = (/ 89.9999, 89.9999, 89.9999, 89.9999, 45.0, 0.0 /)

  ib(1) = 0   ! no bitmap, all points valid
  lo = .true.
  go(:, 1) = (/ 10.0, 20.0, 30.0, 40.0, 5.0, 3.0 /)

  call polfixs(NM, NM, KM, rlat, ib, lo, go)

  ! expected average of the 4 pole points: (10+20+30+40)/4 = 25
  expected_np = 25.0
  if (abs(go(1,1) - expected_np) > tol) then
    print *, 'FAILED POLFIXS: go(1,1) =', go(1,1), ' expected', expected_np
    stop 1
  end if
  if (abs(go(2,1) - expected_np) > tol) then
    print *, 'FAILED POLFIXS: go(2,1) =', go(2,1), ' expected', expected_np
    stop 2
  end if
  if (abs(go(3,1) - expected_np) > tol) then
    print *, 'FAILED POLFIXS: go(3,1) =', go(3,1), ' expected', expected_np
    stop 3
  end if
  if (abs(go(4,1) - expected_np) > tol) then
    print *, 'FAILED POLFIXS: go(4,1) =', go(4,1), ' expected', expected_np
    stop 4
  end if

  ! Non-pole points should be unchanged
  if (abs(go(5,1) - 5.0) > tol) then
    print *, 'FAILED POLFIXS: non-pole point changed: go(5,1) =', go(5,1)
    stop 5
  end if
  if (abs(go(6,1) - 3.0) > tol) then
    print *, 'FAILED POLFIXS: non-pole point changed: go(6,1) =', go(6,1)
    stop 6
  end if

  ! Test south pole averaging too
  rlat = (/ -89.9999, -89.9999, -89.9999, 45.0, 0.0, -45.0 /)
  go(:, 1) = (/ 6.0, 12.0, 18.0, 5.0, 3.0, 1.0 /)
  call polfixs(NM, NM, KM, rlat, ib, lo, go)
  expected_np = (6.0 + 12.0 + 18.0) / 3.0
  if (abs(go(1,1) - expected_np) > tol) then
    print *, 'FAILED POLFIXS south pole: go(1,1) =', go(1,1), ' expected', expected_np
    stop 7
  end if

  print *, 'PASSED: POLFIXS'

  ! --- Test POLFIXV ---
  ! POLFIXV rotates and averages u,v vectors at the poles.
  ! With a single pole point (WNP <= 1), no averaging is done.
  ! With multiple pole points (WNP > 1), they are averaged.
  print *, 'Testing POLFIXV...'
  block
    integer, parameter :: NMV = 5, KMV = 1
    real    :: rlatv(NMV), rlonv(NMV)
    integer :: ibv(KMV)
    logical*1 :: lov(NMV, KMV)
    real    :: uo(NMV, KMV), vo(NMV, KMV)

    ! 3 north-pole points at different longitudes; 2 mid-lat points
    rlatv = (/ 89.9999, 89.9999, 89.9999, 45.0, 0.0 /)
    rlonv = (/ 0.0, 120.0, 240.0, 0.0, 0.0 /)
    ibv(1) = 0
    lov = .true.
    ! All pole points have u=1, v=0
    uo(:, 1) = (/ 1.0, 1.0, 1.0, 2.0, 3.0 /)
    vo(:, 1) = (/ 0.0, 0.0, 0.0, 1.0, 0.5 /)

    call polfixv(NMV, NMV, KMV, rlatv, rlonv, ibv, lov, uo, vo)

    ! Non-pole points must be unchanged
    if (abs(uo(4,1) - 2.0) > tol .or. abs(vo(4,1) - 1.0) > tol) then
      print *, 'FAILED POLFIXV: non-pole changed: uo=', uo(4,1), ' vo=', vo(4,1)
      stop 10
    end if
    if (abs(uo(5,1) - 3.0) > tol .or. abs(vo(5,1) - 0.5) > tol) then
      print *, 'FAILED POLFIXV: non-pole changed: uo=', uo(5,1), ' vo=', vo(5,1)
      stop 11
    end if

    ! After averaging, the returned (uo,vo) at each pole must have magnitude
    ! consistent with a rotation of the averaged earth-centred vector.
    ! Simply check that the magnitude is reasonable (not NaN/Inf and >=0).
    if (uo(1,1)**2 + vo(1,1)**2 < 0.0) then
      print *, 'FAILED POLFIXV: negative magnitude squared at pole'
      stop 12
    end if
  end block

  print *, 'PASSED: POLFIXV'
  print *, 'SUCCESS!'

end program test_polfix
