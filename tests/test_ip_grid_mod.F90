program test_ip_grid_mod
  !! Unit tests for ip_grid_mod.F90, specifically the field_pos() function
  !!
  !! Tests all scanning modes (nscan = 0, 1, 2, 3), boundary conditions,
  !! wraparound scenarios, and position calculations.
  use ip_grid_mod
  use ip_equid_cylind_grid_mod
  use ip_grid_descriptor_mod
  implicit none

  integer :: test_count = 0
  integer :: fail_count = 0

  print *,'Testing ip_grid_mod field_pos() function...'
  print *,''

  ! Test nscan=0 (row-major scanning: x first, then y)
  call test_nscan_0()

  ! Test nscan=1 (column-major scanning: y first, then x)
  call test_nscan_1()

  ! Test nscan=2 (staggered diagonal - complex calculation)
  call test_nscan_2()

  ! Test nscan=3 (staggered diagonal - alternative calculation)
  call test_nscan_3()

  ! Test boundary conditions (out-of-bounds returns 0)
  call test_boundary_conditions()

  ! Test wraparound scenarios
  call test_wraparound()

  print *,''
  print *,'========================================='
  print *,'Tests run:', test_count
  print *,'Tests failed:', fail_count
  print *,'========================================='

  if (fail_count > 0) then
    stop 1
  endif

contains

  subroutine create_test_grid(grid, im, jm, nscan_value, kscan_value, iwrap_val, jwrap1_val, jwrap2_val)
    !! Create a minimal test grid with specified parameters
    !! Note: We manually set grid parameters to avoid complex GRIB initialization
    class(ip_grid), allocatable, intent(out) :: grid
    integer, intent(in) :: im, jm, nscan_value, kscan_value, iwrap_val, jwrap1_val, jwrap2_val
    type(ip_equid_cylind_grid), allocatable :: temp_grid

    allocate(temp_grid)
    temp_grid%im = im
    temp_grid%jm = jm
    temp_grid%nm = im * jm
    temp_grid%nscan = nscan_value
    temp_grid%nscan_field_pos = nscan_value
    temp_grid%kscan = kscan_value
    temp_grid%iwrap = iwrap_val
    temp_grid%jwrap1 = jwrap1_val
    temp_grid%jwrap2 = jwrap2_val
    temp_grid%rerth = 6.3712E6
    temp_grid%eccen_squared = 0.0

    call move_alloc(temp_grid, grid)
  end subroutine create_test_grid

  subroutine test_nscan_0()
    !! Test nscan=0 (row-major: field_pos = ii + (jj-1)*im)
    class(ip_grid), allocatable :: grid
    integer :: pos

    print *,'Test nscan=0 (row-major scanning)...'

    ! Create grid: 5x4 (im=5, jm=4)
    call create_test_grid(grid, 5, 4, 0, 0, 0, 0, 0)

    ! Test corner points
    test_count = test_count + 1
    pos = grid%field_pos(1, 1)
    if (pos /= 1) then
      print *,'  FAIL: (1,1) should be 1, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(5, 1)
    if (pos /= 5) then
      print *,'  FAIL: (5,1) should be 5, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(1, 2)
    if (pos /= 6) then  ! 1 + (2-1)*5
      print *,'  FAIL: (1,2) should be 6, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(5, 4)
    if (pos /= 20) then  ! 5 + (4-1)*5
      print *,'  FAIL: (5,4) should be 20, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(3, 3)
    if (pos /= 13) then  ! 3 + (3-1)*5
      print *,'  FAIL: (3,3) should be 13, got', pos
      fail_count = fail_count + 1
    endif

    deallocate(grid)
    print *,'  PASS: nscan=0 tests completed'
    print *,''

  end subroutine test_nscan_0

  subroutine test_nscan_1()
    !! Test nscan=1 (column-major: field_pos = jj + (ii-1)*jm)
    class(ip_grid), allocatable :: grid
    integer :: pos

    print *,'Test nscan=1 (column-major scanning)...'

    ! Create grid: 5x4 (im=5, jm=4)
    call create_test_grid(grid, 5, 4, 1, 0, 0, 0, 0)

    test_count = test_count + 1
    pos = grid%field_pos(1, 1)
    if (pos /= 1) then
      print *,'  FAIL: (1,1) should be 1, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(1, 4)
    if (pos /= 4) then
      print *,'  FAIL: (1,4) should be 4, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(2, 1)
    if (pos /= 5) then  ! 1 + (2-1)*4
      print *,'  FAIL: (2,1) should be 5, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(5, 4)
    if (pos /= 20) then  ! 4 + (5-1)*4
      print *,'  FAIL: (5,4) should be 20, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(3, 2)
    if (pos /= 10) then  ! 2 + (3-1)*4
      print *,'  FAIL: (3,2) should be 10, got', pos
      fail_count = fail_count + 1
    endif

    deallocate(grid)
    print *,'  PASS: nscan=1 tests completed'
    print *,''

  end subroutine test_nscan_1

  subroutine test_nscan_2()
    !! Test nscan=2 (staggered diagonal with 2*im-1 width)
    !! field_pos = (iif + (jjf-1)*(2*im-1) + 1 - kscan) / 2
    !! where iif = jj + (ii - is1), jjf = jj - (ii - is1) + kscan, is1 = (jm+1-kscan)/2
    class(ip_grid), allocatable :: grid
    integer :: pos

    print *,'Test nscan=2 (staggered diagonal, mass points)...'

    ! Create grid: 6x4 with kscan=0 (mass points)
    call create_test_grid(grid, 6, 4, 2, 0, 0, 0, 0)

    ! For nscan=2 with im=6, jm=4, kscan=0:
    ! is1 = (4 + 1 - 0) / 2 = 2 (integer division)
    ! Test point (2, 2):
    !   iif = 2 + (2 - 2) = 2
    !   jjf = 2 - (2 - 2) + 0 = 2
    !   if (iif >= 1 and iif <= 2*6-1=11 and jjf >= 1 and jjf <= 4)
    !   field_pos = (2 + (2-1)*(2*6-1) + 1 - 0) / 2 = (2 + 11 + 1) / 2 = 7

    test_count = test_count + 1
    pos = grid%field_pos(2, 2)
    if (pos /= 7) then
      print *,'  FAIL: (2,2) should be 7, got', pos
      fail_count = fail_count + 1
    endif

    ! Test point (1, 1):
    !   iif = 1 + (1 - 2) = 0 (out of bounds)
    !   Should return 0
    test_count = test_count + 1
    pos = grid%field_pos(1, 1)
    if (pos /= 0) then
      print *,'  FAIL: (1,1) out of bounds should be 0, got', pos
      fail_count = fail_count + 1
    endif

    deallocate(grid)
    print *,'  PASS: nscan=2 tests completed'
    print *,''

  end subroutine test_nscan_2

  subroutine test_nscan_3()
    !! Test nscan=3 (staggered diagonal alternative)
    !! field_pos = (iif+1)/2 + (jjf-1)*im
    !! where iif = jj + (ii - is1), jjf = jj - (ii - is1) + kscan, is1 = (jm+1-kscan)/2
    class(ip_grid), allocatable :: grid
    integer :: pos

    print *,'Test nscan=3 (staggered diagonal, wind points)...'

    ! Create grid: 6x4 with kscan=1 (wind points - E-stagger)
    call create_test_grid(grid, 6, 4, 3, 1, 0, 0, 0)

    ! For nscan=3 with im=6, jm=4, kscan=1:
    ! is1 = (4 + 1 - 1) / 2 = 2
    ! Test point (2, 2):
    !   iif = 2 + (2 - 2) = 2
    !   jjf = 2 - (2 - 2) + 1 = 3
    !   if (iif >= 1 and iif <= 2*6-1=11 and jjf >= 1 and jjf <= 4)
    !   field_pos = (2+1)/2 + (3-1)*6 = 1 + 12 = 13

    test_count = test_count + 1
    pos = grid%field_pos(2, 2)
    if (pos /= 13) then
      print *,'  FAIL: (2,2) should be 13, got', pos
      fail_count = fail_count + 1
    endif

    ! Test point (3, 3):
    !   iif = 3 + (3 - 2) = 4
    !   jjf = 3 - (3 - 2) + 1 = 3
    !   field_pos = (4+1)/2 + (3-1)*6 = 2 + 12 = 14
    test_count = test_count + 1
    pos = grid%field_pos(3, 3)
    if (pos /= 14) then
      print *,'  FAIL: (3,3) should be 14, got', pos
      fail_count = fail_count + 1
    endif

    deallocate(grid)
    print *,'  PASS: nscan=3 tests completed'
    print *,''

  end subroutine test_nscan_3

  subroutine test_boundary_conditions()
    !! Test that out-of-bounds coordinates return 0
    class(ip_grid), allocatable :: grid
    integer :: pos

    print *,'Test boundary conditions...'

    call create_test_grid(grid, 5, 4, 0, 0, 0, 0, 0)

    ! Test points outside grid bounds
    test_count = test_count + 1
    pos = grid%field_pos(0, 2)  ! i < 1
    if (pos /= 0) then
      print *,'  FAIL: (0,2) should return 0, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(6, 2)  ! i > im
    if (pos /= 0) then
      print *,'  FAIL: (6,2) should return 0, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(2, 0)  ! j < 1
    if (pos /= 0) then
      print *,'  FAIL: (2,0) should return 0, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(2, 5)  ! j > jm
    if (pos /= 0) then
      print *,'  FAIL: (2,5) should return 0, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(-1, -1)  ! Both out of bounds
    if (pos /= 0) then
      print *,'  FAIL: (-1,-1) should return 0, got', pos
      fail_count = fail_count + 1
    endif

    ! Test valid edge points
    test_count = test_count + 1
    pos = grid%field_pos(1, 1)  ! Lower-left corner
    if (pos /= 1) then
      print *,'  FAIL: (1,1) should be 1, got', pos
      fail_count = fail_count + 1
    endif

    test_count = test_count + 1
    pos = grid%field_pos(5, 4)  ! Upper-right corner
    if (pos /= 20) then
      print *,'  FAIL: (5,4) should be 20, got', pos
      fail_count = fail_count + 1
    endif

    deallocate(grid)
    print *,'  PASS: boundary condition tests completed'
    print *,''

  end subroutine test_boundary_conditions

  subroutine test_wraparound()
    !! Test wraparound scenarios (periodic boundary conditions)
    class(ip_grid), allocatable :: grid
    integer :: pos

    print *,'Test wraparound scenarios...'

    ! Test x-wraparound (iwrap > 0, no y-wraparound)
    ! Create grid: 10x5 with x-wraparound at 10
    call create_test_grid(grid, 10, 5, 0, 0, 10, 0, 0)

    ! For nscan=0 with wraparound:
    ! Input (11, 2) wraps to (1, 2) -> field_pos = 1 + (2-1)*10 = 11
    test_count = test_count + 1
    pos = grid%field_pos(11, 2)
    if (pos /= 11) then
      print *,'  FAIL: wraparound (11,2) should be 11, got', pos
      fail_count = fail_count + 1
    endif

    ! Input (0, 2) wraps to (10, 2) -> field_pos = 10 + (2-1)*10 = 20
    test_count = test_count + 1
    pos = grid%field_pos(0, 2)
    if (pos /= 20) then
      print *,'  FAIL: wraparound (0,2) should be 20, got', pos
      fail_count = fail_count + 1
    endif

    ! Input (12, 3) wraps to (2, 3) -> field_pos = 2 + (3-1)*10 = 22
    test_count = test_count + 1
    pos = grid%field_pos(12, 3)
    if (pos /= 22) then
      print *,'  FAIL: wraparound (12,3) should be 22, got', pos
      fail_count = fail_count + 1
    endif

    ! Test with y-wraparound (polar grid scenario)
    deallocate(grid)
    call create_test_grid(grid, 8, 6, 0, 0, 8, 1, 6)

    ! For j < 1 with jwrap1 > 0:
    ! jj = jwrap1 - j, ii = mod(ii-1+iwrap/2, iwrap)+1
    ! Input (1, 0) with jwrap1=1:
    !   jj = 1 - 0 = 1
    !   ii = mod(1-1+8/2, 8)+1 = mod(4, 8)+1 = 5
    !   field_pos = 5 + (1-1)*8 = 5
    test_count = test_count + 1
    pos = grid%field_pos(1, 0)
    if (pos /= 5) then
      print *,'  FAIL: y-wraparound (1,0) should be 5, got', pos
      fail_count = fail_count + 1
    endif

    ! For j > jm with jwrap2 > 0:
    ! jj = jwrap2 - j
    ! Input (1, 7) with jwrap2=6 -> jj = 6 - 7 = -1 (out of bounds)
    test_count = test_count + 1
    pos = grid%field_pos(1, 7)
    if (pos /= 0) then
      print *,'  FAIL: y-wraparound (1,7) out of bounds should be 0, got', pos
      fail_count = fail_count + 1
    endif

    deallocate(grid)
    print *,'  PASS: wraparound tests completed'
    print *,''

  end subroutine test_wraparound

end program test_ip_grid_mod
