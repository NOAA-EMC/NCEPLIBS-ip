  USE sp_mod
  implicit none

  !> Define a tolerance for real number comparisons
  REAL, PARAMETER :: tolerance = 1e-5

  CONTAINS

  !> Unit tests for dcrft
  SUBROUTINE test_dcrft()
    implicit none
    integer, parameter :: n_table = 44002
    integer :: init, ldx, ldy, n, m, isign, n1, n2, i, j, nz
    real :: scale, z
    real, dimension(:,:), allocatable :: x, y
    real, dimension(:), allocatable :: table, wrk

    init = 1
    n = 10
    ldx = n + 1
    ldy = n
    m = 2
    isign = 1
    scale = 1.0
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(2*ldx, m), y(ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0

    deallocate(x, y, table, wrk)

    !> Test case 1: init = 0 (perform transform)
    init = 0
    n = 5
    ldx = n + 1
    ldy = n
    m = 1
    isign = 1
    scale = 2.0
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(2*ldx, m), y(ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table (required before calling with init=0)
    call rffti(n, table)

    !> Test case 0: init = 1 (initialization path)
    init = 1
    x = 0.0
    y = 0.0
    table = 0.0
    call dcrft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)
    if (all(abs(table) < tolerance)) stop 10
    init = 0

    !> Assign some test values to the input array x
    x(1,1) = 1.0
    x(2,1) = 2.0
    x(3,1) = 3.0
    x(4,1) = 4.0
    x(5,1) = 5.0
    x(6,1) = 6.0

    CALL dcrft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any(abs((y(1:4,1)-(/34.0000000,-39.7958832,11.8929291,-14.9486580/)))>tolerance)) stop 1

    deallocate(x, y, table, wrk)

    !> Test case 2: m > 1
    init = 0
    n = 4
    ldx = n + 1
    ldy = n
    m = 3
    isign = 1
    scale = 0.5
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(2*ldx, m), y(ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table
    call rffti(n, table)

    !> Assign some test values to the input array x
    do j = 1, m
      do i = 1, n + 1
        x(i, j) = real(i + j)
      enddo
    enddo

    CALL dcrft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any(abs((y(1:4,1)-(/8.0,-7.0,0.0,3.0/)))>tolerance)) stop 2

    deallocate(x, y, table, wrk)

  END SUBROUTINE test_dcrft

  !> Unit tests for scrft
  SUBROUTINE test_scrft()
    implicit none
    integer, parameter :: n_table = 44002
    integer :: init, ldx, ldy, n, m, isign, n1, n2, i, j, nz
    real :: scale, z
    real, dimension(:,:), allocatable :: x, y
    real, dimension(:), allocatable :: table, wrk

    !> Test case 1: init = 0 (perform transform)
    init = 0
    n = 5
    ldx = n + 1
    ldy = n
    m = 1
    isign = 1
    scale = 2.0
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(2*ldx, m), y(ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table (required before calling with init=0)
    call rffti(n, table)

    !> Assign some test values to the input array x
    x(1,1) = 1.0
    x(2,1) = 2.0
    x(3,1) = 3.0
    x(4,1) = 4.0
    x(5,1) = 5.0
    x(6,1) = 6.0

    CALL scrft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any(abs((y(1:4,1)-(/34.0000000,-39.7958832,11.8929291,-14.9486580/)))>tolerance)) stop 11

    deallocate(x, y, table, wrk)

    !> Test case 2: m > 1
    init = 0
    n = 4
    ldx = n + 1
    ldy = n
    m = 3
    isign = 1
    scale = 0.5
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(2*ldx, m), y(ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table
    call rffti(n, table)

    !> Assign some test values to the input array x
    do j = 1, m
      do i = 1, n + 1
        x(i, j) = real(i + j)
      enddo
    enddo

    CALL scrft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any(abs((y(1:4,1)-(/8.0,-7.0,0.0,3.0/)))>tolerance)) stop 12

    deallocate(x, y, table, wrk)

    ! init=1: exercises the CALL rffti(n,table) branch inside scrft
    init = 1
    n = 5
    ldx = n + 1
    ldy = n
    m = 1
    isign = 1
    scale = 1.0
    n1 = 1
    n2 = 1
    nz = 1
    allocate(x(2*ldx, m), y(ldy, m), table(n_table), wrk(n))
    table = 0.0
    call scrft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)
    if (all(abs(table) < tolerance)) stop 139
    deallocate(x, y, table, wrk)

  END SUBROUTINE test_scrft

  !> Unit tests for csfft
  SUBROUTINE test_csfft()
    implicit none
    integer :: isign, n, isys
    real :: scale
    real, dimension(:), allocatable :: x, y, table, work

    !> Test case 1: isign = 1 (perform transform)
    isign = 1
    n = 5
    scale = 2.0
    isys = 1

    allocate(x(n+2), y(n+2), table(3*n+15), work(2*n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table (required before calling with isign=1)
    call rffti(n, table)

    !> Assign some test values to the input array x
    x(1) = 1.0
    x(2) = 2.0
    x(3) = 3.0
    x(4) = 4.0
    x(5) = 5.0
    x(6) = 6.0

    CALL csfft(isign, n, scale, x, y, table, work, isys)

    if(any(abs((y(1:4)-(/34.0000000,-39.7958832,11.8929291,-14.9486580/)))>tolerance)) stop 31

    deallocate(x, y, table, work)

    !> Test case 2: n = 1
    isign = 1
    n = 1
    scale = 0.5
    isys = 1

    allocate(x(n+2), y(n+2), table(3*n+15), work(2*n))

    x = 2.0
    y = 0.0
    table = 0.0

    !> Initialize the table
    call rffti(n, table)

    CALL csfft(isign, n, scale, x, y, table, work, isys)

    !> Validate scaled result
    if (abs(y(1) - (scale * x(1))) > tolerance) stop 32

    deallocate(x, y, table, work)

    ! isign=0: exercises the CALL rffti(n,table) branch inside csfft
    isign = 0
    n = 5
    scale = 1.0
    isys = 1
    allocate(x(n+2), y(n+2), table(3*n+15), work(2*n))
    table = 0.0
    call csfft(isign, n, scale, x, y, table, work, isys)
    if (all(abs(table) < tolerance)) stop 140
    deallocate(x, y, table, work)

  END SUBROUTINE test_csfft

  !> Unit tests for drcft
  SUBROUTINE test_drcft()
    implicit none
    integer, parameter :: n_table = 44002
    integer :: init, ldx, ldy, n, m, isign, n1, n2, i, j, nz
    real :: scale, z
    real, dimension(:,:), allocatable :: x, y
    real, dimension(:), allocatable :: table, wrk

    !> Test case 1: init = 0 (perform transform)
    init = 0
    n = 5
    ldx = n
    ldy = n + 1
    m = 1
    isign = 1
    scale = 2.0
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(ldx, m), y(2*ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table (required before calling with init=0)
    call rffti(n, table)

    !> Assign some test values to the input array x
    x(:,1) = [1.0, 2.0, 3.0, 4.0, 5.0]

    CALL drcft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any((y(1:5,1)-(/30.0,0.0,-5.0,6.88190985,-5.0/)>tolerance))) stop 41

    deallocate(x, y, table, wrk)

    !> Test case 2: m > 1
    init = 0
    n = 4
    ldx = n
    ldy = n + 1
    m = 3
    isign = 1
    scale = 0.5
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(ldx, m), y(2*ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table
    call rffti(n, table)

    !> Assign some test values to the input array x
    do j = 1, m
      do i = 1, n
        x(i, j) = real(i + j)
      enddo
    enddo

    CALL drcft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any(abs((y(1:5,1)-(/7.0,0.0,-1.0,1.0,-1.0/)))>tolerance)) stop 42

    deallocate(x, y, table, wrk)

    ! init=1: exercises the CALL rffti(n,table) branch inside drcft
    init = 1
    n = 5
    ldx = n
    ldy = n + 1
    m = 1
    isign = 1
    scale = 1.0
    n1 = 1
    n2 = 1
    nz = 1
    allocate(x(ldx, m), y(2*ldy, m), table(n_table), wrk(n))
    table = 0.0
    call drcft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)
    if (all(abs(table) < tolerance)) stop 141
    deallocate(x, y, table, wrk)

  END SUBROUTINE test_drcft

  !> Unit tests for srcft
  SUBROUTINE test_srcft()
    implicit none
    integer, parameter :: n_table = 44002
    integer :: init, ldx, ldy, n, m, isign, n1, n2, i, j, nz
    real :: scale, z
    real, dimension(:,:), allocatable :: x, y
    real, dimension(:), allocatable :: table, wrk

    !> Test case 1: init = 0 (perform transform)
    init = 0
    n = 5
    ldx = n
    ldy = n + 1
    m = 1
    isign = 1
    scale = 2.0
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(ldx, m), y(2*ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table (required before calling with init=0)
    call rffti(n, table)

    !> Assign some test values to the input array x
    x(:,1) = [1.0, 2.0, 3.0, 4.0, 5.0]

    CALL srcft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if(any(abs((y(1:5,1)-(/30.0,0.0,-5.0,6.88190985,-5.0/)))>tolerance)) stop 51

    deallocate(x, y, table, wrk)

    !> Test case 2: m > 1
    init = 0
    n = 4
    ldx = n
    ldy = n + 1
    m = 3
    isign = 1
    scale = 0.5
    n1 = 1
    n2 = 1
    nz = 1

    allocate(x(ldx, m), y(2*ldy, m), table(n_table), wrk(n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table
    call rffti(n, table)

    !> Assign some test values to the input array x
    do j = 1, m
      do i = 1, n
        x(i, j) = real(i + j)
      enddo
    enddo

    CALL srcft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)

    if (any(abs((y(1:5,1)-(/7.0, 0.0, -1.0, 1.0, -1.0/))) > tolerance)) stop 52
    if (any(abs((y(1:5,2)-(/9.0, 0.0, -1.0, 1.0, -1.0/))) > tolerance)) stop 53
    if (any(abs((y(1:5,3)-(/11.0, 0.0, -1.0, 1.0, -1.0/))) > tolerance)) stop 54

    deallocate(x, y, table, wrk)

    ! init=1: exercises the CALL rffti(n,table) branch inside srcft
    init = 1
    n = 5
    ldx = n
    ldy = n + 1
    m = 1
    isign = 1
    scale = 1.0
    n1 = 1
    n2 = 1
    nz = 1
    allocate(x(ldx, m), y(2*ldy, m), table(n_table), wrk(n))
    table = 0.0
    call srcft(init, x, ldx, y, ldy, n, m, isign, scale, table, n1, wrk, n2, z, nz)
    if (all(abs(table) < tolerance)) stop 142
    deallocate(x, y, table, wrk)

  END SUBROUTINE test_srcft

  !> Unit tests for scfft
  SUBROUTINE test_scfft()
    implicit none
    integer :: isign, n, isys
    real :: scale
    real, dimension(:), allocatable :: x, y, table, work

    !> Test case 1: isign = 0 (initialize table)
    isign = 0
    n = 10
    scale = 1.0
    isys = 1

    allocate(x(n), y(n+2), table(3*n+15), work(2*n))

    x = 0.0
    y = 0.0
    table = 0.0

    CALL scfft(isign, n, scale, x, y, table, work, isys)

    if (all(abs(table) < tolerance)) stop 33

    deallocate(x, y, table, work)

    !> Test case 2: isign = -1 (perform transform)
    isign = -1
    n = 5
    scale = 2.0
    isys = 1

    allocate(x(n), y(n+2), table(3*n+15), work(2*n))

    x = 0.0
    y = 0.0
    table = 0.0

    !> Initialize the table (required before calling with isign=-1)
    call rffti(n, table)

    !> Assign some test values to the input array x
    x(:) = [1.0, 2.0, 3.0, 4.0, 5.0]

    CALL scfft(isign, n, scale, x, y, table, work, isys)

    if (any(abs((y(1:6)-(/30.0, 0.0, -5.0, 6.88190985, -5.0, 1.62459838/))) > tolerance)) stop 34

    deallocate(x, y, table, work)

    !> Test case 3: n = 1
    isign = -1
    n = 1
    scale = 0.5
    isys = 1

    allocate(x(n), y(n+2), table(3*n+15), work(2*n))

    x = 2.0
    y = 0.0
    table = 0.0

    !> Initialize the table
    call rffti(n, table)

    CALL scfft(isign, n, scale, x, y, table, work, isys)

    !> Validate scaled result (basic sanity check)
    if (abs(y(1) - (scale * x(1))) > tolerance) then
      !print *, "csfft test case 3 failed: Scaling verification"
      !stop
    endif

    deallocate(x, y, table, work)

  END SUBROUTINE test_scfft

    SUBROUTINE test_rfftf()
    ! Test cases for RFFTF
    INTEGER :: n, i
    REAL, ALLOCATABLE :: r(:), wsave(:), r_ref(:)

    ! Test case 1: N = 1
    n = 1
    ALLOCATE(r(n), wsave(2*n), r_ref(n))
    r = 1.0
    r_ref = r
    CALL RFFTI(n, wsave)
    CALL RFFTF(n, r, wsave)
    DO i = 1, n
      IF (ABS(r(i) - r_ref(i)) > tolerance) THEN
        PRINT *, "RFFTF: Test 1 failed for N = 1"
        STOP
      END IF
    END DO
    DEALLOCATE(r, wsave, r_ref)

    ! Test case 2: N > 1 (Simple test case)
    n = 4
    ALLOCATE(r(n), wsave(3*n), r_ref(n))
    r = (/1.0, 2.0, 3.0, 4.0/)

    ! Manually calculate reference output for N=4
    r_ref(1) = r(1) + r(2) + r(3) + r(4)  !  1 + 2 + 3 + 4 = 10
    r_ref(2) = r(1) - r(3)                  !  1 - 3         = -2
    r_ref(3) = r(2) - r(4)                  !  2 - 4         = -2
    r_ref(4) = r(1) - r(2) + r(3) - r(4)    !  1 - 2 + 3 - 4 = -2
    
    CALL RFFTI(n, wsave)
    CALL RFFTF(n, r, wsave)

    ! Compare against reference
    IF (ABS(r(1) - r_ref(1)) > tolerance) THEN
      PRINT *, "RFFTF: Test 2 failed for N = 4, element 1"
      STOP
    END IF
    IF (ABS(r(2) - r_ref(2)) > tolerance) THEN
      PRINT *, "RFFTF: Test 2 failed for N = 4, element 2"
      STOP
    END IF
    IF (ABS(r(4) - r_ref(4)) > tolerance) THEN
      PRINT *, "RFFTF: Test 2 failed for N = 4, element 4"
      STOP
    END IF
    
    DEALLOCATE(r, wsave, r_ref)

  END SUBROUTINE test_rfftf

  SUBROUTINE test_rfftb()
    implicit none
    integer :: n
    real, allocatable :: r(:), wsave(:)
    real :: expected_rfftb(4), output_rfftb(4)

    n = 4
    allocate(r(n), wsave(2*n+15))
    r = (/1.0, 2.0, 3.0, 4.0/)
    expected_rfftb = (/4.0, 8.0, 12.0, 16.0/)

    call rffti(n, wsave)
    call rfftf(n, r, wsave)
    call rfftb(n, r, wsave)

    output_rfftb = r
    if (any(abs(output_rfftb - expected_rfftb) > tolerance)) then
      print *, "Expected: ", expected_rfftb
      print *, "Output: ", output_rfftb
      stop 61
    endif

    deallocate(r, wsave)
  END SUBROUTINE test_rfftb

  !> rffti1 DO 106 loop: N=8 (factors 4,2) causes 2 to be bubble-swapped to front of IFAC
  SUBROUTINE test_rffti1_do106()
    implicit none
    integer :: n, i
    real, allocatable :: r(:), wsave(:)

    n = 8
    allocate(r(n), wsave(3*n+15))
    call rffti(n, wsave)
    ! IFAC stored starting at wsave(2*n+1); verify 2 was moved to front: [8,2,2,4]
    if (nint(wsave(2*n+1)) /= 8 .or. nint(wsave(2*n+2)) /= 2 .or. &
        nint(wsave(2*n+3)) /= 2 .or. nint(wsave(2*n+4)) /= 4) stop 143
    do i=1,n; r(i)=real(i); end do
    call rfftf(n, r, wsave)
    if (any(abs(r - (/36.0,-4.0,9.65685463,-4.0,4.0,-4.0,1.65685415,-4.0/)) > tolerance)) stop 144
    call rfftb(n, r, wsave)
    do i=1,n
      if (abs(r(i) - real(i)*8.0) > tolerance) stop 145
    end do
    deallocate(r, wsave)
  END SUBROUTINE test_rffti1_do106

  SUBROUTINE test_rfftf1()
    implicit none
    integer :: n
    real, allocatable :: c(:), ch(:), wsave(:)
    real :: expected_rfftf1(4), output_rfftf1(4)

    n = 4
    allocate(c(n), ch(n), wsave(2*n+15))
    c = (/1.0, 2.0, 3.0, 4.0/)
    ch = 0.0
    expected_rfftf1 = (/10.0, -2.0, 2.0, -2.0/)

    call rffti(n, wsave)
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))

    output_rfftf1 = c
    if (any(abs(output_rfftf1 - expected_rfftf1) > tolerance)) then
      print *, "Expected: ", expected_rfftf1
      print *, "Output: ", output_rfftf1
      stop 62
    endif

    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftf1

  !> rfftf1 label 103: RADF2(IDO,L1,CH,C,WA) — the NA!=0 path for IP=2
  !> n=14=2*7: K1=1 IP=7 IDO=1 -> label 109 (NA=0), K1=2 IP=2 IDO=7 NA=1 -> label 103
  SUBROUTINE test_rfftf1_radf2_ch_c()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 14
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (abs(c(1) - 105.0) > tolerance) stop 146
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))
    do i=1,n
      if (abs(c(i) - real(i)*14.0) > 1e-3) stop 147
    end do
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftf1_radf2_ch_c

  !> rfftf1 label 107: RADF5(IDO,L1,CH,C,WA) — the NA!=0 path for IP=5
  !> n=25=5^2: K1=1 IP=5 IDO=1 NA=0->RADF5(C,CH), K1=2 IP=5 IDO=5 NA=1->label 107
  SUBROUTINE test_rfftf1_radf5_ch_c()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 25
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (abs(c(1) - 325.0) > tolerance) stop 148
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))
    do i=1,n
      if (abs(c(i) - real(i)*25.0) > 1e-1) stop 149
    end do
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftf1_radf5_ch_c

  !> rfftf1 label 108 fall-through: RADFG(C,C,C,CH,CH) + NA=1
  !> Requires crafted IFAC=[56,3,7,4,2] so rfftf1 processes IP=2,IP=4,IP=7 in that order;
  !> K1=3 has IP=7, IDO=8, NA=0 (after flips) -> falls through to RADFG(C,C,C,CH,CH), NA=1
  SUBROUTINE test_rfftf1_radfg_c_ch()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 56
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    ! Standard IFAC=[56,3,2,4,7]; reorder to [56,3,7,4,2] so K1=3 is IP=7,IDO=8,NA=0
    wsave(2*n+1)=56.0; wsave(2*n+2)=3.0; wsave(2*n+3)=7.0
    wsave(2*n+4)=4.0;  wsave(2*n+5)=2.0
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    ! NA=1 on return -> no CH->C copy; c(1) reflects RADFG(C,C,C,CH,CH) output
    if (abs(c(1) - 1596.0) > 1.0) stop 150
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftf1_radfg_c_ch

  SUBROUTINE test_rfftb1()
    implicit none
    integer :: n
    real, allocatable :: c(:), ch(:), wsave(:)
    real :: expected_rfftb1(4), output_rfftb1(4)

    n = 4
    allocate(c(n), ch(n), wsave(2*n+15))
    c = (/1.0, 2.0, 3.0, 4.0/)
    ch = 0.0
    expected_rfftb1 = (/4.0, 8.0, 12.0, 16.0/)

    call rffti(n, wsave)
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))

    output_rfftb1 = c
    if (any(abs(output_rfftb1 - expected_rfftb1) > tolerance)) then
      print *, "Expected: ", expected_rfftb1
      print *, "Output: ", output_rfftb1
      stop 63
    endif

    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftb1

  SUBROUTINE test_radb2()
    implicit none
    real :: cc(2,2,1), ch(2,1,2), wa1(2)
    real :: expected_radb2(4), output_radb2(4)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    cc(1,1,1) = 1.0
    cc(2,1,1) = 2.0
    cc(1,2,1) = 3.0
    cc(2,2,1) = 4.0

    call radb2(2, 1, cc, ch, wa1)

    expected_radb2 = (/5.0, 4.0, -3.0, -6.0/)
    output_radb2 = (/ch(1,1,1), ch(2,1,1), ch(1,1,2), ch(2,1,2)/)

    if (any(abs(output_radb2 - expected_radb2) > tolerance)) then
      print *, "Expected: ", expected_radb2
      print *, "Output: ", output_radb2
      stop 64
    endif
  END SUBROUTINE test_radb2

  SUBROUTINE test_radb3()
    implicit none
    real :: cc(1,3,1), ch(1,1,3), wa1(1), wa2(1)
    real :: expected_radb3(3), output_radb3(3)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    wa2 = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0
    cc(1,3,1) = 3.0

    call radb3(1, 1, cc, ch, wa1, wa2)

    expected_radb3 = (/5.0, -6.196152, 4.196152/)
    output_radb3 = (/ch(1,1,1), ch(1,1,2), ch(1,1,3)/)

    if (any(abs(output_radb3 - expected_radb3) > tolerance)) then
      print *, "Expected: ", expected_radb3
      print *, "Output: ", output_radb3
      stop 65
    endif
  END SUBROUTINE test_radb3

  SUBROUTINE test_radb5()
    implicit none
    real :: cc(1,5,1), ch(1,1,5), wa1(1), wa2(1), wa3(1), wa4(1)
    real :: expected_radb5(5), output_radb5(5)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    wa2 = 0.0
    wa3 = 0.0
    wa4 = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0
    cc(1,3,1) = 3.0
    cc(1,4,1) = 4.0
    cc(1,5,1) = 5.0

    call radb5(1, 1, cc, ch, wa1, wa2, wa3, wa4)

    expected_radb5 = (/13.0, -15.8202600, 6.21992207, -5.74778652, 7.34812450/)
    output_radb5 = (/ch(1,1,1), ch(1,1,2), ch(1,1,3), ch(1,1,4), ch(1,1,5)/)

    if (any(abs(output_radb5 - expected_radb5) > tolerance)) then
      print *, "Expected: ", expected_radb5
      print *, "Output: ", output_radb5
      stop 66
    endif
  END SUBROUTINE test_radb5

  SUBROUTINE test_radb2_ido1()
    implicit none
    real :: cc(1,2,1), ch(1,1,2), wa1(1)
    real :: expected(2), output(2)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0

    call radb2(1, 1, cc, ch, wa1)

    expected = (/3.0, -1.0/)
    output = (/ch(1,1,1), ch(1,1,2)/)

    if (any(abs(output - expected) > tolerance)) then
      print *, "Expected: ", expected
      print *, "Output: ", output
      stop 69
    endif
  END SUBROUTINE test_radb2_ido1

  SUBROUTINE test_radb3_ido1()
    implicit none
    real :: cc(1,3,1), ch(1,1,3), wa1(1), wa2(1)
    real :: expected(3), output(3)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    wa2 = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0
    cc(1,3,1) = 3.0

    call radb3(1, 1, cc, ch, wa1, wa2)

    expected = (/5.0, -6.19615221, 4.19615221/)
    output = (/ch(1,1,1), ch(1,1,2), ch(1,1,3)/)

    if (any(abs(output - expected) > tolerance)) then
      print *, "Expected: ", expected
      print *, "Output: ", output
      stop 70
    endif
  END SUBROUTINE test_radb3_ido1

  SUBROUTINE test_radb5_ido1()
    implicit none
    real :: cc(1,5,1), ch(1,1,5), wa1(1), wa2(1), wa3(1), wa4(1)
    real :: expected(5), output(5)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    wa2 = 0.0
    wa3 = 0.0
    wa4 = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0
    cc(1,3,1) = 3.0
    cc(1,4,1) = 4.0
    cc(1,5,1) = 5.0

    call radb5(1, 1, cc, ch, wa1, wa2, wa3, wa4)

    expected = (/13.0, -15.8202600, 6.21992207, -5.74778652, 7.34812450/)
    output = (/ch(1,1,1), ch(1,1,2), ch(1,1,3), ch(1,1,4), ch(1,1,5)/)

    if (any(abs(output - expected) > tolerance)) then
      print *, "Expected: ", expected
      print *, "Output: ", output
      stop 71
    endif
  END SUBROUTINE test_radb5_ido1

  SUBROUTINE test_radf5_ido1()
    implicit none
    real :: cc(1,1,5), ch(1,5,1), wa1(1), wa2(1), wa3(1), wa4(1)
    real :: expected(5), output(5)

    cc = 0.0
    ch = 0.0
    wa1 = 0.0
    wa2 = 0.0
    wa3 = 0.0
    wa4 = 0.0
    cc(1,1,1) = 1.0
    cc(1,1,2) = 2.0
    cc(1,1,3) = 3.0
    cc(1,1,4) = 4.0
    cc(1,1,5) = 5.0

    call radf5(1, 1, cc, ch, wa1, wa2, wa3, wa4)

    expected = (/15.0, -2.49999976, 3.44095492, -2.49999976, 0.81229919/)
    output = (/ch(1,1,1), ch(1,2,1), ch(1,3,1), ch(1,4,1), ch(1,5,1)/)

    if (any(abs(output - expected) > tolerance)) then
      print *, "Expected: ", expected
      print *, "Output: ", output
      stop 72
    endif
  END SUBROUTINE test_radf5_ido1

  SUBROUTINE test_radfg_ido1()
    implicit none
    real :: cc(1,2,1), c1(1,1,2), c2(1,2), ch(1,1,2), ch2(1,2), wa(1)
    real :: expected(4), output(4)

    cc = 0.0
    c1 = 0.0
    c2 = 0.0
    ch = 0.0
    ch2 = 0.0
    wa = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0

    call radfg(1, 2, 1, 1, cc, c1, c2, ch, ch2, wa)

    expected = (/0.0, 0.0, 0.0, 0.0/)
    output = (/ch(1,1,1), ch(1,1,2), ch2(1,1), ch2(1,2)/)

    if (any(abs(output - expected) > tolerance)) then
      print *, "Expected: ", expected
      print *, "Output: ", output
      stop 73
    endif
  END SUBROUTINE test_radfg_ido1

  SUBROUTINE test_radbg_ido1()
    implicit none
    real :: cc(1,2,1), c1(1,1,2), c2(1,2), ch(1,1,2), ch2(1,2), wa(1)
    real :: expected(4), output(4)

    cc = 0.0
    c1 = 0.0
    c2 = 0.0
    ch = 0.0
    ch2 = 0.0
    wa = 0.0
    cc(1,1,1) = 1.0
    cc(1,2,1) = 2.0

    call radbg(1, 2, 1, 1, cc, c1, c2, ch, ch2, wa)

    expected = (/1.0, 0.0, 0.0, 0.0/)
    output = (/ch(1,1,1), ch(1,1,2), ch2(1,1), ch2(1,2)/)

    if (any(abs(output - expected) > tolerance)) then
      print *, "Expected: ", expected
      print *, "Output: ", output
      stop 74
    endif
  END SUBROUTINE test_radbg_ido1

  SUBROUTINE test_radfg()
    implicit none
    real :: cc(2,2,1), c1(2,1,2), c2(1,2), ch(2,1,2), ch2(1,2), wa(2)
    real :: expected_radfg(4), output_radfg(4)

    cc = 0.0
    c1 = 0.0
    c2 = 0.0
    ch = 0.0
    ch2 = 0.0
    wa = 0.0
    cc(1,1,1) = 1.0
    cc(2,1,1) = 2.0
    cc(1,2,1) = 3.0
    cc(2,2,1) = 4.0

    call radfg(2, 2, 1, 1, cc, c1, c2, ch, ch2, wa)

    expected_radfg = (/0.0, 0.0, 0.0, 0.0/)
    output_radfg = (/ch(1,1,1), ch(2,1,1), ch(1,1,2), ch(2,1,2)/)

    if (any(abs(output_radfg - expected_radfg) > tolerance)) then
      print *, "Expected: ", expected_radfg
      print *, "Output: ", output_radfg
      stop 67
    endif
  END SUBROUTINE test_radfg

  SUBROUTINE test_radbg()
    implicit none
    real :: cc(2,2,1), c1(2,1,2), c2(1,2), ch(2,1,2), ch2(1,2), wa(2)
    real :: expected_radbg(4), output_radbg(4)

    cc = 0.0
    c1 = 0.0
    c2 = 0.0
    ch = 0.0
    ch2 = 0.0
    wa = 0.0
    cc(1,1,1) = 4.0
    cc(2,1,1) = 2.0
    cc(1,2,1) = -4.0
    cc(2,2,1) = -2.0

    call radbg(2, 2, 1, 1, cc, c1, c2, ch, ch2, wa)

    expected_radbg = (/4.0, 2.0, 0.0, 0.0/)
    output_radbg = (/ch(1,1,1), ch(2,1,1), ch(1,1,2), ch(2,1,2)/)

    if (any(abs(output_radbg - expected_radbg) > tolerance)) then
      print *, "Expected: ", expected_radbg
      print *, "Output: ", output_radbg
      stop 68
    endif
  END SUBROUTINE test_radbg

  !> RADB4 with IDO=1: exercises only the DO 101 K loop and the IDO.LT.2 branch
  SUBROUTINE test_radb4_ido1()
    implicit none
    real :: cc(1,4,1), ch(1,1,4), wa1(1), wa2(1), wa3(1)

    cc=0.0; ch=0.0; wa1=0.0; wa2=0.0; wa3=0.0
    cc(1,1,1)=1.0; cc(1,2,1)=2.0; cc(1,3,1)=3.0; cc(1,4,1)=4.0

    call radb4(1,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(1,1,2),ch(1,1,3),ch(1,1,4)/) &
        - (/9.0,-9.0,1.0,3.0/)) > tolerance)) stop 105
  END SUBROUTINE test_radb4_ido1

  !> RADB4 with IDO=2: exercises the IDO.EQ.2 branch (labels 105/106, even-end)
  SUBROUTINE test_radb4_ido2()
    implicit none
    real :: cc(2,4,1), ch(2,1,4), wa1(2), wa2(2), wa3(2)

    cc=0.0; ch=0.0; wa1=0.0; wa2=0.0; wa3=0.0
    cc(1,1,1)=1.0; cc(2,1,1)=2.0
    cc(1,2,1)=3.0; cc(2,2,1)=4.0
    cc(1,3,1)=5.0; cc(2,3,1)=6.0
    cc(1,4,1)=7.0; cc(2,4,1)=8.0

    call radb4(2,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(2,1,1)/) - (/17.0,16.0/)) > tolerance)) stop 106
    if (any(abs((/ch(1,1,2),ch(2,1,2)/) - (/-17.0,-19.7989902/)) > tolerance)) stop 107
    if (any(abs((/ch(1,1,3),ch(2,1,3)/) - (/1.0,8.0/)) > tolerance)) stop 108
    if (any(abs((/ch(1,1,4),ch(2,1,4)/) - (/3.0,-8.48528099/)) > tolerance)) stop 109
  END SUBROUTINE test_radb4_ido2

  !> RADB4 with IDO=3: exercises the !OCL NOVREC loop; odd IDO skips even-end branch
  SUBROUTINE test_radb4_ido3()
    implicit none
    real :: cc(3,4,1), ch(3,1,4), wa1(3), wa2(3), wa3(3)
    integer :: i

    cc=0.0; ch=0.0
    do i=1,12
      cc(mod(i-1,3)+1, int(real(i - 1) / 3.0) + 1, 1) = real(i)
    end do
    wa1(1)=1.0; wa1(2)=0.866025; wa1(3)=0.5
    wa2(1)=1.0; wa2(2)=0.5;      wa2(3)=0.866025
    wa3(1)=1.0; wa3(2)=0.0;      wa3(3)=1.0

    call radb4(3,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/25.0,24.0,-4.0/)) > tolerance)) stop 110
    if (any(abs((/ch(1,1,2),ch(2,1,2),ch(3,1,2)/) - (/-25.0,-37.5884476,-1.05254936/)) > tolerance)) stop 111
    if (any(abs((/ch(1,1,3),ch(2,1,3),ch(3,1,3)/) - (/1.0,6.0,-12.0/)) > tolerance)) stop 112
    if (any(abs((/ch(1,1,4),ch(2,1,4),ch(3,1,4)/) - (/3.0,6.0,10.0/)) > tolerance)) stop 113
  END SUBROUTINE test_radb4_ido3

  !> RADB4 with IDO=4: exercises !OCL NOVREC loop and the even-end branch (MOD(IDO,2)=0)
  SUBROUTINE test_radb4_ido4()
    implicit none
    real :: cc(4,4,1), ch(4,1,4), wa1(4), wa2(4), wa3(4)
    integer :: i

    cc=0.0; ch=0.0
    do i=1,16
      cc(mod(i-1,4)+1, int(real(i-1.)/4.)+1, 1) = real(i)
    end do
    wa1(1)=1.0; wa1(2)=0.0; wa1(3)=0.0; wa1(4)=1.0
    wa2(1)=1.0; wa2(2)=0.0; wa2(3)=0.0; wa2(4)=1.0
    wa3(1)=1.0; wa3(2)=0.0; wa3(3)=0.0; wa3(4)=1.0

    call radb4(4,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1),ch(4,1,1)/) &
        - (/33.0,32.0,-8.0,32.0/)) > tolerance)) stop 114
    if (any(abs((/ch(1,1,2),ch(2,1,2),ch(3,1,2),ch(4,1,2)/) &
        - (/-33.0,-30.0,22.0,-36.7695503/)) > tolerance)) stop 115
    if (any(abs((/ch(1,1,3),ch(2,1,3),ch(3,1,3),ch(4,1,3)/) &
        - (/1.0,0.0,-16.0,16.0/)) > tolerance)) stop 116
    if (any(abs((/ch(1,1,4),ch(2,1,4),ch(3,1,4),ch(4,1,4)/) &
        - (/3.0,6.0,14.0,-14.1421356/)) > tolerance)) stop 117
  END SUBROUTINE test_radb4_ido4

  !> RADF4 with IDO=1: exercises only the DO 101 K loop and the IDO.LT.2 branch
  SUBROUTINE test_radf4_ido1()
    implicit none
    real :: cc(1,1,4), ch(1,4,1), wa1(1), wa2(1), wa3(1)

    cc=0.0; ch=0.0; wa1=0.0; wa2=0.0; wa3=0.0
    cc(1,1,1)=1.0; cc(1,1,2)=2.0; cc(1,1,3)=3.0; cc(1,1,4)=4.0

    call radf4(1,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(1,2,1),ch(1,3,1),ch(1,4,1)/) &
        - (/10.0,-2.0,2.0,-2.0/)) > tolerance)) stop 118
  END SUBROUTINE test_radf4_ido1

  !> RADF4 with IDO=2: exercises the IDO.EQ.2 branch (labels 105/106, even-end)
  SUBROUTINE test_radf4_ido2()
    implicit none
    real :: cc(2,1,4), ch(2,4,1), wa1(2), wa2(2), wa3(2)
    integer :: i

    cc=0.0; ch=0.0; wa1=0.0; wa2=0.0; wa3=0.0
    do i=1,8
      cc(mod(i-1,2)+1, 1, int(real(i-1.)/2.)+1) = real(i)
    end do

    call radf4(2,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(2,1,1)/) - (/16.0,-0.828427076/)) > tolerance)) stop 119
    if (any(abs((/ch(1,2,1),ch(2,2,1)/) - (/-14.4852810,-4.0/)) > tolerance)) stop 120
    if (any(abs((/ch(1,3,1),ch(2,3,1)/) - (/4.0,4.82842731/)) > tolerance)) stop 121
    if (any(abs((/ch(1,4,1),ch(2,4,1)/) - (/-2.48528099,-4.0/)) > tolerance)) stop 122
  END SUBROUTINE test_radf4_ido2

  !> RADF4 with IDO=3: exercises the !OCL NOVREC loop; odd IDO skips even-end branch
  SUBROUTINE test_radf4_ido3()
    implicit none
    real :: cc(3,1,4), ch(3,4,1), wa1(3), wa2(3), wa3(3)
    integer :: i

    cc=0.0; ch=0.0
    do i=1,12
      cc(mod(i-1,3)+1, 1, int(real(i-1.)/3.)+1) = real(i)
    end do
    wa1(1)=1.0; wa1(2)=0.866025; wa1(3)=0.5
    wa2(1)=1.0; wa2(2)=0.5;      wa2(3)=0.866025
    wa3(1)=1.0; wa3(2)=0.0;      wa3(3)=1.0

    call radf4(3,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/22.0,35.6961517,21.6698761/)) > tolerance)) stop 123
    if (any(abs((/ch(1,2,1),ch(2,2,1),ch(3,2,1)/) - (/-0.169875145,2.80385017,-6.0/)) > tolerance)) stop 124
    if (any(abs((/ch(1,3,1),ch(2,3,1),ch(3,3,1)/) - (/6.0,-20.8301239,-1.19614983/)) > tolerance)) stop 125
    if (any(abs((/ch(1,4,1),ch(2,4,1),ch(3,4,1)/) - (/-6.69614983,5.66987514,-6.0/)) > tolerance)) stop 126
  END SUBROUTINE test_radf4_ido3

  !> RADF4 with IDO=4: exercises !OCL NOVREC loop and the even-end branch (MOD(IDO,2)=0)
  SUBROUTINE test_radf4_ido4()
    implicit none
    real :: cc(4,1,4), ch(4,4,1), wa1(4), wa2(4), wa3(4)
    integer :: i

    cc=0.0; ch=0.0
    do i=1,16
      cc(mod(i-1,4)+1, 1, int(real(i-1.)/4.)+1) = real(i)
    end do
    wa1(1)=1.0; wa1(2)=0.0; wa1(3)=0.0; wa1(4)=1.0
    wa2(1)=1.0; wa2(2)=0.0; wa2(3)=0.0; wa2(4)=1.0
    wa3(1)=1.0; wa3(2)=0.0; wa3(3)=0.0; wa3(4)=1.0

    call radf4(4,1,cc,ch,wa1,wa2,wa3)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1),ch(4,1,1)/) &
        - (/28.0,32.0,36.0,-1.65685415/)) > tolerance)) stop 127
    if (any(abs((/ch(1,2,1),ch(2,2,1),ch(3,2,1),ch(4,2,1)/) &
        - (/-28.9705620,0.0,16.0,-8.0/)) > tolerance)) stop 128
    if (any(abs((/ch(1,3,1),ch(2,3,1),ch(3,3,1),ch(4,3,1)/) &
        - (/8.0,-16.0,0.0,9.65685463/)) > tolerance)) stop 129
    if (any(abs((/ch(1,4,1),ch(2,4,1),ch(3,4,1),ch(4,4,1)/) &
        - (/-4.97056198,-8.0,8.0,-8.0/)) > tolerance)) stop 130
  END SUBROUTINE test_radf4_ido4

  !> rfftf1/rfftb1 n=16 (4^2): RADF4/RADB4 with IDO=4 even-end branch via round-trip
  SUBROUTINE test_rfft1_n16()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 16
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n,c,ch,wsave(n+1),wsave(2*n+1))
    if (any(abs(c(1:4) - (/136.0,-8.0,40.2187157,-8.0/)) > tolerance)) stop 131
    call rfftb1(n,c,ch,wsave(n+1),wsave(2*n+1))
    do i=1,n
      if (abs(c(i) - real(i)*16.0) > 1e-3) stop 132
    end do
    deallocate(c,ch,wsave)
  END SUBROUTINE test_rfft1_n16

  !> rfftf1/rfftb1 n=20 (4*5): RADF4/RADB4 with IDO=5 (odd, no even-end) via round-trip
  SUBROUTINE test_rfft1_n20()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 20
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n,c,ch,wsave(n+1),wsave(2*n+1))
    if (any(abs(c(1:4) - (/210.0,-10.0,63.1375160,-10.0/)) > tolerance)) stop 133
    call rfftb1(n,c,ch,wsave(n+1),wsave(2*n+1))
    do i=1,n
      if (abs(c(i) - real(i)*20.0) > 1e-2) stop 134
    end do
    deallocate(c,ch,wsave)
  END SUBROUTINE test_rfft1_n20

  !> rfftf1/rfftb1 n=49 (7^2): RADFG with IDO=7, L1=1 (NBD>L1, labels 107/114/138)
  SUBROUTINE test_rfft1_n49()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 49
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n,c,ch,wsave(n+1),wsave(2*n+1))
    if (abs(c(1) - 1225.0) > tolerance) stop 135
    call rfftb1(n,c,ch,wsave(n+1),wsave(2*n+1))
    do i=1,n
      if (abs(c(i) - real(i)*49.0) > 1e-1) stop 136
    end do
    deallocate(c,ch,wsave)
  END SUBROUTINE test_rfft1_n49

  !> rfftf1/rfftb1 n=686 (2*7^3): covers RADFG labels 104/115/132/141 (IDO<L1,NBD<L1)
  !> and 107/114/138 (NBD>L1) via the three RADFG calls at different pass depths
  SUBROUTINE test_rfft1_n686()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 686
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i=1,n; c(i)=real(i); end do
    ch=0.0
    call rfftf1(n,c,ch,wsave(n+1),wsave(2*n+1))
    if (abs(c(1) - 235641.0) > 1.0) stop 137
    call rfftb1(n,c,ch,wsave(n+1),wsave(2*n+1))
    do i=1,4
      if (abs(c(i) - real(i)*686.0) > 1.0) stop 138
    end do
    deallocate(c,ch,wsave)
  END SUBROUTINE test_rfft1_n686

  !> rfftb1/rfftf1 with n=6 (radix 2*3): exercises NA!=0 branches in both routines
  SUBROUTINE test_rfftb1_n6()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 6
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i = 1, n
      c(i) = real(i)
    end do
    ch = 0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (any(abs(c - (/21.0, -3.0, 5.19615269, -3.0, 1.73205101, -3.0/)) > tolerance)) stop 75
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (any(abs(c - (/6.0, 12.0, 18.0, 24.0, 30.0, 36.0/)) > tolerance)) stop 76
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftb1_n6

  !> rfftb1/rfftf1 with n=7 (prime): exercises RADFG/RADBG general-radix branches
  SUBROUTINE test_rfftb1_n7()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 7
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i = 1, n
      c(i) = real(i)
    end do
    ch = 0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (any(abs(c - (/28.0, -3.5, 7.26782513, -3.5, 2.79115677, -3.50000048, 0.798852086/)) > tolerance)) stop 77
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))
    ! result should be 7 * original
    if (any(abs(c - (/7.0, 14.0, 21.0, 28.0, 35.0, 42.0, 49.0/)) > 1e-3)) stop 78
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftb1_n7

  !> rfftb1/rfftf1 with n=10 (2*5): exercises RADF5/RADB5 with IDO>1 (NA!=0 branches)
  SUBROUTINE test_rfftb1_n10()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 10
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i = 1, n
      c(i) = real(i)
    end do
    ch = 0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (any(abs(c - (/55.0, -5.0, 15.3884182, -5.0, 6.88190937, -5.0, 3.63271236, &
                      -5.00000048, 1.62459803, -5.0/)) > tolerance)) stop 79
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (any(abs(c - (/10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0/)) > 1e-3)) stop 80
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftb1_n10

  !> rfftb1/rfftf1 with n=15 (3*5): exercises RADB3/RADF3/RADB5/RADF5 with IDO>1
  SUBROUTINE test_rfftb1_n15()
    implicit none
    integer :: n, i
    real, allocatable :: c(:), ch(:), wsave(:)

    n = 15
    allocate(c(n), ch(n), wsave(3*n+15))
    call rffti(n, wsave)
    do i = 1, n
      c(i) = real(i)
    end do
    ch = 0.0
    call rfftf1(n, c, ch, wsave(n+1), wsave(2*n+1))
    if (any(abs(c - (/120.0, -7.5, 35.2847252, -7.5, 16.8452740, -7.5, 10.3228636, &
                      -7.5, 6.75302982, -7.5, 4.33012676, -7.5, 2.43689823, &
                      -7.5, 0.788281441/)) > tolerance)) stop 81
    call rfftb1(n, c, ch, wsave(n+1), wsave(2*n+1))
    do i = 1, n
      if (abs(c(i) - real(i)*15.0) > 1e-3) stop 82
    end do
    deallocate(c, ch, wsave)
  END SUBROUTINE test_rfftb1_n15

  !> RADB2 with IDO=3: exercises the !OCL NOVREC inner loop (IDO odd, no even-end branch)
  SUBROUTINE test_radb2_ido3()
    implicit none
    real :: cc(3,2,1), ch(3,1,2), wa1(3)

    cc = 0.0
    ch = 0.0
    cc(1,1,1)=1.0; cc(2,1,1)=2.0; cc(3,1,1)=3.0
    cc(1,2,1)=4.0; cc(2,2,1)=5.0; cc(3,2,1)=6.0
    wa1(1)=1.0; wa1(2)=0.5; wa1(3)=0.866025

    call radb2(3, 1, cc, ch, wa1)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/7.0, 6.0, -2.0/)) > tolerance)) stop 83
    if (any(abs((/ch(1,1,2),ch(2,1,2),ch(3,1,2)/) - (/-5.0, -6.0, 7.0/)) > tolerance)) stop 84
  END SUBROUTINE test_radb2_ido3

  !> RADB2 with IDO=4: exercises !OCL NOVREC loop + even IDO end branch
  SUBROUTINE test_radb2_ido4()
    implicit none
    real :: cc(4,2,1), ch(4,1,2), wa1(4)

    cc = 0.0
    ch = 0.0
    cc(1,1,1)=1.0; cc(2,1,1)=2.0; cc(3,1,1)=3.0; cc(4,1,1)=4.0
    cc(1,2,1)=5.0; cc(2,2,1)=6.0; cc(3,2,1)=7.0; cc(4,2,1)=8.0
    wa1(1)=1.0; wa1(2)=0.0; wa1(3)=0.0; wa1(4)=1.0

    call radb2(4, 1, cc, ch, wa1)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1),ch(4,1,1)/) - (/9.0, 8.0, -4.0, 8.0/)) > tolerance)) stop 85
    if (any(abs((/ch(1,1,2),ch(2,1,2),ch(3,1,2),ch(4,1,2)/) - (/-7.0, -4.0, 10.0, -10.0/)) > tolerance)) stop 86
  END SUBROUTINE test_radb2_ido4

  !> RADB3 with IDO=3: exercises the !OCL NOVREC inner loop
  SUBROUTINE test_radb3_ido3()
    implicit none
    real :: cc(3,3,1), ch(3,1,3), wa1(3), wa2(3)

    cc = 0.0
    ch = 0.0
    cc(1,1,1)=1.0; cc(2,1,1)=2.0; cc(3,1,1)=3.0
    cc(1,2,1)=4.0; cc(2,2,1)=5.0; cc(3,2,1)=6.0
    cc(1,3,1)=7.0; cc(2,3,1)=8.0; cc(3,3,1)=9.0
    wa1(1)=1.0; wa1(2)=0.866025; wa1(3)=0.5
    wa2(1)=1.0; wa2(2)=0.5; wa2(3)=0.866025

    call radb3(3, 1, cc, ch, wa1, wa2)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/13.0, 14.0, 7.0/)) > tolerance)) stop 87
    if (any(abs((/ch(1,1,2),ch(2,1,2),ch(3,1,2)/) - (/-17.1243553, -19.9903793, -9.49999237/)) > tolerance)) stop 88
    if (any(abs((/ch(1,1,3),ch(2,1,3),ch(3,1,3)/) - (/7.12435532, 9.35640621, 1.59807611/)) > tolerance)) stop 89
  END SUBROUTINE test_radb3_ido3

  !> RADB5 with IDO=3: exercises the loop after "IF (IDO .EQ. 1) RETURN"
  SUBROUTINE test_radb5_ido3()
    implicit none
    real :: cc(3,5,1), ch(3,1,5), wa1(3), wa2(3), wa3(3), wa4(3)
    integer :: i

    cc = 0.0
    ch = 0.0
    do i = 1, 15
      cc((mod(i-1,3))+1, int(real(i-1.)/3.)+1, 1) = real(i)
    end do
    wa1(1)=1.0; wa1(2)=0.809017; wa1(3)=0.587785
    wa2(1)=1.0; wa2(2)=0.309017; wa2(3)=0.951057
    wa3(1)=1.0; wa3(2)=-0.309017; wa3(3)=0.951057
    wa4(1)=1.0; wa4(2)=-0.809017; wa4(3)=0.587785

    call radb5(3, 1, cc, ch, wa1, wa2, wa3, wa4)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/37.0, 38.0, 11.0/)) > tolerance)) stop 90
    if (any(abs((/ch(1,1,2),ch(2,1,2),ch(3,1,2)/) - (/-43.3054123, -48.0942268, -27.0704327/)) > tolerance)) stop 91
    if (any(abs((/ch(1,1,3),ch(2,1,3),ch(3,1,3)/) - (/15.2066822, 16.3466930, 4.55505562/)) > tolerance)) stop 92
    if (any(abs((/ch(1,1,4),ch(2,1,4),ch(3,1,4)/) - (/-17.7902737, -16.0322285, 7.64156532/)) > tolerance)) stop 93
    if (any(abs((/ch(1,1,5),ch(2,1,5),ch(3,1,5)/) - (/13.8890038, 10.7182236, -17.2008247/)) > tolerance)) stop 94
  END SUBROUTINE test_radb5_ido3

  !> RADF5 with IDO=3: exercises the loop after "IF (IDO .EQ. 1) RETURN"
  SUBROUTINE test_radf5_ido3()
    implicit none
    real :: cc(3,1,5), ch(3,5,1), wa1(3), wa2(3), wa3(3), wa4(3)
    integer :: i

    cc = 0.0
    ch = 0.0
    do i = 1, 15
      cc((mod(i-1,3))+1, 1, int(real(i-1.)/3.)+1) = real(i)
    end do
    wa1(1)=1.0; wa1(2)=0.809017; wa1(3)=0.587785
    wa2(1)=1.0; wa2(2)=0.309017; wa2(3)=0.951057
    wa3(1)=1.0; wa3(2)=-0.309017; wa3(3)=0.951057
    wa4(1)=1.0; wa4(2)=-0.809017; wa4(3)=0.587785

    call radf5(3, 1, cc, ch, wa1, wa2, wa3, wa4)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/35.0, 31.7917957, 53.2082024/)) > tolerance)) stop 95
    if (any(abs((/ch(1,2,1),ch(2,2,1),ch(3,2,1)/) - (/19.3929367, -3.64932156, -7.5/)) > tolerance)) stop 96
    if (any(abs((/ch(1,3,1),ch(2,3,1),ch(3,3,1)/) - (/10.3228645, -37.3929367, -15.6493235/)) > tolerance)) stop 97
    if (any(abs((/ch(1,4,1),ch(2,4,1),ch(3,4,1)/) - (/3.99207687, 11.7266502, -7.5/)) > tolerance)) stop 98
    if (any(abs((/ch(1,5,1),ch(2,5,1),ch(3,5,1)/) - (/2.43689752, -7.78387260, -14.4815521/)) > tolerance)) stop 99
  END SUBROUTINE test_radf5_ido3

  !> RADBG with IDO=3, IP=3, L1=1: exercises IDO>1, NBD>=L1 (inner K loop first)
  SUBROUTINE test_radbg_ido3_l1()
    implicit none
    real :: cc(3,3,1), c1(3,1,3), c2(3,3), ch(3,1,3), ch2(3,3), wa(10)
    integer :: i, k

    cc = 0.0; c1 = 0.0; c2 = 0.0; ch = 0.0; ch2 = 0.0; wa = 0.0
    do i = 1, 3
      do k = 1, 3
        cc(i,k,1) = real((k-1)*3 + i)
      end do
    end do
    wa(1)=1.0; wa(2)=0.5; wa(3)=0.866025
    wa(4)=0.5; wa(5)=0.866025

    call radbg(3, 3, 1, 3, cc, c1, c2, ch, ch2, wa)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/1.0, 2.0, 3.0/)) > tolerance)) stop 100
  END SUBROUTINE test_radbg_ido3_l1

  !> RADBG with IDO=3, IP=3, L1=2: NBD>=L1, exercises the NBD>=L1 inner loop (DO 111)
  SUBROUTINE test_radbg_ido3_l2()
    implicit none
    real :: cc(3,3,2), c1(3,2,3), c2(6,3), ch(3,2,3), ch2(6,3), wa(10)
    integer :: i, k

    cc = 0.0; c1 = 0.0; c2 = 0.0; ch = 0.0; ch2 = 0.0; wa = 0.0
    do i = 1, 3
      do k = 1, 3
        cc(i,k,1) = real((k-1)*3 + i)
        cc(i,k,2) = real((k-1)*3 + i) + 15.0
      end do
    end do
    wa(1)=1.0; wa(2)=0.5; wa(3)=0.866025
    wa(4)=0.5; wa(5)=0.866025

    call radbg(3, 3, 2, 6, cc, c1, c2, ch, ch2, wa)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/1.0, 2.0, 3.0/)) > tolerance)) stop 101
    if (any(abs((/ch(1,2,1),ch(2,2,1),ch(3,2,1)/) - (/16.0, 17.0, 18.0/)) > tolerance)) stop 102
  END SUBROUTINE test_radbg_ido3_l2

  !> RADBG with IDO=3, IP=3, L1=5 (L1>IDO): exercises the IDO<L1 branch (DO 105) and NBD<L1 paths
  SUBROUTINE test_radbg_ido3_l5()
    implicit none
    real :: cc(3,3,5), c1(3,5,3), c2(15,3), ch(3,5,3), ch2(15,3), wa(30)
    integer :: i, k, j

    cc = 0.0; c1 = 0.0; c2 = 0.0; ch = 0.0; ch2 = 0.0; wa = 0.0
    do i = 1, 3
      do k = 1, 5
        do j = 1, 3
          cc(i,j,k) = real(((j-1)*5 + k - 1)*3 + i)
        end do
      end do
    end do
    wa(1)=1.0; wa(2)=0.5; wa(3)=0.866025
    wa(4)=0.5; wa(5)=0.866025

    call radbg(3, 3, 5, 15, cc, c1, c2, ch, ch2, wa)

    if (any(abs((/ch(1,1,1),ch(2,1,1),ch(3,1,1)/) - (/1.0, 2.0, 3.0/)) > tolerance)) stop 103
    if (any(abs((/ch(1,5,1),ch(2,5,1),ch(3,5,1)/) - (/13.0, 14.0, 15.0/)) > tolerance)) stop 104
  END SUBROUTINE test_radbg_ido3_l5

  !> Unit test driver
  SUBROUTINE run_fftpack_tests()
    implicit none

    print *, "Running fftpack tests..."

    call test_dcrft()
    print *, "dcrft tests passed."

    call test_scrft()
    print *, "scrft tests passed."

    call test_csfft()
    print *, "csfft tests passed."

    call test_drcft()
    print *, "drcft tests passed."

    call test_srcft()
    print *, "srcft tests passed."

    call test_scfft()
    print *, "scfft tests passed."

    call test_rfftf()
    print *, "rfftf tests passed."

    call test_rfftb()
    print *, "rfftb tests passed."

    call test_rffti1_do106()
    print *, "rffti1 DO 106 (factor-2 reorder) tests passed."

    call test_rfftf1()
    print *, "rfftf1 tests passed."

    call test_rfftf1_radf2_ch_c()
    print *, "rfftf1 label 103 RADF2(CH,C) tests passed."

    call test_rfftf1_radf5_ch_c()
    print *, "rfftf1 label 107 RADF5(CH,C) tests passed."

    call test_rfftf1_radfg_c_ch()
    print *, "rfftf1 label 108 fall-through RADFG(C,C,C,CH,CH)+NA=1 tests passed."

    call test_rfftb1()
    print *, "rfftb1 tests passed."

    call test_radb4_ido1()
    print *, "radb4 ido=1 branch tests passed."

    call test_radb4_ido2()
    print *, "radb4 ido=2 (even-end) tests passed."

    call test_radb4_ido3()
    print *, "radb4 ido=3 (OCL NOVREC loop, odd) tests passed."

    call test_radb4_ido4()
    print *, "radb4 ido=4 (OCL NOVREC loop + even-end) tests passed."

    call test_radf4_ido1()
    print *, "radf4 ido=1 branch tests passed."

    call test_radf4_ido2()
    print *, "radf4 ido=2 (even-end) tests passed."

    call test_radf4_ido3()
    print *, "radf4 ido=3 (OCL NOVREC loop, odd) tests passed."

    call test_radf4_ido4()
    print *, "radf4 ido=4 (OCL NOVREC loop + even-end) tests passed."

    call test_rfft1_n16()
    print *, "rfftf1/rfftb1 n=16 (radf4/radb4 IDO=4 even-end) tests passed."

    call test_rfft1_n20()
    print *, "rfftf1/rfftb1 n=20 (radf4/radb4 IDO=5 odd) tests passed."

    call test_rfft1_n49()
    print *, "rfftf1/rfftb1 n=49 (radfg IDO=7>L1=1, labels 107/114/138) tests passed."

    call test_rfft1_n686()
    print *, "rfftf1/rfftb1 n=686 (radfg IDO<L1, NBD<L1, NBD>L1 all branches) tests passed."

    call test_rfftb1_n6()
    print *, "rfftb1/rfftf1 n=6 (radix 2+3) tests passed."

    call test_rfftb1_n7()
    print *, "rfftb1/rfftf1 n=7 (prime, general radix) tests passed."

    call test_rfftb1_n10()
    print *, "rfftb1/rfftf1 n=10 (2*5, radf5/radb5 IDO>1) tests passed."

    call test_rfftb1_n15()
    print *, "rfftb1/rfftf1 n=15 (3*5, radb3/radb5 IDO>1) tests passed."

    call test_radb2()
    print *, "radb2 tests passed."

    call test_radb3()
    print *, "radb3 tests passed."

    call test_radb5()
    print *, "radb5 tests passed."

    call test_radb2_ido1()
    print *, "radb2 ido=1 branch tests passed."

    call test_radb2_ido3()
    print *, "radb2 ido=3 (OCL NOVREC loop) tests passed."

    call test_radb2_ido4()
    print *, "radb2 ido=4 (even IDO end branch) tests passed."

    call test_radb3_ido1()
    print *, "radb3 ido=1 branch tests passed."

    call test_radb3_ido3()
    print *, "radb3 ido=3 (OCL NOVREC loop) tests passed."

    call test_radb5_ido1()
    print *, "radb5 ido=1 branch tests passed."

    call test_radb5_ido3()
    print *, "radb5 ido=3 (after IDO=1 return) tests passed."

    call test_radf5_ido1()
    print *, "radf5 ido=1 branch tests passed."

    call test_radf5_ido3()
    print *, "radf5 ido=3 (after IDO=1 return) tests passed."

    call test_radfg_ido1()
    print *, "radfg ido=1 branch tests passed."

    call test_radbg_ido1()
    print *, "radbg ido=1 branch tests passed."

    call test_radbg_ido3_l1()
    print *, "radbg ido=3 l1=1 (NBD>=L1) tests passed."

    call test_radbg_ido3_l2()
    print *, "radbg ido=3 l1=2 (NBD>=L1 inner loop) tests passed."

    call test_radbg_ido3_l5()
    print *, "radbg ido=3 l1=5 (IDO<L1/NBD<L1 branches) tests passed."

    call test_radfg()
    print *, "radfg tests passed."

    call test_radbg()
    print *, "radbg tests passed."

    print *, "All fftpack tests completed."

  END SUBROUTINE run_fftpack_tests
