! This is a test for the NCEPLIBS-ip library.
!
! Tests the SPVAR subroutine: compute variance by total wavenumber
! from a spectral coefficient array.
!
! Hang Lei, Aug 2026

program test_spvar
  use sp_mod
  implicit none

  integer, parameter :: iromb = 0, maxwv = 3
  integer, parameter :: ncoef = (maxwv+1)*((iromb+1)*maxwv+2)
  real :: q(ncoef)
  real :: qvar(0:(iromb+1)*maxwv)
  integer :: n
  real, parameter :: tol = 1.0e-6

  print *, 'Testing SPVAR...'

  ! Test 1: all-zero field -> all variances zero
  q = 0.0
  call spvar(iromb, maxwv, q, qvar)
  do n = 0, (iromb+1)*maxwv
    if (abs(qvar(n)) > tol) then
      print *, 'FAILED test 1: qvar(', n, ') =', qvar(n), ' expected 0'
      stop 1
    end if
  end do

  ! Test 2: unit field in the first wavenumber (n=0, l=0) coefficient
  ! Q(1) is the real part of (l=0, n=0); QVAR(0) = 0.5*Q(1)^2
  q = 0.0
  q(1) = 2.0
  call spvar(iromb, maxwv, q, qvar)
  if (abs(qvar(0) - 0.5*4.0) > tol) then
    print *, 'FAILED test 2: qvar(0) =', qvar(0), ' expected', 0.5*4.0
    stop 2
  end if
  ! Higher wavenumbers should still be zero
  do n = 1, (iromb+1)*maxwv
    if (abs(qvar(n)) > tol) then
      print *, 'FAILED test 2: qvar(', n, ') =', qvar(n), ' expected 0'
      stop 3
    end if
  end do

  ! Test 3: variances must be non-negative for any field
  q = 1.0
  call spvar(iromb, maxwv, q, qvar)
  do n = 0, (iromb+1)*maxwv
    if (qvar(n) < 0.0) then
      print *, 'FAILED test 3: negative variance qvar(', n, ') =', qvar(n)
      stop 4
    end if
  end do

  print *, 'SUCCESS!'

end program test_spvar
