! This is a test from the NCEPLIBS-sp project.
!
! This test tests the ncpus() function.
!
! Kyle Gerheiser
program test_ncpus
  use sp_mod

  implicit none

  integer :: n

  n = ncpus()
#ifndef OPENMP
  if (n .ne. 1) stop 2
#endif
  
end program test_ncpus
