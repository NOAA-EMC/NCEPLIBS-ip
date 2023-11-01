! This is a test from the NCEPLIBS-sp project.
!
! This test tests the ncpus() function.
!
! Kyle Gerheiser
program test_ncpus
  implicit none

  integer :: n, ncpus

  n = ncpus()
#ifndef OPENMP
  if (n .ne. 1) stop 2
#endif
  
end program test_ncpus
