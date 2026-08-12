program test_ip_version
   use ip_mod, only: ip_version
   implicit none

   character(len=256) :: ctemp
   character(len=:), allocatable :: expected_version
   integer :: status

   if (command_argument_count() /= 1) then
      write(*, '(a)') "Usage: test_ip_version_fortran EXPECTED_VERSION"
      error stop 1
   endif

   call get_command_argument(1, ctemp, status=status)
   if (status /= 0) then
      write(*, '(a)') "Unable to determine command-line argument length"
      error stop 1
   endif

   expected_version=trim(ctemp)

   if (expected_version /= ip_version) then
      write(*, '(a)') "Version mismatch"
      write(*, '(2a)') "Expected: ", expected_version
      write(*, '(2a)') "Actual:   ", ip_version
      error stop 1
   endif

   write(*, '(2a)') "NCEPLIBS-ip version: ", ip_version

end program test_ip_version
