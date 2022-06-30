! This is a test for the NCEPLBS-ip library.
!
! Kyle Gerheiser June, 2021
program test_scalar_grib2
  use input_data_mod_grib2
  use interp_mod_grib2
  implicit none

  integer :: num_args, len, status
  character(len=32) :: grid_type, interp_opt

  num_args = command_argument_count()
  if (num_args /= 2) then
     print *, "Two command line arguments expected: grid, interpolation scheme"
     error stop
  end if

  call get_command_argument(1, grid_type, len, status)
  call get_command_argument(2, interp_opt, len, status)

  call read_input_data()
  call interp(grid_type, interp_opt)
end program test_scalar_grib2
