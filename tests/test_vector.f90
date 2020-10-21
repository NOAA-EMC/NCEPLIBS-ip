program test_vector
  use input_data_mod
  use interp_mod
  implicit none

  integer :: num_args, i, len, status
  character(len=32) :: grid_type, interp_opt

  num_args = command_argument_count()
  if (num_args /= 2) then
     print *, "Two command line arguments expected: grid, interpolation scheme"
     error stop
  end if

  call get_command_argument(1, grid_type, len, status)
  call get_command_argument(2, interp_opt, len, status)
  
  call read_vector_input_data()
  call interp_vector(grid_type, interp_opt)
end program test_vector


