program scalar_driver

  !-----------------------------------------------------------------
  ! Interpolate a global lat/lon grid of scalars to several
  ! grids of various projections using all ipolates 
  ! interpolation options.
  !-----------------------------------------------------------------

  use scalar_get_input_data

  implicit none

  call read_input_data

  call scalar_interp

  !print*,"- NORMAL TERMINATION"

  stop
end program scalar_driver
