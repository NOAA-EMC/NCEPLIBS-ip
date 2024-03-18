!> @file
!! @brief Top-level module to export interpolation routines and constants.
!! @author Kyle Gerheiser

!> Top-level module to export interpolation routines and constants.
!! @author Kyle Gerheiser
module ip_interpolators_mod
  use bilinear_interp_mod
  use bicubic_interp_mod
  use budget_interp_mod
  use neighbor_interp_mod
  use spectral_interp_mod
  use neighbor_budget_interp_mod
  implicit none

  !> @param Constant to choose BILINEAR interpolation method
  integer,parameter,public :: bilinear_interp_id=0
  !> @param Constant to choose BICUBIC interpolation method
  integer,parameter,public :: bicubic_interp_id=1
  !> @param Constant to choose NEIGBOR interpolation method
  integer,parameter,public :: neighbor_interp_id=2
  !> @param Constant to choose BUDGET interpolation method
  integer,parameter,public :: budget_interp_id=3
  !> @param Constant to choose SPECTRAL interpolation method
  integer,parameter,public :: spectral_interp_id=4
  !> @param Constant to choose NEIGBOR_BUDGET interpolation method
  integer,parameter,public :: neighbor_budget_interp_id=6

contains

endmodule ip_interpolators_mod

