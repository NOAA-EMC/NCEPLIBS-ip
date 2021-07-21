!> @file
!! @brief Top-level module to export interpolation routines and constants.
!! @author Kyle Gerheiser

!> Top-level module to export interpolation routines and constants.
!! @author Kyle Gerheiser
module ip_interpolators_mod
  implicit none

  use bilinear_interp_mod
  use bicubic_interp_mod
  use budget_interp_mod
  use neighbor_interp_mod
  use spectral_interp_mod
  use neighbor_budget_interp_mod

  !> @param Constant to choose BILINEAR interpolation method
  integer, parameter, public :: BILINEAR_INTERP_ID = 0
  !> @param Constant to choose BICUBIC interpolation method
  integer, parameter, public :: BICUBIC_INTERP_ID = 1
  !> @param Constant to choose NEIGBOR interpolation method
  integer, parameter, public :: NEIGHBOR_INTERP_ID = 2
  !> @param Constant to choose BUDGET interpolation method
  integer, parameter, public :: BUDGET_INTERP_ID = 3
  !> @param Constant to choose SPECTRAL interpolation method
  integer, parameter, public :: SPECTRAL_INTERP_ID = 4
  !> @param Constant to choose NEIGBOR_BUDGET interpolation method
  integer, parameter, public :: NEIGHBOR_BUDGET_INTERP_ID = 6

contains


end module ip_interpolators_mod

