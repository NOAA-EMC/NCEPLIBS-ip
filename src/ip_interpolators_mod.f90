module ip_interpolators_mod
  ! re-export specific interpolation routines
  ! use bilinear_interpolator_scalar_mod
  ! use bicubic_interpolator_scalar_mod
  ! use neighbor_interpolator_scalar_mod
  ! use budget_interpolator_scalar_mod
  ! use spectral_interpolator_scalar_mod
  ! use neighbor_budget_interpolator_scalar_mod

  implicit none

  integer, parameter, public :: BILINEAR_INTERP_ID = 0
  integer, parameter, public :: BICUBIC_INTERP_ID = 1
  integer, parameter, public :: NEIGHBOR_INTERP_ID = 2
  integer, parameter, public :: BUDGET_INTERP_ID = 3
  integer, parameter, public :: SPECTRAL_INTERP_ID = 4
  integer, parameter, public :: NEIGHBOR_BUDGET_INTERP_ID = 6
end module ip_interpolators_mod

