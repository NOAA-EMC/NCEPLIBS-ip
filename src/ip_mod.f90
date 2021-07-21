!> @file
!! @brief Top-level module for the ip library.
!! @author Kyle Gerheiser

!> Top-level module for the ip library which re-exports public routines such as ipolates, ipolatev, and gdswzd.
module ip_mod

  ! Make these constants public to everyone instead of
  ! using numbers directly
  use ip_interpolators_mod, only: BILINEAR_INTERP_ID, &
       BICUBIC_INTERP_ID, &
       NEIGHBOR_INTERP_ID, &
       BUDGET_INTERP_ID, &
       SPECTRAL_INTERP_ID, &
       NEIGHBOR_BUDGET_INTERP_ID
       
  use ipolates_mod
  use ipolatev_mod
  use gdswzd_mod
end module ip_mod
