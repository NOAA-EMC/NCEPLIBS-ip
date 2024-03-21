!> @file
!! @brief Top-level module for the ip library.
!! @author Kyle Gerheiser

!> Top-level module for the ip library which re-exports public routines such as ipolates, ipolatev, and gdswzd.
module ip_mod

  ! Make these constants public to everyone instead of
  ! using numbers directly
  use ip_interpolators_mod,only:bilinear_interp_id, &
                                 bicubic_interp_id, &
                                 neighbor_interp_id, &
                                 budget_interp_id, &
                                 spectral_interp_id, &
                                 neighbor_budget_interp_id

  use ipolates_mod
  use ipolatev_mod
  use gdswzd_mod
endmodule ip_mod
