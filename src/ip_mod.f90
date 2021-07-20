!> @file
!! @brief Top-level module for the ip library.
!! @author Kyle Gerheiser

!> Top-level module for the ip library which re-exports public routines such as ipolates, ipolatev, and gdswzd.
module ip_mod
  use ipolates_mod
  use ipolatev_mod
  use gdswzd_mod
end module ip_mod
