!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains a set of dummy interfaces for compiling the MOM6 DA
! driver code. These interfaces are not finalized and will be replaced by
! supported
! interfaces at some later date.
!
! 4/5/18
! matthew.harrison@noaa.gov
! feiyu.lu@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module eakf_oda_mod

  ! FMS shared modules
  use time_manager_mod, only : time_type, get_time
  use mpp_domains_mod, only : domain2d

  ! ODA Modules
  use oda_types_mod, only : ocean_profile_type, ocean_control_struct, grid_type
  use kdtree, only : kd_root

  implicit none

  public ensemble_filter

contains

  subroutine ensemble_filter(Prior, Posterior, Profiles, kdroot, Domain, oda_grid)
    type(ocean_control_struct), pointer, intent(in) :: Prior
    type(ocean_control_struct), pointer, intent(inout) :: Posterior
    type(ocean_profile_type), pointer, intent(in) :: Profiles
    type(kd_root), pointer, intent(inout) :: kdroot
    type(domain2d), intent(in) :: Domain
    type(grid_type), pointer, intent(in) :: oda_grid

    Posterior => Prior
    return

  end subroutine ensemble_filter

end module eakf_oda_mod
