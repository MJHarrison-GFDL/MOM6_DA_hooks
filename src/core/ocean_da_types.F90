 module ocean_da_types_mod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains a set of data structures and interfaces for compiling the MOM6 DA
! driver code.
! Contact: Matthew.Harrison@noaa.gov and Feiyu.Lu@noaa.goy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef MAX_LEVS_FILE_
#define MAX_LEVS_FILE_ 1000
#endif

#ifndef MAX_LINKS_
#define MAX_LINKS_ 1
#endif

  use time_manager_mod, only : time_type

  implicit none

  private

  integer, parameter, public :: MAX_LEVELS_FILE = MAX_LEVS_FILE_ !< Controls record length for optimal storage
  integer, parameter, public :: MAX_LINKS = MAX_LINKS_ !< Maximum number of records per profile for storage for profiles
  integer, parameter, public :: UNKNOWN = 0

  integer, save, public :: TEMP_ID = 1
  integer, save, public :: SALT_ID = 2
  integer, save, public ::  SSH_ID = 3
  integer, save, public ::    U_ID = 4
  integer, save, public ::    V_ID = 5
  real, parameter, public :: MISSING_VALUE = -1.e10
  integer, save, public :: ODA_PFL = 1
  integer, save, public :: ODA_XBT = 2
  integer, save, public :: ODA_MRB = 3
  integer, save, public :: ODA_OISST = 4

  public :: copy_profile

!> Type for ocean state in DA space (same decomposition and vertical grid)
  type, public :: OCEAN_CONTROL_STRUCT
     integer :: ensemble_size
     real, pointer, dimension(:,:,:) :: SSH=>NULL() !<sea surface height (m) across ensembles
     real, pointer, dimension(:,:,:,:) :: h=>NULL() !<layer thicknesses (m or kg) across ensembles
     real, pointer, dimension(:,:,:,:) :: T=>NULL() !<layer potential temperature (degC) across ensembles
     real, pointer, dimension(:,:,:,:) :: S=>NULL() !<layer salinity (psu or g kg-1) across ensembles
     real, pointer, dimension(:,:,:,:) :: U=>NULL() !<layer zonal velocity (m s-1) across ensembles
     real, pointer, dimension(:,:,:,:) :: V=>NULL() !<layer meridional velocity (m s-1) across ensembles
     integer, dimension(:), pointer :: id_t=>NULL(), id_s=>NULL()  !< diagnostic IDs for temperature and salinity
     integer, dimension(:), pointer :: id_u=>NULL(), id_v=>NULL()     !< diagnostic IDs for zonal and meridional velocity
     integer, dimension(:), pointer :: id_ssh=>NULL()  !< diagnostic IDs for SSH
  end type OCEAN_CONTROL_STRUCT

  type, public :: ocean_profile_type
     integer :: inst_type !< A numeric code indicating the type of instrument (e.g. ARGO drifter, CTD, ...)
     logical :: initialized=.false. !< a True value indicates that this profile has been allocated for use
     logical :: colocated=.true. !< a True value indicated that the measurements of (num_variables) data are colocated in space-time
     integer :: ensemble_size !< size of the ensemble of model states used in association with this profile
     integer :: num_variables !< number of measurement types associated with this profile.
     integer, pointer, dimension(:) :: var_id !< variable ids are defined by the ocean_types module (e.g. TEMP_ID, SALT_ID)
     integer :: platform !< platform types are defined by platform class (e.g. MOORING, DROP, etc.) and instrument type (XBT, CDT, etc.)
     integer :: levels !< number of levels in the current profile
     integer :: basin_mask !<1:Southern Ocean, 2:Atlantic Ocean, 3:Pacific Ocean,
                           !! 4:Arctic Ocean, 5:Indian Ocean, 6:Mediterranean Sea, 7:Black Sea,
                           !!  8:Hudson Bay, 9:Baltic Sea, 10:Red Sea, 11:Persian Gulf
     integer :: profile_flag !< an overall flag for the profile
     real :: lat, lon !< latitude and longitude (degrees E and N)
     logical :: accepted !< logical flag to disable a profile
     integer :: nlinks !< number of links used to construct the profile (when reading from disk)
     type(time_type) :: time_window !< The time window associated with this profile [s]
     real, pointer, dimension(:) :: obs_error  !< The observation error by variable
     real  :: loc_dist   !< The impact radius of this observation (m)
     type(ocean_profile_type), pointer :: next=>NULL() !< all profiles are stored as linked list.
     type(ocean_profile_type), pointer :: prev=>NULL()
     type(ocean_profile_type), pointer :: cnext=>NULL() ! current profiles are stored as linked list.
     type(ocean_profile_type), pointer :: cprev=>NULL()
     integer :: nbr_xi, nbr_yi ! nearest neighbor model gridpoint for the profile
     real :: nbr_dist ! distance to nearest neighbor model gridpoint
     logical :: compute !< profile is within current compute domain
     real, dimension(:,:), pointer :: depth => NULL() !< depth of measurement [m]
     real, dimension(:,:), pointer :: data => NULL() !< data by variable type
     integer, dimension(:,:), pointer :: flag => NULL() !< flag by depth and variable type
     real, dimension(:,:,:), pointer :: forecast => NULL() !< ensemble member first guess
     real, dimension(:,:,:), pointer :: analysis => NULL() !< ensemble member analysis
     type(forward_operator_type), pointer :: obs_def => NULL() !< observation forward operator
     type(time_type) :: time !< profile FMS time type
     real :: i_index, j_index !< model longitude and latitude indices respectively
     real, dimension(:,:), pointer :: k_index !< model depth indices
     type(time_type) :: tdiff !< difference between model time and observation time
     character(len=128) :: filename
  end type ocean_profile_type

  type, public :: forward_operator_type
     integer :: num
     integer, dimension(2) :: state_size !< for
     integer, dimension(:), pointer :: state_var_index !< for flattened data
     integer, dimension(:), pointer :: i_index !< i-dimension index
     integer, dimension(:), pointer :: j_index !< j-dimension index
     real, dimension(:), pointer :: coef
  end type forward_operator_type

!> Grid information for ODA purposes, including arrays of
! lat, lon, depth, thickness, basin and land mask
  type, public :: grid_type
     real, pointer, dimension(:,:) :: x=>NULL(), y=>NULL()
     real, pointer, dimension(:,:,:) :: z=>NULL()
     real, pointer, dimension(:,:,:) :: h=>NULL()
     real, pointer, dimension(:,:) :: basin_mask => NULL()
     real, pointer, dimension(:,:,:) :: mask => NULL()
     real, pointer, dimension(:,:) :: bathyT => NULL()
     logical :: tripolar_N
     integer :: ni, nj, nk
  end type grid_type

contains

  subroutine copy_profile(Prof,Cprof)
    type(ocean_profile_type), pointer :: Prof
    type(ocean_profile_type), pointer :: Cprof

    if (.not. associated(CProf)) allocate(CProf)

    !print *,'prof%variables, prof%levels=',Prof%variable,Prof%levels
    !print *,'CProf%variable=',CProf%variable
    CProf%num_variables = Prof%num_variables
    allocate(CProf%var_id(Prof%num_variables))
    CProf%var_id = Prof%var_id
    CProf%platform = Prof%platform
    CProf%levels = Prof%levels
    CProf%basin_mask = Prof%basin_mask
    CProf%profile_flag = Prof%profile_flag
    CProf%lon = Prof%lon
    CProf%lat = Prof%lat
    CProf%accepted = Prof%accepted
    CProf%time_window = Prof%time_window
    allocate(CProf%obs_error(Prof%num_variables))
    CProf%obs_error = Prof%obs_error
    CProf%nbr_xi = Prof%nbr_xi
    CProf%nbr_yi = Prof%nbr_yi
    CProf%nbr_dist = Prof%nbr_dist
    CProf%compute = Prof%compute
    if (Prof%levels>0) then
       allocate(CProf%depth(Prof%num_variables,Prof%levels))
       CProf%depth=Prof%depth
       allocate(CProf%data(Prof%num_variables,Prof%levels))
       CProf%data=Prof%data
       allocate(CProf%k_index(Prof%num_variables,Prof%levels))
       CProf%k_index=Prof%k_index
    endif
    CProf%i_index = Prof%i_index
    CProf%j_index = Prof%j_index
    CProf%tdiff = Prof%tdiff
  end subroutine copy_profile

end module ocean_da_types_mod
