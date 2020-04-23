!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! infrastructure for ocean data assimilation and model-observation misfit diagnostics.
! It relies heavily on the Flexible  Modeling System (FMS) which provides
! interfaces to communication and I/O layers.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_da_core_mod


  use mpp_domains_mod, only : domain2d
  use fms_mod, only : file_exist, read_data
  use fms_mod, only : open_namelist_file, check_nml_error, close_file
  use fms_mod, only : error_mesg, FATAL, NOTE
#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#endif
  use mpp_mod, only : mpp_sum, stdout, stdlog, mpp_sync_self
  use mpp_mod, only : mpp_pe, mpp_root_pe, mpp_error
  use mpp_io_mod, only : mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_NETCDF
  use mpp_io_mod, only : mpp_get_atts, mpp_get_info, mpp_get_fields, mpp_read, axistype, fieldtype, mpp_get_axes
  use mpp_io_mod, only : mpp_get_times
  use mpp_io_mod, only : mpp_get_axis_data, mpp_get_field_name
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : domain2d, mpp_get_global_domain
  use time_manager_mod, only : time_type, set_time, get_date
  use time_manager_mod, only : decrement_time, increment_time
  use time_manager_mod, only : operator( <= ), operator( >= ), operator(*), operator(/)
  use time_manager_mod, only : operator( + ), operator( - ), operator( > ), operator ( < )
  use time_manager_mod, only : get_time
  use get_cal_time_mod, only : get_cal_time
  use axis_utils_mod, only : frac_index
  use horiz_interp_type_mod, only: horiz_interp_type
  use horiz_interp_bilinear_mod, only : horiz_interp_bilinear_new
  use constants_mod, only : DEG_TO_RAD
  ! ODA_tools modules
  use ocean_da_types_mod, only : ocean_profile_type, grid_type
  use ocean_da_types_mod, only : forward_operator_type
  use ocean_da_types_mod, only : TEMP_ID, SALT_ID, MISSING_VALUE
  use ocean_da_types_mod, only : ODA_PFL, ODA_XBT, ODA_MRB, ODA_OISST
  use ocean_da_types_mod, only : UNKNOWN, MAX_LEVELS_FILE, MAX_LINKS
  use ocean_da_types_mod, only : copy_profile, ocean_control_struct
  use kdtree, only : kd_root, kd_search_nnearest, kd_init
  use loc_and_dist_mod, only : within_domain

  implicit none
  private

  public :: ocean_da_core_init
  public :: get_profiles, show_profiles, kd_root
  public :: apply_forward_operator

  ! Parameters
  integer, parameter :: num_obs_types = 5 !< number of observation types, change when adding new observation types
  integer, parameter :: PROFILE_FILE = 1
  integer, parameter :: SURFACE_FILE = 2
  integer, parameter :: ARGO_FILE = 3
  integer, parameter :: MOORING_FILE = 4
  integer, parameter :: VIRTUAL_MOORING_FILE = 5
  type(time_type) , dimension(num_obs_types) :: time_window !< time window for DROP, MOORING and SATELLITE data respectively

  ! ocean_obs_nml variables
  integer :: max_levels = 1000 !< maximium number of levels for a single profile (nd)
  integer :: max_neighbors=1 !< maximum number of neighbors for forward operator (nd)
  real, dimension(num_obs_types) :: obs_sbound = -89.0 !< observation domain limits by observation type (degrees)
  real, dimension(num_obs_types) :: obs_nbound = 89.0 !< observation domain limits by observation type (degrees)
  real :: depth_cut = 2000.0 !< depth below which observations are ignored (m)
  real, dimension(num_obs_types) :: data_window = 24.0 !< data window limit by observation type (hours)
  integer, dimension(num_obs_types) :: sec_offset = 0 !< offset to observation time (sec)
  integer, dimension(num_obs_types) :: day_offset = 0 !< offset to observation time (days)
  real, dimension(num_obs_types) :: temp_error = 0.001 !< error in temperature observations (degK)
  real, dimension(num_obs_types) :: salt_error = 0.001 !< error in salinity observations (psu)
  real, dimension(num_obs_types) :: temp_dist = 200.0e3 !< lateral impact radius for temperature observations (m)
  real, dimension(num_obs_types) :: salt_dist = 200.0e3 !< lateral impact radius for temperature observations (m)
  integer :: obs_days_plus = 0.0 !< use observations within time window ending at the initial model time plus (days)
  integer :: obs_days_minus = 0.0 !< use observations within time window starting at the initial model time minus (days)
  integer :: max_files = 100000 !< maximum number of files
  real :: shelf_depth = 0.0 !< seafloor depth below which to stop assimilating data in the column (m)
  real :: small_xy=1.e-10 !< A smallish number (m)

  namelist /ocean_obs_nml/ max_levels, max_neighbors,obs_sbound, obs_nbound, depth_cut, &
          data_window, sec_offset, day_offset, temp_error, salt_error, &
          temp_dist, salt_dist, shelf_depth, max_files, obs_days_minus, obs_days_plus

contains

!> This subroutine initializes the core ocean da module
  subroutine ocean_da_core_init(Domain, global_grid, local_grid, Profiles, model_time, kdroot)
    type(domain2d), pointer, intent(in) :: Domain
    type(grid_type), pointer, intent(in) :: global_grid
    type(grid_type), pointer, intent(in) :: local_grid
    type(ocean_profile_type), pointer :: Profiles
    type(time_type), intent(in) :: model_time
    type(kd_root), pointer :: kdroot

    type obs_entry_type
       character(len=128) :: filename
       character(len=32)  :: file_type
       integer               :: time_interval
    end type obs_entry_type


    type(time_type)  :: time_s, time_e
    integer :: data_seconds
    integer :: i, j, n, ni, nj
    integer :: isc,iec,jsc,jec
    integer :: isd,ied,jsd,jed
    integer :: halox, haloy, blk
    integer :: ioun, io_status, ierr
    integer :: stdout_unit, stdlog_unit
    integer :: nfiles, nrecs, unit
    integer :: i_offset, j_offset ! translates local to global indices
    integer, dimension(:), allocatable :: filetype
    integer, allocatable, dimension(:) :: lon1d,lat1d
    real, allocatable, dimension(:) :: glon1d, glat1d
    character(len=256) :: record
    character(len=128), dimension(:), allocatable :: input_files
    integer, dimension(:), allocatable :: time_interval
    type(obs_entry_type) :: tbl_entry

    tbl_entry%filename=''
    tbl_entry%file_type=''
    tbl_entry%time_interval=0

    stdout_unit = stdout()
    stdlog_unit = stdlog()
    ioun = open_namelist_file()
    read(UNIT=ioun, NML=ocean_obs_nml, IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_obs_nml')
    call close_file(ioun)
    write (UNIT=stdlog_unit, NML=ocean_obs_nml)
    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, NML=ocean_obs_nml)
    end if
    ! Allocate filetype* and input_files* variables
    allocate(filetype(max_files), input_files(max_files),time_interval(max_files))
    filetype = -1
    input_files = ''
    time_s = decrement_time(model_time, 0, obs_days_minus)
    time_e = increment_time(model_time, 0, obs_days_plus)
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    i_offset = isc - 1 !< offset to reference isc to 1
    j_offset = jsc - 1 !< offset to reference jsc to 1
    isc = isc - isd + 1 ; iec = iec - isd + 1 ; jsc = jsc - jsd + 1 ; jec = jec - jsd + 1
    ied = ied - isd + 1 ; jed = jed - jsd + 1 ; isd = 1 ; jsd = 1
    halox = isc - isd
    haloy = jsc - jsd
    blk = (jed-jsd+1)*(ied-isd+1)
    ni = global_grid%ni; nj = global_grid%nj
    allocate(lon1d(blk))
    allocate(lat1d(blk))
    allocate(glon1d(blk))
    allocate(glat1d(blk))
    if (associated(kdroot)) then
       call error_mesg('ocean_da_core_mod: ocean_da_core_init', 'kdroot already allocated',FATAL)
    endif
    allocate(kdroot)
    do i = isd, ied; do j = jsd, jed
      lon1d((j-1)*ied+i) = i ! i-index  of data in blk
      lat1d((j-1)*ied+i) = j ! j-index  of data in blk
      glon1d((j-1)*ied+i) = local_grid%x(i,j) ! longitude coordinate of data in blk
      glat1d((j-1)*ied+i) = local_grid%y(i,j) ! latitude coordinate of data in blk
    enddo; enddo
    call kd_init(kdroot, glon1d, glat1d)
    ! time window for DROP, MOORING and SATELLITE data respectively
    ! will be available from namelist
    do i=1,num_obs_types
      data_seconds = data_window(i) * 3600
      time_window(i) = set_time(data_seconds,0)
    enddo
    nfiles = 0
    nrecs=0
    call mpp_open(unit, 'ocean_obs_table', action=MPP_RDONLY)
    read_obs: do while ( nfiles <= max_files )
       read (UNIT=unit, FMT='(A)', IOSTAT=io_status) record
       if ( io_status < 0 ) then
          exit read_obs
       else if ( io_status > 0 ) then
          cycle read_obs
       else
          nrecs = nrecs + 1
          if ( record(1:1) == '#' ) cycle read_obs
          read ( UNIT=record, FMT=*, IOSTAT=io_status ) tbl_entry
          if ( io_status < 0 ) then
             print *,tbl_entry%filename
             print *,tbl_entry%file_type
             print *,tbl_entry%time_interval
             call error_mesg('ocean_da_core_mod::init_observations', 'error in obs_table entry format', FATAL)
          else if ( io_status > 0 ) then
             cycle read_obs
          else
             nfiles = nfiles + 1
             write (UNIT=stdout_unit, FMT='("Obs filename:",A," type:",A)') tbl_entry%filename,tbl_entry%file_type
             input_files(nfiles) = tbl_entry%filename
             select case ( trim(tbl_entry%file_type) )
             case ('profiles')
                filetype(nfiles) = PROFILE_FILE
             case ('argo')
                filetype(nfiles) = ARGO_FILE
             case ('mooring')
                filetype(nfiles) = MOORING_FILE
             case('virtual_mooring')
                filetype(nfiles) = VIRTUAL_MOORING_FILE
                if (tbl_entry%time_interval<=0) then
                   print *,'virtual mooring time interval must be > 0'
                endif
                time_interval(nfiles)=tbl_entry%time_interval ! input time interval in model units [s]
             case default
                print *,tbl_entry
                call error_mesg('ocean_da_core_mod::init_observations', 'error in obs_table entry format', FATAL)
             end select
          end if
       end if
    end do read_obs
    if ( nfiles > max_files ) then
       call error_mesg('ocean_da_core_mod::init_observations', 'number of obs files exceeds max_files parameter', FATAL)
    end if
    CALL mpp_close(unit)
    if( .not. associated(Profiles) ) allocate(Profiles)
    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("TEMP_ID = ",I5)') TEMP_ID
    end if
    do n=1, nfiles
      select case ( filetype(n) )
!      case (ARGO_FILE)
!         call open_argo_dataset(Profiles, Domain, global_grid, &
!              trim(input_files(n)), time_s, time_e)
!      case (MOORING_FILE)
!         call open_mooring_dataset(Profiles, Domain, global_grid, &
!              trim(input_files(n)), time_s, time_e)
      case (VIRTUAL_MOORING_FILE)
         call open_virtual_mooring_dataset(Profiles, Domain, global_grid, &
              local_grid, kdroot, trim(input_files(n)), time_s, time_e, time_interval(n))
      case default
         call error_mesg('ocean_da_core_mod::init_observations', 'filetype not currently supported for temp_obs', FATAL)
      end select
    end do

    ! Deallocate local arrays before exiting routine
    deallocate(glat1d, glon1d)
    deallocate(filetype, input_files)

  end subroutine ocean_da_core_init


  ! subroutine open_argo_dataset(Profiles, Domain, T_grid, &
  !                 filename, time_start, time_end, localize)
  !   type(ocean_profile_type), pointer :: Profiles !< This is a recursive list of profiles
  !   type(domain2d), pointer, intent(in) :: Domain !< MOM grid type for the local domain
  !   type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
  !   character(len=*), intent(in) :: filename !< filename containing profile data
  !   type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
  !   logical, intent(in), optional :: localize !< Localize the observations to the current computational domain

  !   real :: lon, lat, time, rlink, var_type
  !   integer :: ni, nj, nk
  !   real :: ri0, rj0
  !   real, dimension(MAX_LEVELS_FILE) :: depth, data, data_s
  !   type(ocean_profile_type), pointer :: Prof

  !   integer :: unit, ndim, nvar, natt, nstation, max_profiles
  !   integer :: stdout_unit
  !   integer :: inst_type, var_id
  !   integer :: station_count, station_link
  !   integer :: num_levs, m, k, kk, i, j, i0, j0, k0, nlevs, a, nn, nlinks
  !   integer :: yr, mon, day, hr, min, sec
  !   integer :: ii, jj

  !   logical :: data_is_local, localize_data, cont
  !   logical :: data_in_period, data_in_compute
  !   real, dimension(MAX_LEVELS_FILE) :: flag

  !   character(len=32) :: fldname, axisname, time_units
  !   character(len=138) :: emsg_local

  !   type(time_type) :: obs_time, profile_time
  !   type(axistype), pointer :: depth_axis, station_axis
  !   type(axistype), allocatable, dimension(:), target :: axes
  !   type(fieldtype), allocatable, dimension(:), target :: fields
  !   type(fieldtype), pointer :: field_lon, field_lat, field_time, field_depth
  !   type(fieldtype), pointer :: field_t, field_s, field_link, field_var_type

  !   integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  !   integer :: isg, ieg, jsg, jeg, halox, haloy, lon_len, blk
  !   integer :: profile_count = 0, temp_count = 0, num_variables = 0
  !   type(horiz_interp_type) :: Interp
  !   real :: lon_out(1, 1), lat_out(1, 1)
  !   real :: lat_bound = 59.0
  !   integer :: inds(1), r_num
  !   real :: dist(1), frac_lon, frac_lat, frac_k
  !   real, dimension(6) :: coef
  !   integer, dimension(8) :: state_index

  !   Prof=>Profiles
  !   do while (associated(Prof%next))
  !     Prof=>Prof%next
  !   end do
  !   if ( PRESENT(localize) ) then
  !      localize_data = localize
  !   else
  !      localize_data = .true.
  !   end if
  !   ni = T_grid%ni; nj = T_grid%nj; nk = T_grid%nk
  !   call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
  !   call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
  !   call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
  !   lon_len = ied-isd+1
  !   blk = (jed-jsd+1)*lon_len
  !   stdout_unit = stdout()
  !   inst_type = ODA_PFL

  !   call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
  !   call mpp_get_info(unit, ndim, nvar, natt, nstation)
  !   write (UNIT=stdout_unit, FMT='("Opened argo dataset: ",A)') trim(filename)
  !   if (nstation .EQ. 0) then
  !     write(UNIT=stdout_unit, FMT='("There are ZERO records in this dataset.")')
  !     call mpp_close(unit)
  !     return
  !   end if
  !   ! get axis information
  !   allocate(axes(ndim))
  !   call mpp_get_axes(unit, axes)
  !   do i=1, ndim
  !      call mpp_get_atts(axes(i), name=axisname)
  !      select case ( trim(axisname) )
  !      case ('depth_index')
  !         depth_axis => axes(i)
  !      case ('station_index')
  !         station_axis => axes(i)
  !      end select
  !   end do

  !   ! get field information
  !   allocate(fields(nvar))
  !   call mpp_get_fields(unit, fields)
  !   field_lon=>NULL();field_lat=>NULL();field_time=>NULL()
  !   field_t=>NULL();field_s=>NULL();field_var_type=>NULL()
  !   field_depth=>NULL();field_link=>NULL()
  !   do i=1, nvar
  !     call mpp_get_atts(fields(i), name=fldname)
  !     select case (trim(lowercase(fldname)))
  !     case ('longitude')
  !       field_lon => fields(i)
  !     case ('latitude')
  !       field_lat => fields(i)
  !     case ('juld')
  !       field_time => fields(i)
  !     case ('temp')
  !       field_t => fields(i)
  !     case ('psal')
  !       field_s => fields(i)
  !     case ('pres')
  !       field_depth => fields(i)
  !     end select
  !   end do

  !   call mpp_get_atts(depth_axis, len=nlevs)

  !   if ( nlevs > MAX_LEVELS_FILE ) then
  !      call error_mesg('ocean_da_core_mod::open_profile_dataset', 'increase parameter MAX_LEVELS_FILE', FATAL)
  !   else if (nlevs < 1) then
  !      call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Value of nlevs is less than 1.', FATAL)
  !   end if

  !   if ( .NOT.ASSOCIATED(field_t) .and. .not. ASSOCIATED(field_s) ) then
  !      call error_mesg('ocean_da_core_mod::open_profile_dataset',&
  !           & 'profile dataset not used because data not needed for Analysis', NOTE)
  !      return
  !   end if

  !   write(UNIT=stdout_unit, FMT='("There are ",I8," records in this dataset.")') nstation
  !   write(UNIT=stdout_unit, FMT='("Searching for profiles . . .")')

  !   call mpp_get_atts(field_time, units=time_units)

  !   station_count = 1
  !   profile_count = 0
  !   cont = .true.

  !   do while ( cont )
  !      depth = MISSING_VALUE  ! snz add
  !      data = MISSING_VALUE   ! snz add
  !      num_variables=0
  !      call mpp_read(unit, field_lon, lon, tindex=station_count)
  !      call mpp_read(unit, field_lat, lat, tindex=station_count)
  !      call mpp_read(unit, field_time, time, tindex=station_count)
  !      if (.not. associated(field_t) .and. .not. associated(field_s)) then
  !         station_count = station_count + 1
  !         if ( station_count .gt. nstation ) cont = .false.
  !         cycle
  !      end if
  !      data_is_local = .false.
  !      data_in_period = .false.
  !      if ( lon .lt. 0.0 ) lon = lon + 360.0
  !      if ( lon .gt. 360.0 ) lon = lon - 360.0
  !      if ( lon .gt. 60.0 ) lon = lon - 360.0
  !      if ( lat < obs_sbound(inst_type) .or. lat > obs_nbound(inst_type) ) then
  !        station_count = station_count + 1
  !        if ( station_count .gt. nstation ) cont = .false.
  !        cycle
  !      end if
  !      obs_time = get_cal_time(time, time_units, 'julian')
  !      profile_time = increment_time(obs_time, sec_offset(inst_type),day_offset(inst_type))
  !      if ( profile_time >= time_start .and. profile_time <= time_end ) data_in_period = .true.
  !      if ( .not. data_in_period ) then
  !        station_count = station_count + 1
  !        if ( station_count .gt. nstation ) cont = .false.
  !        cycle
  !      end if
  !      if ( localize_data ) then
  !        call kd_search_nnearest(kdroot, lon, lat, &
  !                1, inds, dist, r_num, .false.)
  !        data_is_local = within_domain(lon1d(inds(1)), lat1d(inds(1)), isd+1, ied-1, jsd+1, jed-1, ni, nj)
  !        data_in_compute = within_domain(lon1d(inds(1)), lat1d(inds(1)), isc, iec, jsc, jec, ni, nj)
  !      else
  !        data_is_local = .true.
  !      end if
  !      if (.not. data_is_local) then
  !         station_count = station_count + 1
  !         if ( station_count .gt. nstation ) cont = .false.
  !         cycle
  !      end if
  !      profile_count = profile_count + 1
  !      call mpp_read(unit, field_depth, depth(1:nlevs), tindex=station_count)
  !      if ( associated(field_t) ) then
  !        call mpp_read(unit, field_t, data(1,1:nlevs), tindex=station_count)
  !      else if ( associated(field_s) ) then
  !        call mpp_read(unit, field_s, data(2,1:nlevs), tindex=station_count)
  !      end if
  !      num_levs = 0
  !      do k=1, MAX_LEVELS_FILE
  !         flag(k) = 1.0
  !         if ( depth(k) > depth_cut ) depth(k) = MISSING_VALUE
  !         if ( data(k) .eq. 99999. .or. depth(k) .eq. MISSING_VALUE) then !< ARGO _FillValue
  !            flag(k) = 0.0
  !         else
  !            num_levs = num_levs + 1
  !         end if
  !      end do

  !      if ( num_levs == 0 ) then
  !         station = station + 1
  !         if ( station_count .gt. nstation ) cont = .false.
  !         cycle
  !      end if

  !      ! allocate profile structure content and put in data
  !      Prof%num_variables=0
  !      if (associated(field_t)) Prof%num_variables=Prof%num_variables+1
  !      if (associated(field_s)) Prof%num_variables=Prof%num_variables+1
  !      allocate(Prof%depth(Prof%num_variables,num_levs));Prof%depth(:,:)=MISSING_VALUE
  !      allocate(Prof%var_id(Profnum_variables))
  !      allocate(Prof%data(Prof%num_variables,num_levs));Prof%data(:,:)=MISSING_VALUE
  !      allocate(Prof%flag(Prof%num_variables,num_levs));Prof%flag(:,:)= 0
  !      allocate(Prof%obs_error(Prof%num_variables));Prof%obs_error(:)= 0.
  !      allocate(Prof%loc_dist(Prof%num_variables));Prof%loc_dist(:)= 0.
  !      Prof%inst_type = inst_type
  !      Prof%levels = num_levs
  !      Prof%lat = lat; Prof%lon = lon
  !      Prof%nbr_xi = lon1d(inds(1)); Prof%nbr_yi = lat1d(inds(1))
  !      Prof%nbr_dist = dist(1)
  !      Prof%compute = data_in_compute
  !      Prof%time_window = time_window(inst_type)
  !      k=1
  !      if ( associated(field_t) ) then
  !         Prof%obs_error(k) = temp_error(inst_type)
  !         Prof%loc_dist(k) = temp_dist(inst_type)
  !         Prof%var_id(k) = TEMP_ID
  !         k=k+1
  !      endif
  !      if ( associated(field_s )) then
  !         Prof%obs_error(k) = salt_error(inst_type)
  !         Prof%loc_dist(k) = salt_dist(inst_type)
  !         Prof%var_id(k) = SAlT_ID
  !         k=k+1
  !      end if
  !      do k=1, MAX_LEVELS_FILE
  !        Prof%depth(k) = depth(k)
  !        do m=1,Prof%num_variables
  !          Prof%data(m,k) = data(m,k)
  !          Prof%flag(m,k) = flag(m,k)
  !        end do
  !      enddo
  !      Prof%time = profile_time
  !      i0 = lon1d(inds(1)); j0 = lat1d(inds(1))
  !      Prof%i_index = i0
  !      Prof%j_index = j0
  !      Prof%accepted = .true.
  !      if (i0 < 1 .or. j0 < 1) then
  !         Prof%accepted = .false.
  !      else
  !         Prof%basin_mask = T_grid%basin_mask(i0,j0)
  !      end if
  !      if ( Prof%accepted ) then ! check surface land-sea mask (nearest point) and depth of ocean
  !         if (T_grid%mask(i0,j0,1) == 0.0 ) then
  !            Prof%accepted = .false.
  !         end if
  !         if (T_grid%bathyT(i0,j0) < shelf_depth ) then
  !            Prof%accepted = .false.
  !         end if
  !      end if ! check surface land-sea mask (nearest point) and depth of ocean

  !      if ( Prof%accepted ) then ! determine vertical position and check mask at depth
  !         allocate(Prof%k_index(Prof%levels))
  !         do k=1, Prof%levels
  !            Prof%k_index(k) = frac_index(Prof%depth(k), (/T_grid%z(i0,j0,:)/))
  !            if ( Prof%k_index(k) < 1.0 ) then
  !               if ( Prof%depth(k) < T_grid%z(i0,j0,1) ) then
  !                  Prof%k_index(k) = 0.0
  !               else if ( Prof%depth(k) > T_grid%z(i0,j0,nk) ) then
  !                   Prof%k_index(k) = real(nk)
  !                   Prof%flag(k)= 0.0
  !               end if
  !            end if
  !            if ( Prof%k_index(k) > real(nk) ) then
  !               call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is greater than nk', FATAL)
  !            else if ( Prof%k_index(k) < 0.0 ) then
  !               call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is less than 0', FATAL)
  !            end if
  !            k0 = floor(Prof%k_index(k))

  !            if ( k0 >= 1 ) THEN
  !                  do m=1,Prof%num_variables
  !                    if ( Prof%flag(m,k) > 0.0 ) then ! flag
  !                       if ( T_grid%mask(i0,j0,k0) == 0.0 .or. T_grid%mask(i0,j0,min(k0+1,T_grid%nk)) == 0.0) then
  !                          Prof%flag(m,k) = 0.0
  !                       end if

  !                       if ( Prof%data(m,k) == MISSING_VALUE .or. Prof%depth(m,k) == MISSING_VALUE ) then
  !                          Prof%flag(m,k) = 0.0
  !                       end if
  !                    end if ! flag
  !                  enddo
  !            end if
  !          end do
  !      end if ! determine vertical position and check mask at depth

  !      if ( Prof%accepted ) then ! calculate forward operator indices and weights
  !         state_size=(/lon_len,lat_len/)
  !         call assign_forward_operator(max_neighbors,state_size, dist, Prof)
  !      endif ! calculate forward operator indices and weights

  !      if ( station_count .gt. nstation ) cont = .false.
  !      allocate(Prof%next) ! allocate next profile and link it to current one
  !      Prof%next%prev=>Prof
  !      Prof=>Prof%next
  !   end do

  !   call mpp_sync_self()
  !   call mpp_close(unit)
  ! end subroutine open_argo_dataset

  ! subroutine open_mooring_dataset(Profiles, Domain, T_grid, &
  !                 filename, time_start, time_end, obs_variable, localize)
  !   type(ocean_profile_type), pointer :: Profiles
  !   !< This is an unstructured recursive list of profiles
  !   !< which are either within the localized domain corresponding
  !   !< to the Domain argument, or the global profile list
  !   type(domain2d), pointer, intent(in) :: Domain !< MOM grid type for the local domain
  !   type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
  !   character(len=*), intent(in) :: filename !< filename containing profile data
  !   type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
  !   integer, intent(in), optional :: obs_variable !< If present, then extract corresponding data
  !   !< from file, otherwise, extract all available data which.
  !   logical, intent(in), optional :: localize !< Localize the observations to the current computational domain

  !   real :: lon, lat, time
  !   integer :: nlat, nlon, ndepth, ntime
  !   integer :: ni, nj, nk
  !   real :: ri0, rj0
  !   real, allocatable, dimension(:) :: data, quality
  !   real, allocatable, dimension(:) :: flag
  !   type(ocean_profile_type), pointer :: Prof

  !   integer :: unit, ndim, nvar, natt, max_profiles
  !   integer :: stdout_unit
  !   integer :: inst_type, var_id
  !   integer :: num_levs, k, t1, kk, i, j, i0, j0, k0, nlevs, a, nn, nlinks
  !   integer :: yr, mon, day, hr, min, sec
  !   integer :: ii, jj

  !   logical :: data_is_local, localize_data, data_in_period

  !   character(len=32) :: fldname, axisname, time_units
  !   character(len=138) :: emsg_local

  !   type(time_type) :: mooring_time, obs_time
  !   type(axistype), pointer :: lon_axis, lat_axis, time_axis, depth_axis
  !   type(axistype), allocatable, dimension(:), target :: axes
  !   type(fieldtype), allocatable, dimension(:), target :: fields
  !   type(fieldtype), pointer :: field_t, field_quality

  !   real, allocatable, dimension(:) :: lons, lats, depths, times
  !   real, allocatable, dimension(:,:,:,:) :: mooring_obs, mooring_quality
  !   integer :: mooring_size(4)

  !   integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
  !   integer :: isg, ieg, jsg, jeg, halox, haloy, lon_len, blk
  !   integer :: mooring_count = 0
  !   type(horiz_interp_type) :: Interp
  !   real :: lon_out(1, 1), lat_out(1, 1)
  !   real :: lat_bound = 59.0
  !   integer :: inds(1), r_num
  !   real :: dist(1), frac_lon, frac_lat, frac_k
  !   real, dimension(6) :: coef
  !   integer, dimension(8) :: state_index

  !   Prof=>Profiles
  !   do while (associated(Prof%next))
  !     Prof=>Prof%next
  !   end do

  !   if ( PRESENT(localize) ) then
  !      localize_data = localize
  !   else
  !      localize_data = .true.
  !   end if

  !   ni = T_grid%ni; nj = T_grid%nj; nk = T_grid%nk
  !   call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
  !   call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
  !   call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
  !   lon_len = ied-isd+1
  !   blk = (jed-jsd+1)*lon_len
  !   stdout_unit = stdout()

  !   inst_type = ODA_MRB
  !   var_id = obs_variable

  !   call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
  !   call mpp_get_info(unit, ndim, nvar, natt, ntime)

  !   write (UNIT=stdout_unit, FMT='("Opened mooring dataset: ",A)') trim(filename)

  !   ! get axis information
  !   allocate(axes(ndim))
  !   call mpp_get_axes(unit, axes)
  !   do i=1, ndim
  !      call mpp_get_atts(axes(i), name=axisname)
  !      select case ( trim(axisname) )
  !      case ('lon')
  !         lon_axis => axes(i)
  !      case ('lat')
  !         lat_axis => axes(i)
  !      case ('depth')
  !         depth_axis => axes(i)
  !      case ('time')
  !         time_axis => axes(i)
  !      end select
  !   end do

  !   call mpp_get_atts(lon_axis,len=nlon)
  !   call mpp_get_atts(lat_axis,len=nlat)
  !   call mpp_get_atts(depth_axis,len=ndepth)
  !   call mpp_get_atts(time_axis,len=ntime)
  !   call mpp_get_atts(time_axis, units=time_units)

  !   allocate(lons(nlon), lats(nlat), depths(ndepth), times(ntime))
  !   allocate(mooring_obs(ntime,ndepth,nlat,nlon))
  !   allocate(mooring_quality(ntime,ndepth,nlat,nlon))
  !   allocate(data(ndepth),quality(ndepth))
  !   allocate(flag(ndepth))

  !   call mpp_get_axis_data(lon_axis, lons)
  !   call mpp_get_axis_data(lat_axis, lats)
  !   call mpp_get_axis_data(depth_axis, depths)
  !   call mpp_get_axis_data(time_axis, times)

  !   ! get field information
  !   allocate(fields(nvar))
  !   call mpp_get_fields(unit, fields)
  !   field_t=>NULL()
  !   do i=1, nvar
  !     call mpp_get_atts(fields(i), name=fldname)
  !     select case (trim(fldname))
  !     case ('T_20')
  !       field_t => fields(i)
  !     case ('QT_5020')
  !       field_quality => fields(i)
  !     end select
  !   end do

  !   call mpp_get_atts(field_t, siz=mooring_size)

  !   call mpp_read(unit, field_t, mooring_obs)
  !   call mpp_read(unit, field_quality, mooring_quality)
  !   do t1=1, ntime
  !     data_in_period = .false.
  !     time = times(t1)
  !     obs_time = get_cal_time(time, time_units, 'julian')
  !     mooring_time = increment_time(obs_time, sec_offset(inst_type),day_offset(inst_type))

  !     if ( mooring_time >= time_start .and. mooring_time <= time_end ) data_in_period = .true.
  !     if ( .not. data_in_period ) cycle

  !     do j=1, nlat
  !       do i=1, nlon
  !         lon = lons(i)
  !         lat = lats(j)
  !         data = MISSING_VALUE   ! snz add
  !         data_is_local = .false.

  !         if ( lon .lt. 0.0 ) lon = lon + 360.0
  !         if ( lon .gt. 360.0 ) lon = lon - 360.0
  !         if ( lon .gt. 60.0 ) lon = lon - 360.0

  !         if ( lat < obs_sbound(inst_type) .or. lat > obs_nbound(inst_type) ) cycle

  !         if ( localize_data ) then
  !           call kd_search_nnearest(kdroot, lon, lat, &
  !                   1, inds, dist, r_num, .false.)
  !           data_is_local = within_domain(lon1d(inds(1)), lat1d(inds(1)), &
  !                   isd+1, ied-1, jsd+1, jed-1, ni, nj)
  !         else
  !           data_is_local = .true.
  !         end if

  !         if (.not. data_is_local) cycle

  !         mooring_count = mooring_count + 1
  !         !type_count(inst_type) = type_count(inst_type)+1

  !         num_levs = 0
  !         data=mooring_obs(t1,:,j,i)
  !         quality=mooring_quality(t1,:,j,i)
  !         do k=1, ndepth
  !            flag(k) = 1.0
  !            if ( data(k)>50 .or. data(k)<-10 .or. quality(k)<0.5 .or. quality(k)>2.5) then
  !               flag(k) = 0.0
  !            else
  !               num_levs = num_levs + 1
  !            end if
  !         end do
  !         if (num_levs == 0) cycle

  !         !! allocate profile structure content and put in data
  !         allocate(Prof%depth(num_levs))
  !         allocate(Prof%data(num_levs))
  !         allocate(Prof%flag(num_levs))
  !         Prof%variable = var_id
  !         Prof%inst_type = inst_type
  !         Prof%levels = num_levs
  !         Prof%lat = lat; Prof%lon = lon
  !         Prof%nbr_xi = lon1d(inds(1)); Prof%nbr_yi = lat1d(inds(1))
  !         Prof%nbr_dist = dist(1)
  !         Prof%time_window = time_window(inst_type)
  !         Prof%impact_levels = impact_levels(inst_type)
  !         Prof%temp_to_salt = temp_to_salt(inst_type)
  !         Prof%salt_to_temp = salt_to_temp(inst_type)
  !         Prof%obs_error = temp_error(inst_type)
  !         Prof%loc_dist = temp_dist(inst_type)
  !         Prof%time = mooring_time

  !         kk = 1
  !         do k=1, ndepth
  !            if ( flag(k) > 0. ) then
  !              if ( kk > Prof%levels ) then
  !               call error_mesg('ocean_da_core_mod::open_profile_dataset',&
  !                    & 'Loop value "kk" is greater than profile levels', FATAL)
  !              end if
  !              Prof%depth(kk) = depths(k)
  !              Prof%data(kk) = data(k)
  !              Prof%flag(kk) = flag(k)
  !              kk = kk + 1
  !            end if
  !         end do

  !         if ( lat < lat_bound ) then ! calculate interpolation coefficients
  !            ri0 = frac_index(lon, T_grid%x(:,jsg))
  !            rj0 = frac_index(lat, T_grid%y(isg,:))
  !            i0 = floor(ri0)
  !            j0 = floor(rj0)
  !            if ( i0 > ieg .or. j0 > jeg ) then
  !               write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') mpp_pe(), i0, j0
  !               call error_mesg('ocean_da_core_mod::open_profile_dataset',&
  !                    & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
  !            end if
  !            Prof%i_index = ri0
  !            Prof%j_index = rj0
  !         else ! tripolar grids
  !            lon_out(1,1) = lon*DEG_TO_RAD
  !            lat_out(1,1) = lat*DEG_TO_RAD
  !            call horiz_interp_bilinear_new (Interp, T_grid%x*DEG_TO_RAD, T_grid%y*DEG_TO_RAD,&
  !                 & lon_out, lat_out, new_search=.true., no_crash_when_not_found=.true.)

  !            if ( Interp%wti(1,1,2) < 1.0 ) then
  !               i0 = Interp%i_lon(1,1,1)
  !            else
  !               i0 = Interp%i_lon(1,1,2)
  !            end if
  !            if ( Interp%wtj(1,1,2) < 1.0 ) then
  !               j0 = Interp%j_lat(1,1,1)
  !            else
  !               j0 = Interp%j_lat(1,1,2)
  !            end if
  !            if ( i0 > ieg .or. j0 > jeg ) then
  !               write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
  !               call error_mesg('ocean_da_core_mod::open_profile_dataset',&
  !                    & 'For tripolar grids, either i0 > ieg or j0 > jeg', FATAL)
  !            end if
  !            if ( Interp%wti(1,1,2) < 1.0 ) then
  !               Prof%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
  !            else
  !               Prof%i_index =Interp%i_lon(1,1,2)
  !            end if
  !            if (Interp%wtj(1,1,2) < 1.0) then
  !               Prof%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
  !            else
  !               Prof%j_index =Interp%j_lat(1,1,2)
  !            end if
  !         end if ! interpolation coefficients

  !         Prof%accepted = .true.

  !         if (i0 < 1 .or. j0 < 1) then
  !            Prof%accepted = .false.
  !         else
  !            Prof%basin_mask = T_grid%basin_mask(lon1d(inds(1)),lat1d(inds(1)))
  !         end if

  !         if ( Prof%accepted ) then ! check surface land-sea mask and depth of ocean around profile location
  !            if ( i0 /= ieg .and. j0 /= jeg ) then
  !               if (T_grid%mask(i0,j0,1) == 0.0 .or.&
  !                    & T_grid%mask(i0+1,j0,1) == 0.0 .or.&
  !                    & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
  !                    & T_grid%mask(i0+1,j0+1,1) == 0.0 ) then
  !                  Prof%accepted = .false.
  !               end if
  !               if (T_grid%bathyT(i0,j0) < shelf_depth .or.&
  !                    & T_grid%bathyT(i0+1,j0) < shelf_depth .or.&
  !                    & T_grid%bathyT(i0,j0+1) < shelf_depth .or.&
  !                    & T_grid%bathyT(i0+1,j0+1) < shelf_depth ) then
  !                  Prof%accepted = .false.
  !               end if
  !            else if ( i0 == ieg .and. j0 /= jeg ) then
  !               if (T_grid%mask(i0,j0,1) == 0.0 .or.&
  !                    & T_grid%mask(1,j0,1) == 0.0 .or.&
  !                    & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
  !                    & T_grid%mask(1,j0+1,1) == 0.0 ) then
  !                  Prof%accepted = .false.
  !               end if
  !               if (T_grid%bathyT(i0,j0) < shelf_depth .or.&
  !                    & T_grid%bathyT(1,j0) < shelf_depth .or.&
  !                    & T_grid%bathyT(i0,j0+1) < shelf_depth .or.&
  !                    & T_grid%bathyT(1,j0+1) < shelf_depth ) then
  !                  Prof%accepted = .false.
  !               end if
  !            else if ( i0 /= ieg .and. j0 == jeg ) then
  !               if ( T_grid%mask(i0,j0,1) == 0.0 .or. T_grid%mask(i0+1,j0,1) == 0.0 ) then
  !                  Prof%accepted = .false.
  !               end if
  !               if ( T_grid%bathyT(i0,j0) < shelf_depth .or.&
  !                     & T_grid%bathyT(i0+1,j0) < shelf_depth ) then
  !                  Prof%accepted = .false.
  !               end if
  !            else
  !               if ( T_grid%mask(i0,j0,1) == 0.0 ) then
  !                  Prof%accepted = .false.
  !               end if
  !               if ( T_grid%bathyT(i0,j0) < shelf_depth ) then
  !                  Prof%accepted = .false.
  !               end if
  !            end if
  !         end if ! check surface land-sea mask and depth of ocean

  !         if ( Prof%accepted ) then ! determine vertical position and check mask at depth
  !            allocate(Prof%k_index(Prof%levels))
  !            do k=1, Prof%levels
  !               Prof%k_index(k) = frac_index(Prof%depth(k), (/T_grid%z(i0,j0,:)/))
  !               if ( Prof%k_index(k) < 1.0 ) then
  !                  if ( Prof%depth(k) < T_grid%z(i0,j0,1) ) then
  !                     Prof%k_index(k) = 0.0
  !                  else if ( Prof%depth(k) > T_grid%z(i0,j0,nk) ) then
  !                      Prof%k_index(k) = real(nk)
  !                      Prof%flag(k)= 0.
  !                  end if
  !               end if
  !               if ( k > 3 ) then ! thinning the profile observations to a maximum of 3 within each layer
  !                  if (floor(Prof%k_index(k)) == floor(Prof%k_index(k-3))) then
  !                      Prof%flag(k)=0.0
  !                  end if
  !               end if
  !               if ( Prof%k_index(k) > real(nk) ) then
  !                  call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is greater than nk', FATAL)
  !               else if ( Prof%k_index(k) < 0.0 ) then
  !                  call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is less than 0', FATAL)
  !               end if
  !               k0 = floor(Prof%k_index(k))

  !               if ( k0 >= 1 ) THEN ! snz add
  !                  if ( Prof%flag(k) > 0. ) then ! flag
  !                     if ( i0 /= ieg .and. j0 /= jeg ) then
  !                        if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
  !                             & T_grid%mask(i0+1,j0,k0) == 0.0 .or.&
  !                             & T_grid%mask(i0,j0+1,k0) == 0.0 .or.&
  !                             & T_grid%mask(i0+1,j0+1,k0) == 0.0 ) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     else if ( i0 == ieg .and. j0 /= jeg ) then
  !                        if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
  !                             & T_grid%mask(1,j0,k0) == 0.0 .or.&
  !                             & T_grid%mask(i0,j0+1,k0) == 0.0 .or.&
  !                             & T_grid%mask(1,j0+1,k0) == 0.0) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     else if ( i0 /= ieg .and. j0 == jeg ) then
  !                        if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
  !                             & T_grid%mask(i0+1,j0,k0) == 0.0) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     else
  !                        if ( T_grid%mask(i0,j0,k0) == 0.0 ) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     end if

  !                     if ( i0 /= ieg .and. j0 /= jeg) then
  !                        if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(i0+1,j0,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(i0,j0+1,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(i0+1,j0+1,k0+1) == 0.0 ) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     else if ( i0 == ieg .and. j0 /= jeg ) then
  !                        if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(1,j0,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(i0,j0+1,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(1,j0+1,k0+1) == 0.0) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     else if ( i0 /= ieg .and. j0 == jeg ) then
  !                        if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
  !                             & T_grid%mask(i0+1,j0,k0+1) == 0.0) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     else
  !                        if ( T_grid%mask(i0,j0,k0+1) == 0.0 ) then
  !                           Prof%flag(k) = 0.0
  !                        end if
  !                     end if

  !                     if ( Prof%data(k) == MISSING_VALUE &
  !                             .or. Prof%depth(k) == MISSING_VALUE ) then
  !                        Prof%flag(k) = 0.0
  !                     end if
  !                  end if ! flag
  !               end if ! snz add
  !            end do
  !         end if ! determine vertical position and check mask at depth

  !         if ( Prof%accepted ) then ! calculate forward operator indices and weights
  !           allocate(Prof%obs_def(Prof%levels))
  !           ii = i0; jj = j0
  !           frac_lat = Prof%j_index - jj
  !           frac_lon = Prof%i_index - ii

  !           coef(1) = (1.0 - frac_lon) * (1.0 - frac_lat)
  !           coef(2) = frac_lon * (1.0 - frac_lat)
  !           coef(3) = (1.0 - frac_lon) * frac_lat
  !           coef(4) = frac_lon * frac_lat

  !           if ( ied > ni .and. ii < isd ) ii = ii + ni
  !           if ( isd < 1 .and. ii > ied ) ii = ii - ni

  !           do k=1, Prof%levels
  !             k0 = floor(Prof%k_index(k))
  !             frac_k = Prof%k_index(k) - k0

  !             if ( k0 == 0 ) then
  !               state_index(1) = (jj-jsd)*lon_len + ii-isd + 1
  !               state_index(2) = (jj-jsd)*lon_len + ii-isd + 2
  !               state_index(3) = (jj-jsd+1)*lon_len + ii-isd + 1
  !               state_index(4) = (jj-jsd+1)*lon_len + ii-isd + 2
  !               state_index(5) = state_index(1)
  !               state_index(6) = state_index(2)
  !               state_index(7) = state_index(3)
  !               state_index(8) = state_index(4)
  !             else if (k0 == nk ) then
  !               state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+1
  !               state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+2
  !               state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+1
  !               state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+2
  !               state_index(5) = state_index(1)
  !               state_index(6) = state_index(2)
  !               state_index(7) = state_index(3)
  !               state_index(8) = state_index(4)
  !             else
  !               state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 1
  !               state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 2
  !               state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 1
  !               state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 2
  !               state_index(5) = k0*blk + (jj-jsd)*lon_len + ii-isd + 1
  !               state_index(6) = k0*blk + (jj-jsd)*lon_len + ii-isd + 2
  !               state_index(7) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 1
  !               state_index(8) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 2
  !             end if

  !             if ( frac_lon == 0.0 ) then
  !               state_index(2) = state_index(1)
  !               state_index(4) = state_index(3)
  !               state_index(6) = state_index(5)
  !               state_index(8) = state_index(7)
  !             end if

  !             if ( frac_lat == 0.0 ) then
  !               state_index(3) = state_index(1)
  !               state_index(4) = state_index(2)
  !               state_index(7) = state_index(5)
  !               state_index(8) = state_index(6)
  !             end if

  !             coef(5) = 1.0 - frac_k
  !             coef(6) = frac_k

  !             if ( frac_k == 0.0 ) then
  !               state_index(5) = state_index(1)
  !               state_index(6) = state_index(2)
  !               state_index(7) = state_index(3)
  !               state_index(8) = state_index(4)
  !             end if

  !             call def_forward_operator(8, state_index(1:8), coef(1:6), Prof%obs_def(k))
  !           end do
  !         endif ! calculate forward operator indices and weights

  !         allocate(Prof%next) ! allocate next profile and link it to current one
  !         Prof%next%prev=>Prof
  !         Prof=>Prof%next
  !       end do
  !     end do
  !   end do

  !   call mpp_sync_self()
  !   call mpp_close(unit)
  ! end subroutine open_mooring_dataset

  subroutine open_virtual_mooring_dataset(Profiles, Domain, global_grid, &
                  local_grid, kdroot, filename, time_start, time_end, dtime)
    type(ocean_profile_type), pointer :: Profiles !< This is an unstructured recursive list of profiles
                                                  !< which are either within the localized domain corresponding
                                                  !< to the Domain argument, or the global profile list
    type(domain2d), pointer, intent(in) :: Domain !< MOM grid type for the local domain
    type(grid_type), pointer, intent(in) :: global_grid !< MOM grid type for the local domain
    type(grid_type), pointer, intent(in) :: local_grid  !< MOM grid type for the local domain
    type(kd_root), pointer :: kdroot
    character(len=*), intent(in) :: filename            !< The filename containing profile data
    type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
    integer, intent(in) :: dtime                        !< time interval for virtual moorings (s)

    real :: lon, lat, time
    integer :: nlat, nlon, ndepth, ntime
    integer :: ni, nj, nk
    real :: ri0, rj0
    real, allocatable, dimension(:) :: data, quality
    real, allocatable, dimension(:) :: flag
    type(ocean_profile_type), pointer :: Prof
    integer :: unit, ndim, nvar, natt, max_profiles
    integer :: stdout_unit
    integer :: inst_type, num_vars
    integer :: num_levs, k, t1, kk, i, j, i0, j0, k0, nlevs, a, n
    integer, allocatable, dimension(:) :: var_ids
    integer :: yr, mon, day, hr, min, sec
    integer :: ii, jj
    logical :: data_is_local, localize_data, data_in_period
    character(len=32) :: fldname, axisname, time_units
    character(len=138) :: emsg_local
    type(time_type) :: mooring_time, obs_time
    type(axistype), pointer :: lon_axis, lat_axis, time_axis, depth_axis
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_t, field_s
    integer, allocatable, dimension(:) :: lon1d,lat1d
    real, allocatable, dimension(:) :: lons, lats, depths, times
    real, allocatable, dimension(:,:,:,:) :: mooring_obs
    integer :: mooring_size(4)
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isg, ieg, jsg, jeg, halox, haloy, lon_len, lat_len, blk
    integer :: mooring_count = 0, r_num
    type(horiz_interp_type) :: Interp
    real :: lon_out(1, 1), lat_out(1, 1)
    real :: lat_bound = 90.
    real :: frac_lon, frac_lat, frac_k
    real, dimension(max_neighbors) :: coef, dist
    integer, dimension(2) :: state_size
    integer, dimension(max_neighbors) :: state_index, inds
    integer :: dy, sc,ticks, dy2, sc2, tot_sec, tot_sec2
    real :: sum_coef
    integer :: seconds_per_day=86400


    Prof=>Profiles
    ! Advance to end of profile list
    do while (associated(Prof%next))
      Prof=>Prof%next
    end do
    if (Prof%initialized) then
       allocate(Prof%next)
       call update_profile_ptr(Prof,Prof%next)
    endif
    ni = global_grid%ni; nj = global_grid%nj; nk = global_grid%nk
    num_vars = 2 ! hard-coded for now
    allocate(var_ids(num_vars))
    var_ids(1)=TEMP_ID
    var_ids(2)=SALT_ID
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
    isc = isc - isd + 1 ; iec = iec - isd + 1 ; jsc = jsc - jsd + 1 ; jec = jec - jsd + 1
    ied = ied - isd + 1 ; jed = jed - jsd + 1 ; isd = 1 ; jsd = 1
    lon_len = ied-isd+1
    lat_len = jed-jsd+1
    blk = lat_len*lon_len
    stdout_unit = stdout()
    ni = global_grid%ni; nj = global_grid%nj
    allocate(lon1d(blk))
    allocate(lat1d(blk))
    do i = isd, ied; do j = jsd, jed
      lon1d((j-1)*ied+i) = i
      lat1d((j-1)*ied+i) = j
    enddo; enddo
    inst_type = ODA_MRB
    call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_SINGLE, action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    write (UNIT=stdout_unit, FMT='("Opened mooring dataset: ",A)') trim(filename)
    ! get axis information
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case ( trim(axisname) )
       case ('lon')
          lon_axis => axes(i)
       case ('lat')
          lat_axis => axes(i)
       case ('depth')
          depth_axis => axes(i)
       end select
    end do
    call mpp_get_atts(lon_axis,len=nlon)
    call mpp_get_atts(lat_axis,len=nlat)
    call mpp_get_atts(depth_axis,len=ndepth)
    call get_time(time_start,sc,dy)
    tot_sec=dy*seconds_per_day + sc
    call get_time(time_end,sc,dy)
    tot_sec2=dy*seconds_per_day + sc
    if (dtime<=0) then
       print *,'time interval must be positive'
       return
    endif
    ntime=ceiling(real(tot_sec2-tot_sec)/real(dtime))
    print *,'ntime=',ntime,tot_sec2,tot_sec,dtime
    if (nlon .ne. nlat) call error_mesg('ocean_da_core_mod::open_virtual_mooring_dataset', 'len(lon) not equal to len(lat)', FATAL)
    allocate(lons(nlon), lats(nlat), depths(ndepth), times(ntime))
    allocate(mooring_obs(ntime,ndepth,nlat,nlon))
    allocate(data(ndepth),flag(ndepth)); data(:)=0.0; flag(:)=1.0
    call mpp_get_axis_data(lon_axis, lons)
    call mpp_get_axis_data(lat_axis, lats)
    call mpp_get_axis_data(depth_axis, depths)
    !call mpp_get_axis_data(time_axis, times)
    mooring_count=0
    do j=1, nlat
      lat = lats(j)
      lon = lons(j)
      data_is_local = .false.
      if ( lon .lt. 0.0 ) lon = lon + 360.0
      if ( lon .gt. 360.0 ) lon = lon - 360.0
      call kd_search_nnearest(kdroot, lon, lat, &
             max_neighbors, inds, dist, r_num, .false.)
      data_is_local = within_domain(lon1d(inds(1)), lat1d(inds(1)), &
             isd+1, ied-1, jsd+1, jed-1, ni, nj)
      if (.not. data_is_local) cycle
      i0 = lon1d(inds(1)); j0 = lat1d(inds(1))
      !print *,'Profile location i,j= ',i0,j0
      Prof%accepted=.true.
      Prof%colocated=.true.
      if (local_grid%mask(i0,j0,1) == 0.0  ) then
         print *,'rejecting adjacent to land'
         Prof%accepted = .false.
      end if
      if (local_grid%bathyT(i0,j0) < shelf_depth  ) then
         print *,'rejecting below shelf depth'
         Prof%accepted = .false.
      end if
      if (.not. Prof%accepted) cycle
      !! allocate profile structure content and put in data
      call allocate_virtual_profile(Prof, ndepth, num_vars, var_ids, inst_type, lon, lat, i0, j0, dist(1), depths, data, flag)
      mooring_time=time_start
      state_size=(/lon_len,lat_len/)
      do t1=1, ntime
        data_in_period = .false.
        Prof%filename = trim(filename)
        if ( mooring_time >= time_start .and. mooring_time <= time_end ) data_in_period = .true.
        if ( .not. data_in_period ) cycle
        Prof%time=mooring_time
        call assign_forward_operator(max_neighbors,state_size, lon1d, lat1d, inds, dist, Prof)
        mooring_count = mooring_count + 1
        Prof%initialized=.true.
        if (t1<ntime) then
           allocate(Prof%next) ! allocate next profile and link it to current one
           call copy_virtual_profile(Prof,Prof%next)
           call update_profile_ptr(Prof,Prof%next)
        endif
        mooring_time = increment_time(mooring_time, dtime,0)
      end do
    end do

    call show_profiles(Prof)
    deallocate(lon1d,lat1d)
    call mpp_sync_self()
    call mpp_close(unit)

  contains

    subroutine allocate_virtual_profile(Prof, ndepth, nvar, var_ids, inst_type, lon, lat, i0, j0, dist, depths, data, flag)
      type(ocean_profile_type), pointer :: Prof !< Pointer to profile data structure associated with
                                                !! current profile in the list
      integer :: ndepth                         !< The number of depth levels in the profile
      integer, intent(in)       :: nvar         !< The number of variables associated with the current profile
      integer, dimension(nvar) ::  var_ids      !< list of unique variable identifiers
      integer :: inst_type                      !< instrument type integer identifer
      real    :: lat                            !< Profile latitude
      real    :: lon                            !< Profile longitude
      integer :: i0, j0                         !< nearest indices corresponding to lat,lon arrays
      real    :: dist                           !< distance between profile and nearest grid point
      integer :: n
      real, dimension(ndepth) :: depths         !< profile depth (m)
      real, dimension(ndepth) :: data           !< profile tracer concentration  (degC or psu)
      real, dimension(ndepth) :: flag           !< profile quality flag

      allocate(Prof%depth(nvar,ndepth))
      allocate(Prof%data(nvar,ndepth))
      allocate(Prof%flag(nvar,ndepth))
      Prof%flag(:,:)=1.0
      Prof%num_variables = nvar
      allocate(Prof%var_id(nvar)); Prof%var_id=var_ids
      Prof%inst_type = inst_type
      Prof%levels = ndepth
      Prof%lat = lat; Prof%lon = lon
      Prof%i_index = i0; Prof%j_index = j0
      Prof%nbr_dist = dist
      Prof%basin_mask = local_grid%basin_mask(i0,j0)
      Prof%time_window = time_window(inst_type)
      allocate(Prof%obs_error(nvar))
      Prof%obs_error(1) = temp_error(VIRTUAL_MOORING_FILE)
      Prof%obs_error(2) = salt_error(VIRTUAL_MOORING_FILE)
      Prof%loc_dist = temp_dist(inst_type)
      do k=1, Prof%levels
        do n=1,Prof%num_variables
          Prof%data(n,k)=0.0
          Prof%depth(n,k) = depths(k)
        enddo
      end do
      allocate(Prof%k_index(nvar,ndepth))
      !print *,'ndepth, Prof%depth=',ndepth, Prof%depth(1:Prof%levels)
      !print *,'local_grid%z=',local_grid%z(i0,j0,:)
      do n=1,Prof%num_variables
        do k=1, Prof%levels
          Prof%k_index(n,k) = frac_index(Prof%depth(n,k), (/0.,local_grid%z(i0,j0,:)/)) -1
          if ( Prof%k_index(n,k) < 1.0 ) then
             if ( Prof%depth(n,k) < local_grid%z(i0,j0,1) ) then
                Prof%k_index(n,k) = 1.0
             else if ( Prof%depth(n,k) > local_grid%z(i0,j0,nk) ) then
                Prof%k_index(n,k) = real(nk)
                Prof%flag(n,k)= 0.
             end if
          end if
          if ( Prof%k_index(n,k) > real(nk) ) then
             call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is greater than nk', FATAL)
          else if ( Prof%k_index(n,k) < 0.0 ) then
             call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is less than 0', FATAL)
          end if
          k0 = floor(Prof%k_index(n,k))
          if ( k0 >= 1 ) THEN
             if ( Prof%flag(n,k) > 0. ) then ! flag
                if ( local_grid%mask(i0,j0,k0) == 0.0 .or.&
                   & local_grid%mask(i0,j0,k0+1) == 0.0 ) Prof%flag(n,k) = 0.0
             end if
          end if
        enddo
      enddo

    end subroutine allocate_virtual_profile

    subroutine copy_virtual_profile(Prof,Prof_next)
      type(ocean_profile_type), pointer :: Prof !< Pointer to profile data structure associated with
                                                      !< current profile in the list
      type(ocean_profile_type), pointer :: Prof_next  !< Pointer to profile data structure associated with
                                                      !< next profile in the list

      integer :: ndepth !< number of depth levels in current profile
      integer :: nvar
      if (.not. ASSOCIATED(Prof_next)) then
        call error_mesg('ocean_da_core_mod: copy_virtual_profile', 'next profile not associated!',FATAL)
      endif

      ndepth=Prof%levels
      nvar=Prof%num_variables
      Prof_next%num_variables = Prof%num_variables
      allocate(Prof_next%var_id(Prof%num_variables))
      Prof_next%var_id = Prof%var_id
      Prof_next%levels = Prof%levels
      allocate(Prof_next%depth(nvar,ndepth))
      allocate(Prof_next%data(nvar,ndepth))
      allocate(Prof_next%flag(nvar,ndepth))
      allocate(Prof_next%k_index(nvar,ndepth))
      Prof_next%initialized=.true.
      Prof_next%flag = Prof%flag
      Prof_next%depth = Prof%depth
      Prof_next%time = Prof%time
      Prof_next%data = Prof%data
      Prof_next%inst_type = Prof%inst_type
      Prof_next%levels = ndepth
      Prof_next%k_index = Prof%k_index
      Prof_next%lat = Prof%lat; Prof_next%lon = Prof%lon
      Prof_next%i_index = Prof%i_index; Prof_next%j_index = Prof%j_index
      Prof_next%nbr_dist = Prof%nbr_dist
      Prof_next%basin_mask = Prof%basin_mask
      Prof_next%time_window = Prof%time_window
      allocate(Prof_next%obs_error(nvar))
      Prof_next%obs_error = Prof%obs_error
      Prof_next%loc_dist = Prof%loc_dist
      Prof_next%accepted = Prof%accepted
      Prof_next%colocated = Prof%colocated
      Prof_next%cnext => NULL();Prof_next%cprev => NULL()

    end subroutine copy_virtual_profile

  end subroutine open_virtual_mooring_dataset

  ! get profiles obs relevant to current analysis interval
  subroutine get_profiles(model_time, Profiles, Cprof)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), pointer :: Profiles
    type(ocean_profile_type), pointer :: Cprof

    type(ocean_profile_type), pointer :: Prof
    integer :: nprof
    integer :: i, yr, mon, day, hr, min, sec
    integer :: stdout_unit
    integer :: dy,sc,ticks
    type(ocean_profile_type), pointer :: PREVIOUS=>NULL()
    type(time_type) :: tdiff

    nprof = 0
    stdout_unit = stdout()
    if (.not. associated(Profiles)) return
    Prof=>Profiles
    do while (associated(Prof%prev))
      Prof=>Prof%prev
    enddo

    write (UNIT=stdout_unit, FMT='("Gathering profiles for current analysis time")')
    call get_date(model_time, yr, mon, day, hr, min, sec)
    write (UNIT=stdout_unit, FMT='("Current YYYY/MM/DD = ",I4,"/",I2,"/",I2,",",I2,":",I2)') yr, mon, day, hr, min

    Cprof=>NULL()
    nprof=0
    do while(associated(Prof))
       if (.not. Prof%initialized) exit
       if ( Prof%time <= model_time ) then
          tdiff = model_time - Prof%time
       else
          tdiff = Prof%time - model_time
       end if

       !call get_time(Prof%time_window,sc,dy,ticks)
       !print *,'in get profiles: Profile time window= ',dy,sc,ticks
       !call get_time(Prof%time,sc,dy,ticks)
       !print *,'in get_profiles: profile time = ',dy,sc,ticks
       !call get_time(model_time,sc,dy,ticks)
       !print *,'in get_profiles: Model time= ',dy,sc,ticks
       !call get_time(tdiff,sc,dy,ticks)
       !print *,'in get_profiles: Model time diff= ',dy,sc,ticks
       ! no tdiff criteria for monthly mean data like
       ! but tdiff criteria has to be set for daily data
       if ( tdiff <= Prof%time_window .and. Prof%accepted ) then
          !current_type(Prof%inst_type) = current_type(Prof%inst_type) + 1
          Prof%tdiff = tdiff
          if (.not.associated(CProf)) then
             CProf=>Prof
             Cprof%cprev=>NULL()
          else
             Cprof%cnext=>Prof
             PREVIOUS=>Cprof
             Cprof=>Prof
             Cprof%cprev=>PREVIOUS
          endif
          nprof=nprof+1
       end if
       if (associated(Prof%next)) then
          Prof=>Prof%next
       else
          Prof=>NULL()
       endif
     end do
     print *,'Total number of profiles = ',nprof


     return

  end subroutine get_profiles

  subroutine show_profiles(Profiles,current)
    type(ocean_profile_type), pointer :: Profiles
    logical, optional, intent(in) :: current
    type(ocean_profile_type), pointer :: Prof
    logical :: cur
    integer :: sc,dy
    cur=.false.
    if (present(current)) cur=current
    if (.not. associated(Profiles)) return

    Prof=>Profiles

    if (cur) then
       do while (associated(Prof%cprev))
         Prof=>Prof%cprev
       enddo
    else
       do while (associated(Prof%prev))
         Prof=>Prof%prev
       enddo
    endif

    do while (associated(Prof) )
      if (Prof%initialized) then
         call get_time(Prof%time,sc,dy)
        print *,'Profile (day,sec,levels,lon,lat)=,',dy,sc,Prof%levels,Prof%lon,Prof%lat
        if (cur) then
           Prof=>Prof%cnext
        else
           Prof=>Prof%next
        endif
     else
        exit
     endif
    enddo
  end subroutine show_profiles

!=======================================================================
! Assigns a forward operator type to an observation
  subroutine assign_forward_operator(r_num, state_size, ilons, jlats, inds, dist, Prof)
    integer, intent(in) :: r_num  !< the number of neighbor points
    integer, dimension(2), intent(in) :: state_size
    integer, dimension(state_size(1)*state_size(2)), intent(in) :: ilons
    integer, dimension(state_size(1)*state_size(2)), intent(in) :: jlats
    integer, dimension(r_num), intent(in) :: inds
    real, dimension(r_num), intent(in) :: dist
    type(ocean_profile_type), pointer :: Prof

    real, dimension(r_num) :: coef
    real :: sum_coef
    integer :: k
    integer, dimension(r_num) :: i_index
    integer, dimension(r_num) :: j_index

    sum_coef=0.0
    do k=1,r_num
      i_index(k) = ilons(inds(k))
      j_index(k) = jlats(inds(k))
      coef(k)=1.0/(dist(k)+small_xy)
      sum_coef=sum_coef+coef(k)
    enddo
    do k=1,r_num
      coef(k)=coef(k)/sum_coef
    enddo
    allocate(Prof%obs_def)
    call def_forward_operator(r_num, state_size, i_index, j_index, coef, Prof%obs_def)

  end subroutine assign_forward_operator

  !=======================================================================
  ! Puts FO interpolation coefficients and state variable indices into an FO_type data structure.
  subroutine def_forward_operator(num_state, state_size, i_ind, j_ind, coef, obs_def)
    integer, intent(in) :: num_state
    integer, dimension(2), intent(in) :: state_size
    integer, dimension(num_state), intent(in) :: i_ind, j_ind
    real, dimension(num_state), intent(in) :: coef
    type(forward_operator_type), intent(inout) :: obs_def

    integer :: i

    ! Set aside storage for defining this ob
    obs_def%num = num_state
    allocate(obs_def%coef(num_state))
    allocate(obs_def%i_index(num_state), obs_def%j_index(num_state))
    obs_def%state_size=state_size
    ! Load the state variable index and coefficient for each state variable
    do i = 1, num_state
      obs_def%i_index(i) = i_ind(i)
      obs_def%j_index(i) = j_ind(i)
    end do

    do i=1, num_state
       obs_def%coef(i) = coef(i)
    end do
  end subroutine def_forward_operator

  subroutine apply_forward_operator(Prof, ocean_CS)
    type(ocean_profile_type), pointer :: Prof
    type(ocean_control_struct), pointer :: ocean_CS
    logical :: have_temp, have_salt, have_u, have_v

    integer :: m,n, p, k
    integer :: i_index, j_index, k0, k1
    integer :: itemp, isalt
    real :: frac_cell


    if (.not. ASSOCIATED(Prof)) then
       print *,'apply_forward_operator: profile not associated!'
       return
    endif

    have_temp=.false.;have_salt=.false.
    have_u=.false.;have_v=.false.

    do n=1, Prof%num_variables
      if (Prof%var_id(n) == TEMP_ID) then
         have_temp=.true.
         itemp=n
      else if (Prof%var_id(n) == SALT_ID) then
         have_salt=.true.
         isalt=n
      endif
    enddo

    Prof%ensemble_size=ocean_CS%ensemble_size
    if (associated(ocean_CS%T)) then
       if ( .not. associated(Prof%forecast)) then
          allocate(Prof%forecast(Prof%ensemble_size,Prof%num_variables,Prof%levels))
          Prof%forecast(:,:,:)=0.0
       endif
    endif


    if (have_temp) then
       do p = 1, Prof%ensemble_size
         do n = 1, Prof%obs_def%num
           i_index=Prof%obs_def%i_index(n)
           j_index=Prof%obs_def%j_index(n)
           do k=1,Prof%levels
             k0=floor(Prof%k_index(itemp,k)) ! use linear interpolation in the vertical
             k1=ceiling(Prof%k_index(itemp,k))
             frac_cell=Prof%k_index(itemp,k)-k0
             if (frac_cell<0. .or. frac_cell>=1.0) then
                call mpp_error(FATAL,'fractional distance outside bounds [0,1)')
             endif
             Prof%forecast(p,itemp,k)=Prof%forecast(p,itemp,k)+Ocean_CS%T(i_index,j_index,k0,p)*Prof%obs_def%coef(n)*(1.-frac_cell)&
                    + Ocean_CS%T(i_index,j_index,k1,p)*Prof%obs_def%coef(n)*(frac_cell)
           enddo
         enddo
       enddo
    endif

    if (have_salt) then
       do p = 1, Prof%ensemble_size
         do n = 1, Prof%obs_def%num
           i_index=Prof%obs_def%i_index(n)
           j_index=Prof%obs_def%j_index(n)
           do k=1,Prof%levels
             k0=floor(Prof%k_index(isalt,k)) ! use linear interpolation in the vertical
             k1=ceiling(Prof%k_index(isalt,k))
             frac_cell=Prof%k_index(isalt,k)-k0
             if (frac_cell<0. .or. frac_cell>=1.0) then
                call mpp_error(FATAL,'fractional distance outside bounds [0,1)')
             endif
             Prof%forecast(p,isalt,k)=Prof%forecast(p,isalt,k)+Ocean_CS%S(i_index,j_index,k0,p)*Prof%obs_def%coef(n)*(1.-frac_cell)&
                    + Ocean_CS%S(i_index,j_index,k1,p)*Prof%obs_def%coef(n)*(frac_cell)
           enddo
         enddo
       enddo
    endif

  end subroutine apply_forward_operator


  subroutine update_profile_ptr(Prof,Next)
    type(ocean_profile_type), pointer :: Prof
    type(ocean_profile_type), pointer :: Next

    Prof%next%prev=>Prof
    Prof=>Next
    Prof%initialized = .true.

  end subroutine update_profile_ptr

end module ocean_da_core_mod
