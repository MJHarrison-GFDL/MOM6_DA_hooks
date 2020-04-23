module write_ocean_obs_mod

use mpp_io_mod, only : fieldtype, axistype, mpp_open, mpp_write_meta, mpp_write, mpp_close
use mpp_io_mod, only : MPP_OVERWR, MPP_APPEND, MPP_WRONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
use mpp_mod, only : mpp_error, FATAL
use ocean_da_types_mod, only : MISSING_VALUE, TEMP_ID, SALT_ID, U_ID, V_ID, SSH_ID
use ocean_da_types_mod, only : ocean_profile_type, MAX_LEVELS_FILE
use time_manager_mod, only : time_type, get_time, set_date, operator ( - )

implicit none
private


integer, parameter :: ref_yr=1900, ref_mon=1, ref_day=1, &
        ref_hr=0, ref_min=0, ref_sec=0,max_files=1000

integer :: ref_seconds, ref_days, chid, wmo_id
!integer, save :: station_num(max_files), unit_num(max_files), nfiles
type(time_type) :: ref_time, time
logical :: module_is_initialized=.false.

public :: open_profile_file, write_profile, close_profile_file, write_ocean_obs_init, &
          profile_name_generator

type, private :: profile_file_type
   integer :: unit
   logical :: have_temp
   type(axistype) :: time_axis
   type(fieldtype) :: depth_field
   type(fieldtype) :: lon_field
   type(fieldtype) :: lat_field
   type(fieldtype) :: time_field
   type(fieldtype) :: temp_err_field
   type(fieldtype) :: temp_field
   type(fieldtype) :: temp_flag_field
   type(fieldtype) :: temp_forecast_field
   type(fieldtype)  :: temp_analysis_field
   logical :: have_salt
   type(fieldtype) :: salt_err_field
   type(fieldtype) :: salt_field
   type(fieldtype) :: salt_flag_field
   type(fieldtype)  :: salt_forecast_field
   type(fieldtype)  :: salt_analysis_field
   logical :: have_u
   type(fieldtype) :: u_err_field
   type(fieldtype) :: u_field
   type(fieldtype) :: u_flag_field
   type(fieldtype)  :: u_forecast_field
   type(fieldtype)  :: u_analysis_field
   logical :: have_v
   type(fieldtype) :: v_err_field
   type(fieldtype) :: v_field
   type(fieldtype) :: v_flag_field
   type(fieldtype)  :: v_forecast_field
   type(fieldtype)  :: v_analysis_field
end type profile_file_type
#include <netcdf.inc>

contains

function open_profile_file(Prof, name, include_posterior)

  type(ocean_profile_type), intent(in) :: Prof
  character(len=*), intent(in) :: name !< The output file name (string).
  logical, optional, intent(in) :: include_posterior !< Write the posterior corrections as well as the priors
  type(profile_file_type) :: open_profile_file

  integer :: i, unit, m, station
  integer :: threading, fileset
  character(len=128) :: units, time_units
  real, dimension(1) :: arr !< temporary storage for station axis (hard-coded to 1 station per file)
                            !! stations share the same time axis.
  real, dimension(Prof%levels) :: array
  real, dimension(Prof%ensemble_size) :: array2
  type(axistype) :: depth_axis,ensemble_axis
  type(axistype) :: station_axis
  logical :: have_temp, have_salt, have_u, have_v
  logical :: include_post
  type(fieldtype) :: lon_field, lat_field, time_field, depth_field

  threading=MPP_MULTI
  fileset=MPP_SINGLE

  have_temp=.false.;have_salt=.false.;have_u=.false.;have_v=.false.

  include_post=.true.
  if (present(include_posterior)) include_post=include_posterior

  do i=1,Prof%num_variables
    if (Prof%var_id(i)==TEMP_ID) have_temp=.true.
    if (Prof%var_id(i)==SALT_ID) have_salt=.true.
    if (Prof%var_id(i)==U_ID) have_u=.true.
    if (Prof%var_id(i)==V_ID) have_v=.true.
  enddo

  ref_time = set_date(ref_yr, ref_mon, ref_day, ref_hr, ref_min, ref_sec)
  call get_time(ref_time, ref_seconds, ref_days)
  print *,'Calling mpp_open'
  call mpp_open(unit, trim(name), action=MPP_OVERWR, form=MPP_NETCDF,&
                threading=threading, fileset=fileset)

  open_profile_file%unit = unit


  !call mpp_write_meta(unit,depth_axis,'depth_index','none','depth index',&
  !                  cartesian='Z',sense=-1)!,data=(/(float(i),i=1,MAX_LEVELS_FILE)/))
  !pgf90 complains about the above. This is a compiler bug. Workaround:
  ! currently, variables associated with a profile share a common depth axis
  array = Prof%depth(1,:)
  call mpp_write_meta(unit,depth_axis,'depth','m','depth',&
       cartesian='Z',sense=-1,data=array)

  array2 = (/(float(i),i=1,Prof%ensemble_size)/)
  call mpp_write_meta(unit,ensemble_axis,'ensemble','none','ensemble index',&
                      cartesian='N',sense=1,data=array2)

  arr(1)=1.0 ! hard-coding for now (may want to have groups of simultaneous profiles for analysis purposes)
  call mpp_write_meta(unit,station_axis,'station','none','station index',&
                      cartesian='N',sense=1,data=arr)

  call mpp_write_meta(unit,lon_field,(/station_axis/),&
          'longitude','degrees_E','longitude',min=-301.0,max=61.0)
  call mpp_write_meta(unit,lat_field,(/station_axis/),&
          'latitude','degrees_N','latitude',min=-91.0,max=91.0)

  call mpp_write_meta(unit,open_profile_file%time_axis,'ntime','none',&
                      'time index', cartesian='T',sense=1)

  write(time_units,'(a,i4.4,a,i2.2,a,i2.2,a)')  'days since ',ref_yr,'-',ref_mon,'-',ref_day,' 00:00:00'

  call mpp_write_meta(unit,open_profile_file%time_field,(/open_profile_file%time_axis/),'time',trim(time_units),'time')

  open_profile_file%have_temp=.false.
  open_profile_file%have_salt=.false.
  open_profile_file%have_u=.false.
  open_profile_file%have_v=.false.

  if (have_temp) then
     open_profile_file%have_temp=.true.
     units='degC'
     print *,'Writing temp metadata'
     call mpp_write_meta(unit,open_profile_file%temp_err_field,(/station_axis,open_profile_file%time_axis/),&
          'thetao_obs_error',trim(units),'obs potential temperature error',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%temp_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'thetao_obs',trim(units),'obs potential temperature',min=-10.0,max=50.0,missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%temp_flag_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'thetao_obs_flag','none','observation level flag',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%temp_forecast_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
          'thetao_prior',trim(units),'model first-guess potential temperature',min=-10.0,max=50.0,missing=MISSING_VALUE)
     if (include_post) then
        call mpp_write_meta(unit,open_profile_file%temp_analysis_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
             'thetao_posterior',trim(units),'analysis potential temperature',min=-10.0,max=50.0,missing=MISSING_VALUE)
     endif
  endif

  if (have_salt) then
     open_profile_file%have_salt=.true.
     units='psu'
     print *,'Writing temp metadata'
     call mpp_write_meta(unit,open_profile_file%salt_err_field,(/station_axis,open_profile_file%time_axis/),&
          'so_obs_error',trim(units),'obs salinity error',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%salt_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'so_obs',trim(units),'obs salinity',min=-10.0,max=50.0,missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%salt_flag_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'so_obs_flag','none','observation level flag',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%salt_forecast_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
          'so_prior',trim(units),'model first-guess salinity',min=-10.0,max=50.0,missing=MISSING_VALUE)
     if (include_post) then
        call mpp_write_meta(unit,open_profile_file%salt_analysis_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
             'so_posterior',trim(units),'analysis salinity',min=-10.0,max=50.0,missing=MISSING_VALUE)
     endif
  endif

  if (have_u) then
     open_profile_file%have_u=.true.
     units='m s-1'
     call mpp_write_meta(unit,open_profile_file%u_err_field,(/station_axis,open_profile_file%time_axis/),&
          'uo_obs_error',trim(units),'obs zonal velocity error',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%u_flag_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'uo_obs_flag','none','observation level flag',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%u_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'uo_obs',trim(units),'obs zonal velocity',min=-10.0,max=50.0,missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%u_forecast_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
          'uo_prior',trim(units),'model first-guess zonal velocity',min=-10.0,max=10.0,missing=MISSING_VALUE)
     if (include_post) then
        call mpp_write_meta(unit,open_profile_file%u_analysis_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
             'uo_posterior',trim(units),'analysis zonal velocity',min=-10.0,max=10.0,missing=MISSING_VALUE)
     endif
  endif

  if (have_v) then
     open_profile_file%have_v=.true.
     units='m s-1'
     call mpp_write_meta(unit,open_profile_file%v_err_field,(/station_axis,open_profile_file%time_axis/),&
          'vo_obs_error',trim(units),'obs meridional velocity error',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%v_flag_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'vo_obs_flag','none','observation level flag',missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%v_field,(/station_axis,depth_axis,open_profile_file%time_axis/),&
          'vo_obs',trim(units),'obs meridional velocity',min=-10.0,max=50.0,missing=MISSING_VALUE)
     call mpp_write_meta(unit,open_profile_file%v_forecast_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
          'vo_prior',trim(units),'model first-guess meridional velocity',min=-10.0,max=10.0,missing=MISSING_VALUE)
     if (include_post) then
        call mpp_write_meta(unit,open_profile_file%v_analysis_field,(/ensemble_axis,station_axis,depth_axis,open_profile_file%time_axis/),&
             'vo_posterior',trim(units),'analysis meridional velocity',min=-10.0,max=10.0,missing=MISSING_VALUE)
     endif
  endif

!  print *,'Writing depth metadata'
!  call mpp_write_meta(unit,depth_field,(/depth_axis/),&
!          'depth','meters','depth of obs',min=0.0,max=7000.0,missing=MISSING_VALUE)
!  print *,'Writing ensemble'
  !     print *,'Writing temp metadata'
  call mpp_write(unit, depth_axis)
  call mpp_write(unit, ensemble_axis)
  call mpp_write(unit, station_axis)
  ! only one station currently
  station=1
  call mpp_write(unit, lon_field,Prof%lon,real(station))
  call mpp_write(unit, lat_field,Prof%lat,real(station))
!  call mpp_write()
!  call mpp_write(prof_file%unit,prof_file%lat_field,profile%lat,station)
!  call mpp_write(prof_file%unit,prof_file%depth_field,depth,station)

!  call mpp_write(prof_file%unit,prof_file%time_field,days_since,station)


end function open_profile_file

function profile_name_generator(profile)
  type(ocean_profile_type), pointer :: profile
  character(len=128) :: profile_name_generator
  integer :: dy,sec,len_fnam
  character(len=128) fnam

  call get_time(profile%time,sec,dy)
  len_fnam=len(trim(profile%filename))
  if (profile%filename(len_fnam-1:len_fnam)=='nc') then
     fnam=profile%filename(1:len_fnam-3)
  else
     fnam=profile%filename(1:len_fnam)
  endif
  write(profile_name_generator,'(a,".",I6.6,".",I6.6,".nc")') trim(fnam),dy,sec
  return
end function profile_name_generator

subroutine write_profile(profile)
  type(ocean_profile_type), pointer:: profile

  real, allocatable, dimension(:) :: depth
  real, allocatable, dimension(:,:) :: obs_data
  real, allocatable, dimension(:,:,:) :: forecast, anal
  integer :: levels, secs, days, i, j, n, m, k
  real :: days_since, station, nt
  real, allocatable, dimension(:,:) :: flag
  integer :: findex
  logical :: debug=.false.
  logical :: f_avail=.false. , a_avail=.false.
  integer :: unit
  integer :: temp_index, salt_index
  character(len=128) :: fname
  type(profile_file_type) :: prof_file

  if (.not. associated(profile)) return

  fname = profile_name_generator(Profile)
  print *,'Opening profile file ',fname
  prof_file = open_profile_file(Profile, fname)
  print *, 'profile opened ',prof_file%unit
  !  station_num(findex)=station_num(findex)+1
  !  station=station_num(findex)
  if (allocated(obs_data)) deallocate(obs_data)
  allocate(obs_data(profile%num_variables,profile%levels))
  if (allocated(depth)) deallocate(depth)
  allocate(depth(profile%levels))
  if (allocated(forecast)) deallocate(forecast)
  allocate(forecast(profile%ensemble_size,profile%num_variables,profile%levels))
  if (allocated(anal)) deallocate(anal)
  allocate(anal(profile%ensemble_size,profile%num_variables,profile%levels))
  if (allocated(flag)) deallocate(flag)
  allocate(flag(profile%num_variables,profile%levels))
  if(associated(profile%forecast)) f_avail = .true.
  if(associated(profile%analysis)) a_avail = .true.
  !print *,'Forecase/analysis avail = ',f_avail, a_avail
  levels = profile%levels
  obs_data=MISSING_VALUE; depth=MISSING_VALUE
  forecast=MISSING_VALUE; anal=MISSING_VALUE
  flag=MISSING_VALUE

  if (profile%colocated) then
     depth(1:levels)=profile%depth(1,1:levels)
  else
     print *,'non-colocated profile data not implemented yet'
     call abort()
  endif

  time = profile%time - ref_time
  call get_time(time, secs, days)
  days_since = days + secs/86400.
  !! Hard-coding these variables here since currently writing only one record per file
  !! with a single station
  station=1
  nt=1
  call mpp_write(prof_file%unit,prof_file%time_field,days_since,nt)



  do m=1,profile%ensemble_size
    do n=1,profile%num_variables
      if(f_avail) forecast(m,n,1:levels)=profile%forecast(m,n,1:levels)
      if(a_avail) anal(m,n,1:levels)=profile%analysis(m,n,1:levels)
    enddo
  enddo

  do n=1,profile%num_variables
    flag(n,1:levels)=profile%flag(n,1:levels)
    obs_data(n,1:levels)=profile%data(n,1:levels)
    if (profile%var_id(n) == TEMP_ID) temp_index=n
    if (profile%var_id(n) == SALT_ID) salt_index=n
  enddo


  if (prof_file%have_temp) then
     call mpp_write(prof_file%unit,prof_file%temp_field,obs_data(temp_index,:),nt)
     call mpp_write(prof_file%unit,prof_file%temp_err_field,profile%obs_error(temp_index),nt)
     if(f_avail) then
        call mpp_write(prof_file%unit,prof_file%temp_forecast_field,forecast(:,temp_index,:),nt)
     endif
     if(a_avail) call mpp_write(prof_file%unit,prof_file%temp_analysis_field,anal(:,temp_index,:),nt)
  endif
  if (prof_file%have_salt) then
     call mpp_write(prof_file%unit,prof_file%salt_field,obs_data(salt_index,:),nt)
     call mpp_write(prof_file%unit,prof_file%salt_err_field,profile%obs_error(salt_index),nt)
     if(f_avail) call mpp_write(prof_file%unit,prof_file%salt_forecast_field,forecast(:,salt_index,:),nt)
     if(a_avail) call mpp_write(prof_file%unit,prof_file%salt_analysis_field,anal(:,salt_index,:),nt)
  endif



  call close_profile_file(prof_file%unit)

end subroutine write_profile

subroutine close_profile_file(unit)

  integer, intent(in) :: unit

  call mpp_close(unit)

end subroutine close_profile_file

subroutine write_ocean_obs_init()

  module_is_initialized=.true.


  return

end subroutine write_ocean_obs_init

end module write_ocean_obs_mod
