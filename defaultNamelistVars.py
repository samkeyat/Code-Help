import datetime
from collections import namedtuple


######################### Default namelist.wps variables #######################
# when generating, main script will automatically append commas to end of every line

# &share
wrf_core		= "ARW"
max_dom			= "3"
#date form = YYYY-MM-DD_HH:MM:SS
start_date		= ''
end_date		= ''
interval_seconds	= 21600
io_form_geogrid		= 2

# &geogrid

parent_id 		= "  0,   1,   1"
parent_grid_ratio 	= "  1,   4,   3"
i_parent_start		= "  1,  33,  64"
j_parent_start		= "  1,  45,  31"
e_we			= "130,  61, 145"
e_sn			= "110,  49, 121"
geog_data_res		= "default"
dx			= 30000
dy			= 30000
map_proj		= "lambert"
ref_lat			=   0.5
ref_long		= -84.0
truelat1		= - 8.0
truelat2		= - 2.0
stand_lon		= -84.0
geog_data_path		= "\'/cm/shared/uniol/software/WRFGeodata/3.8/geog\'" # appended single quotes to this line

# &ungrib

out_format		= "\'WPS\'"
prefix			= "\'FILE\'"

# &metgrid
fg_name			= "\'FILE\'"
io_form_metgrid		= 2


##################### Default namelist.input variables #########################

#&time_control
run_days                            = 0
run_hours                           = 6
run_minutes                         = 0
run_seconds                         = 0

start_year                          = '' # "2017, 2017, 2017"
start_month                         = '' # "  01,   01,   01"
start_day                           = '' # "  01,   01,   01"
start_hour                          = '' # "  12,   12,   12"
start_minute                        = '' # "  00,   00,   00"
start_second                        = '' # "  00,   00,   00"

end_year                            = '' # "2017, 2017, 2017"
end_month                           = '' # "  01,   01,   01"
end_day                             = '' # "  01,   01,   01"
end_hour                            = '' # "  18,   18,   18"
end_minute                          = '' # "  00,   00,   00"
end_second                          = '' # "  00,   00,   00"
#interval_seconds                    = 21600
input_from_file                     = "\'.true.\',\'.true.\',\'.true.\'" 
history_interval                    = "60,  60,   60"
frames_per_outfile                  = "1000, 1000, 1000"
restart                             = "\'.false.\'"
restart_interval                    = 5000
io_form_history                     = 2
io_form_restart                     = 2
io_form_input                       = 2
io_form_boundary                    = 2
debug_level                         = 0
#/

#&domains
time_step                           = 90
time_step_fract_num                 = 0
time_step_fract_den                 = 1
max_dom                             = 3
#e_we                                = 130,   61,    145,
#e_sn                                = 110,   49,    121,
#same as in namelist.wps
e_vert                              = "30,    30,    30"
p_top_requested                     = 5000
num_metgrid_levels                  = 33
num_metgrid_soil_levels             = 4
dx                                  = "30000, 7500,  10000"
dy                                  = "30000, 7500,  10000"
grid_id                             = "1,     2,     3"
#parent_id                           = 0,     1,     1,
#i_parent_start                      = 1,     33,    64,
#j_parent_start                      = 1,     45,    31,
#parent_grid_ratio                   = 1,     4,     3,
parent_time_step_ratio              = "1,     1,     1"
feedback                            = 1
smooth_option                       = 0
#/

#&physics
mp_physics                          = "3,     3,     3"
ra_lw_physics                       = "1,     1,     1"
ra_sw_physics                       = "1,     1,     1"
radt                                = "30,    30,    30"
sf_sfclay_physics                   = "1,     1,     1"
sf_surface_physics                  = "2,     2,     2"
bl_pbl_physics                      = "1,     1,     1"
bldt                                = "0,     0,     0"
cu_physics                          = "1,     1,     0"
cudt                                = "5,     5,     5"
isfflx                              = 1
ifsnow                              = 1
icloud                              = 1
surface_input_source                = 1
num_soil_layers                     = 4
num_land_cat                        = 21
sf_urban_physics                    = "0,     0,     0"
#/

#&fdda
#/

#&dynamics
w_damping                           = 0
diff_opt                            = "1,      1,      1"
km_opt                              = "4,      4,      4"
diff_6th_opt                        = "0,      0,      0"
diff_6th_factor                     = "0.12,   0.12,   0.12"
base_temp                           = 290.
damp_opt                            = 0
zdamp                               = "5000.,  5000.,  5000."
dampcoef                            = "0.2,    0.2,    0.2"
khdif                               = "0,      0,      0"
kvdif                               = "0,      0,      0"
non_hydrostatic                     = "\'.true.\',\'.true.\',\'.true.\'"
moist_adv_opt                       = "1,      1,      1"
scalar_adv_opt                      = "1,      1,      1"
#/

#&bdy_control
spec_bdy_width                      = 5
spec_zone                           = 1
relax_zone                          = 4
specified                           = "\'.true.\',\'.false.\',\'.false.\'"
nested                              = "\'.false.\',\'.true.\',\'.true.\'"
#/

#&grib2
#/

#&namelist_quilt
nio_tasks_per_group = 0
nio_groups = 1
#/


############################ RUN_TIME FUNCTIONS ################################

#def run_duration(timedelta duration):
#	r = namedtuple('run_duration',['days','hours','minutes','seconds'])
#	r.days = duration.days
#	r.hours = duration.seconds//3600
#	r.minutes = (duration.seconds%3600)//60
#	r.seconds = (duration.seconds%3600)%60
#	return r

#def run_days(run_duration):
#	run_days = run_duration.days
#	print("run_days: ",str(run_days))
#	return run_days

#def run_hours(run_duration):
#	run_hours = run_duration.days
#	print("run_hours: ",str(run_hours))
#	return run_hours

#def run_minutes(run_duration):
#	run_minutes = int(end_min) - int(start_min)
#	print("run_minutes: ",str(run_minutes))
#	return run_minutes

#def run_seconds(run_duration):
#	run_seconds = int(end_sec) - int(start_sec)
#	print("run_seconds: ",str(run_seconds))
#	return run_seconds

########################## START_DATE FUNCTIONS ################################

def start_date(YYYY,MM,DD,HH,Min,SS):	
	start_date = datetime.datetime(int(YYYY), int(MM), int(DD), int(HH), int(Min), int(SS))
	#t = start_date.strftime('%Y-%m-%d_%H:%M:%S')
	#print("start_date: ",t)
	return start_date

def start_year(YYYY):
	start_year = str(YYYY)+', '+str(YYYY)+', '+str(YYYY)
	print("start_year: ",start_year)
	return start_year

def start_month(MM):
	start_month = '  '+str(MM)+',   '+str(MM)+',   '+str(MM)
	print("start_month: ",start_month)
	return start_month

def start_day(DD):
	start_day = '  '+str(DD)+',   '+str(DD)+',   '+str(DD)
	print("start_day: ", start_day)
	return start_day

def start_hour(HH):
	start_hour = '  '+str(HH)+',   '+str(HH)+',   '+str(HH)
	print("start_hour: ", start_hour)
	return start_hour

def start_minute(mm):
	start_minute = '  '+str(mm)+',   '+str(mm)+',   '+str(mm)
	print("start_minute: ",start_minute)
	return start_minute

def start_second(SS):
	start_second = '  '+str(SS)+',   '+str(SS)+',   '+str(SS)
	print("start_second: ",start_second)
	return start_second


############################ END_DATE FUNCTIONS ################################

def end_date(YYYY,MM,DD,HH,Min,SS):
	end_date = datetime.datetime(int(YYYY), int(MM), int(DD), int(HH), int(Min), int(SS))
	#t = end_date.strftime('%Y-%m-%d_%H:%M:%S')
	#print("end_date: ",t)
	return end_date

def end_year(YYYY):
	end_year = str(YYYY)+', '+str(YYYY)+', '+str(YYYY)
	print("end_year: ",end_year)
	return end_year

def end_month(MM):
	end_month = '  '+str(MM)+',   '+str(MM)+',   '+str(MM)
	print("end_month: ",end_month)
	return end_month

def end_day(DD):
	end_day = '  '+str(DD)+',   '+str(DD)+',   '+str(DD)
	print("end_day: ", end_day)
	return end_day

def end_hour(HH):
	end_hour = '  '+str(HH)+',   '+str(HH)+',   '+str(HH)
	print("end_hour: ", end_hour)
	return end_hour

def end_minute(mm):
	end_minute = '  '+str(mm)+',   '+str(mm)+',   '+str(mm)
	print("end_minute: ",end_minute)
	return end_minute

def end_second(SS):
	end_second = '  '+str(SS)+',   '+str(SS)+',   '+str(SS)
	print("end_second: ",end_second)
	return end_second

def interval_seconds(r_days,r_hours,r_minutes,r_seconds):
	interval_seconds = int(r_days)*86400 + int(r_hours)*3600 + int(r_minutes)*60 + int(r_seconds)
	print("interval_seconds: ",str(interval_seconds))
	return interval_seconds

