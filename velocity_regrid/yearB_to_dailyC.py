import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import sys
import os

half_len_data = nc.Dataset('CM2p6_half_len.nc')
lat_T = half_len_data['grid_y_T'][0:2107]
lon_T = half_len_data['grid_x_T'][...]

vel_sample = nc.Dataset('surf_vel_sample.nc')
lat_U = vel_sample['yu_ocean'][0:2107]
lon_U = vel_sample['xu_ocean'][...]

xu_c = np.zeros(np.shape(lon_U))
yu_c = np.zeros(np.shape(lat_U))

xv_c = np.zeros(np.shape(lon_U))
yv_c = np.zeros(np.shape(lat_U))

#fills the C grid U data  
for i in range(len(lon_U)):
	xu_c[i] = lon_U[i-1] 
for j in range(len(lat_U)):
	yu_c[j] = lat_T[j]
        
#creates C grid V points everywhere above the southernmost Tcell row  
for i in range(len(lon_U)):
	xv_c[i] = lon_T[i] 
for j in range(1,len(lat_U)):
	yv_c[j] = lat_U[j-1]
        
#adds in a row of C grid V points beneath the existing Tcell row 
yv_c[0] = lat_U[0] + (lat_U[0] - lat_U[1])

#the regrid loop 
path2file = ''
filename = ''
path = path2file+filename
vel_data = nc.Dataset(path)
time = vel_data['time'][:]

for tind in range(len(time)):
	u_B = vel_data['usurf'][tind, 0,0:2107,:] #omits polar region
	v_B = vel_data['vsurf'][tind, 0,0:2107,:] #omits polar region 

	u_B = np.transpose(u_B) # convert CM2.6 output from (lat, lon) indexing to (lon, lat) 
	v_B = np.transpose(v_B) # convert CM2.6 output from (lat, lon) indexing to (lon, lat)

	u_B = ma.filled(u_B, fill_value=0.0) #change fill value from -1e20 to 0.0 
	v_B = ma.filled(v_B, fill_value=0.0) #change fill value from -1e20 to 0.0 

	du_Eface_N = half_len_data['ds_21_22_T'][0:2107, :] #omits polar region
	du_Eface_S = half_len_data['ds_20_21_T'][0:2107, :] #omits polar region 

	du_Eface_N = np.transpose(du_Eface_N) 
	du_Eface_S = np.transpose(du_Eface_S) 

	u_C = np.zeros(np.shape(u_B))
	v_C = np.zeros(np.shape(v_B))

	LON_POINTS = 3600
	LAT_POINTS = 2107

	for i in range(LON_POINTS):
	    for j in range(LAT_POINTS):
	        v_C[i,j] = (v_B[i-1, j-1] + v_B[i, j-1])/2
	        
	for i in range(LON_POINTS):
	    for j in range(1,LAT_POINTS):
	        u_C[i,j] = ( (du_Eface_N[i-1, j]*u_B[i-1, j] + du_Eface_S[i-1, j]*u_B[i-1, j-1]) \
	                    /(du_Eface_N[i, j-1] + du_Eface_S[i, j-1] ) )  

	out_filename_u = filename[-27:-35]+'_'+str(time[tind])+'_'+'U.nc'
	print 'Saving reformat to: '+out_filename_u

	netcdf_file = nc.Dataset(out_filename_u, 'w', format='NETCDF4')
	netcdf_file.createDimension('time', None)
	netcdf_file.createDimension('yu_c', len(yu_c))
	netcdf_file.createDimension('xu_c', len(xu_c))

	# variables
	ti = netcdf_file.createVariable('time', 'f4', ('time',))
	ti.units = 'days since 0001-01-01 00:00:00'
	ti.calender_type = 'JULIAN'
	ti.calender = 'JULIAN'

	yu_c_var = netcdf_file.createVariable('yu_c', 'f4', ('yu_c',))
	yu_c_var.units = 'degrees_N'
	yu_c_var.long_name = 'C grid U point latitude'

	xu_c_var = netcdf_file.createVariable('xu_c', 'f4', ('xu_c',))
	xu_c_var.units = 'degrees_E'
	xu_c_var.long_name = 'C grid U point longitude'

	u_var = netcdf_file.createVariable('u', 'f4', ('yu_c', 'xu_c',))
	u_var.units = 'm/s'
	u_var.long_name = 'zonal velocity on C grid'

	ti[:] = time[tind]
	yu_c = yu_c[:]
	xu_c = xu_c[:]
	u_var[...] = u_C 

	netcdf_file.close()

	out_filename_v = filename[-27:-35]+'_'+str(time[tind])+'_'+'V.nc'
	print 'Saving reformat to: '+out_filename_v

	netcdf_file = nc.Dataset(out_filename_v, 'w', format='NETCDF4')
	netcdf_file.createDimension('time', None)
	netcdf_file.createDimension('yv_c', len(yv_c))
	netcdf_file.createDimension('xv_c', len(xv_c))

	# variables
	ti = netcdf_file.createVariable('time', 'f4', ('time',))
	ti.units = 'days since 0001-01-01 00:00:00'
	ti.calender_type = 'JULIAN'
	ti.calender = 'JULIAN'

	yv_c_var = netcdf_file.createVariable('yv_c', 'f4', ('yv_c',))
	yv_c_var.units = 'degrees_N'
	yv_c_var.long_name = 'C grid V point latitude'

	xv_c_var = netcdf_file.createVariable('xv_c', 'f4', ('xv_c',))
	xv_c_var.units = 'degrees_E'
	xv_c_var.long_name = 'C grid V point longitude'

	v_var = netcdf_file.createVariable('u', 'f4', ('yv_c', 'xv_c',))
	v_var.units = 'm/s'
	v_var.long_name = 'meridional velocity on C grid'

	ti[:] = time[tind] 
	yv_c = yv_c[:]
	xv_c = xv_c[:]
	v_var[...] = v_C

	netcdf_file.close()