import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import sys
import os
import argparse


parser = argparse.ArgumentParser(
        description='generate daily surface C grid velocites files from a years worth of B grid daily velocities ')

parser.add_argument('--C_grid_spec', metavar='C_grid_spec', default='/archive/net/CM2.6/regrid/CM2p6_half_len.nc',
                    help='the netcdf file containing the C grid half-face lengths')

parser.add_argument('--B_vel', type=str, metavar='B_vel',
                    help='netcdf file containing a years worth of B grid daily surface velocities')

parser.add_argument('--out', type=str, metavar='out', 
                    help='the output folder for C grid daily surface velocity files')

args = parser.parse_args()

# storagepath = '/archive/net/CM2.6/regrid/'
# out_folder = '/archive/net/CM2.6/regrid/02020101/'

half_len_data = nc.Dataset(args.C_grid_spec)
lat_T = half_len_data.variables['grid_y_T'][0:2107]
lon_T = half_len_data.variables['grid_x_T'][...]
du_Eface_N = half_len_data.variables['ds_21_22_T'][0:2107, :] #omits polar region
du_Eface_S = half_len_data.variables['ds_20_21_T'][0:2107, :] #omits polar region 
du_Eface_N = np.transpose(du_Eface_N) 
du_Eface_S = np.transpose(du_Eface_S) 
half_len_data.close()


# path2file = '/archive/Richard.Slater/CM2.6/CM2.6_A_Control-1860_V03/history/'
# filename = '02020101.ocean_minibling_surf_field.nc'

vel_data = nc.Dataset(args.B_vel)
lat_U = vel_data.variables['yu_ocean'][0:2107]
lon_U = vel_data.variables['xu_ocean'][...]
vel_data.close()

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
vel_data = nc.Dataset(args.B_vel)
time = vel_data.variables['time'][:]

for tind in range(len(time)):
	u_B = vel_data.variables['usurf'][tind, 0:2107,:] #omits polar region
	v_B = vel_data.variables['vsurf'][tind, 0:2107,:] #omits polar region 

	u_B = np.transpose(u_B) # convert CM2.6 output from (lat, lon) indexing to (lon, lat) 
	v_B = np.transpose(v_B) # convert CM2.6 output from (lat, lon) indexing to (lon, lat)

	u_B = ma.filled(u_B, fill_value=0.0) #change fill value from -1e20 to 0.0 
	v_B = ma.filled(v_B, fill_value=0.0) #change fill value from -1e20 to 0.0 

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

	out_file_u = args.out+'U'+str(tind).zfill(10)+'.nc'
	print 'Saving reformat to: '+out_file_u

	netcdf_file = nc.Dataset(out_file_u, 'w', format='NETCDF4')
	netcdf_file.createDimension('time', 1)
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

	u_var = netcdf_file.createVariable('u', 'f4', ('xu_c', 'yu_c',))
	u_var.units = 'm/s'
	u_var.long_name = 'zonal velocity on C grid'

	ti[0] = time[tind]
	yu_c = yu_c[:]
	xu_c = xu_c[:]
	u_var[...] = u_C 

	netcdf_file.close()

	out_filename_v = args.out+'V'+str(tind).zfill(10)+'.nc'
	print 'Saving reformat to: '+out_filename_v

	netcdf_file = nc.Dataset(out_filename_v, 'w', format='NETCDF4')
	netcdf_file.createDimension('time', 1)
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

	v_var = netcdf_file.createVariable('u', 'f4', ('xv_c', 'yv_c',))
	v_var.units = 'm/s'
	v_var.long_name = 'meridional velocity on C grid'

	ti[0] = time[tind] 
	yv_c = yv_c[:]
	xv_c = xv_c[:]
	v_var[...] = v_C

	netcdf_file.close()
