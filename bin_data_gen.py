#!/usr/bin/env python
import numpy as np
import sys
import os
import netCDF4 as nc
import bisect
import argparse
import math

"""Generates bathymetry mask and delY bin files from CM2p6 grid specification netcdf file.
   Notes: The bathymetry mask input netcdf file is assumed to have the format: 1 = Ocean, 0 = Land.
          The program uses the Tcell height to calculate delY and assumes the grid is perfectly spherical
          (R_EARTH is constant so that delY = Tcell_height/R_Earth). To graphically verify output, use
          the --plot option.
"""


def write_field(fname, data):
    """writes out data to a .bin file based on MITgcm verification tests"""
    print 'wrote to file: '+fname
    if sys.byteorder == 'little':
        data.byteswap(True)
    fid = open(fname, "wb")
    data.tofile(fid)
    fid.close()


def nearest_index(array, value):
    """returns index in array of nearest value in O(log n)"""
    y = bisect.bisect(array, value)
    if (y == len(array)) or (abs(array[y-1] - value) < abs(array[y] - value)):
        return y-1
    else:
        return y


parser = argparse.ArgumentParser(
        description='Generate 2D delY.bin and bathy.bin from CM2.6 grid spec file')

parser.add_argument('--input_nc', metavar='input_nc',
                    help='the netcdf file containing lat, lon, delY, and land mask')

parser.add_argument('--lat', type=str, metavar='lat_key', default='grid_y_T',
                    help='the name of geographic latitude variable in input_nc')

parser.add_argument('--lon', type=str, metavar='lon_key', default='grid_x_T',
                    help='the name of geographic longitude variable in input_nc')

parser.add_argument('--landmask', type=str, metavar='mask_key', default='wet',
                    help='the name of land mask variable in netcdf file, assumes land = 0, water = 1')

parser.add_argument('--Tcell_height', type=str, metavar='Tcell_height', default='ds_10_12_T',
                    help='the name of Tcell height variable in input_nc')

parser.add_argument('--out_dir', default='./', metavar='out_dir',
                    help='the folder for the .bin file ouputs')

parser.add_argument('--plot',
                    help='saves a plot to aid in verifying bathymetry and delY',
                    action='store_true')

parser.add_argument('--maxlat', default=65.0, metavar='maxlat',
                    help='the upper latitude limit of output grid (all data above is ignored)')

args = parser.parse_args()

ncpath = args.input_nc  # path to file containing grid data
if ncpath is not None:
    ncdata = nc.Dataset(ncpath)
else:
    raise ValueError('no input netcdf file specified (input_nc = None)')

out_dir = args.out_dir
Tcell_height_key = args.Tcell_height
mask_key = args.landmask
lat_key = args.lat
lon_key = args.lon
cutoff_value = args.maxlat

print "Using the following variable assignments: "
print "Tcell_height_key = ", Tcell_height_key
print "mask_key = ", mask_key
print "lat_key = ", lat_key
print "lon_key = ", lon_key

lat = ncdata.variables[lat_key][...]
lon = ncdata.variables[lon_key][...]

lat_max = np.shape(lat)[0]
cutoff_index = nearest_index(lat[:], cutoff_value)-1  # for CM2p6, index =  2107 value = 64.9732
# polar region to throw away exists between latitude interval [cutoff_index, lat_max]
lat = np.delete(lat, range(cutoff_index, lat_max))  # removes polar region

bathmask = ncdata.variables[mask_key][0:cutoff_index, :]
Z0 = -100.0  # arbitrary depth of ocean layer
R_EARTH = 6.371*10**6   # equatorial radius of the earth in meters
depth = (Z0*bathmask).astype(dtype='float32')
depth = depth.transpose()  # MITgcm format is indexed by (lon, lat)

Tcell_height = ncdata.variables[Tcell_height_key][0:cutoff_index,0]
delY = Tcell_height/R_EARTH*(180/math.pi) # units are degrees
delY = delY.astype(dtype='float32') 

if np.shape(delY)[0] != np.shape(bathmask)[0]:
    raise ValueError('Tcell height lat points different than bathymetry mask lat points')

if not os.path.exists(out_dir):  # create folder for outputs
    os.makedirs(out_dir)

if args.plot:
    try:
        import matplotlib.pyplot as plt
    except:
        raise ValueError('cannot import matplotlib')
    # plot depth and delY fields
    plt.subplot(2, 1, 1)
    plt.title('Depth [m]')
    plt.pcolormesh(lon, lat, depth.transpose(), cmap='bwr')
    plt.colorbar()
    plt.subplot(2, 1, 2)
    plt.title('delY [degrees]')
    plt.scatter(lat, delY,  marker='+')
    plt.xlabel("Latitude (deg)")
    plt.ylabel("delY [degrees]")
    plt.savefig(out_dir+"grid_verification.png")
    print "plot saved to", out_dir+"grid_verification.png"

# write out depth as bathy.bin
bath_file = out_dir+'bathy.bin'
write_field(bath_file, depth)

# write out delY as delY.bin
delY_file = out_dir+'delY.bin'
write_field(delY_file, delY)
