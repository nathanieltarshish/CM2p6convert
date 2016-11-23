#!/usr/bin/env python
import numpy as np
import sys
import os
import netCDF4 as nc
import bisect
import argparse
import math

def write_field(fname, data):
    """writes out data to a .bin file, based on MITgcm verification tests"""
    print 'wrote to file: '+fname
    if sys.byteorder == 'little':
        data.byteswap(True)
    fid = open(fname, "wb")
    data.tofile(fid)
    fid.close()

parser = argparse.ArgumentParser(
        description='convert folder of daily surface C grid velocites to folder of .data  ')

parser.add_argument('--input', metavar='input', 
                    help='the input folder of netcdf files ')

parser.add_argument('--out', type=str, metavar='out', 
                    help='the output folder of .data of files')

args = parser.parse_args()


files = os.listdir(args.input)
files = np.sort(files)

for file in files:
	if not file.endswith('nc'):
		continue

	netcdf_file = nc.Dataset(args.input+file)
	outfile = args.out+file[0]+'.'+file[1:len(file)-3]+'.data'

	if file[0] == 'U':
		u = netcdf_file.variables['u'][...]
		write_field(outfile, u)
	
	if file[0] == 'V':
		v = netcdf_file.variables['u'][...]
		write_field(outfile, v)

