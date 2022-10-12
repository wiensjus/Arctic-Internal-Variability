# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:38:41 2021

@author: Justin
"""

import netCDF4 as nc
import numpy as np
import numpy.matlib
import os
import pickle

# File with the changes in climate: (ts, temp) (TS,T,Q)
changefile = 'demodata/changefields.nc'

# File with initial surface SW downwelling and net radiative fields (for calculating
# albedo).
basefile = 'demodata/basefields.nc'

###############################################

# Read in coordinate info
data_dir = "/export/data/CESM-LENS/proc/cam5-kernels/"
data_dir2 = "/export/data/CESM-LENS/"
ds1 = nc.Dataset(data_dir + "/" + 'kernels/PS.nc')
ds2 = nc.Dataset(data_dir + "/" + 'kernels/t.kernel.nc') # Gaussian weights for the CESM grid

lat = ds1['lat'][:]
lon = ds1['lon'][:]
gw = ds2['gw'][:]
lev = ds2['lev'][:]

# Make an area weighting matrix
weight = np.tile(gw.T, [len(lon), 1]) 
weight = np.divide(weight, np.nansum(weight))
weight = weight[:,-32:]
weight_mon = np.tile(weight.T, [12, 1, 1]).T
print(np.shape(weight_mon))

#FLNT Data
filenames1b = os.listdir(data_dir2 + "/" + 'FLNT')

filenames1a = []
for i in range(len(filenames1b)):
    if '.nc' in filenames1b[i]:
        filenames1a.append(filenames1b[i])
        
#T Data
filenames2b = os.listdir(data_dir2 + "/" + 'T')

filenames2a = []
for i in range(len(filenames2b)):
    if '.nc' in filenames2b[i]:
        filenames2a.append(filenames2b[i])
        
data1 = []
for i in range(len(filenames1a)):
    data1.append(nc.Dataset(data_dir2 + "/" + 'FLNT' + "/" + filenames1a[i]))
    
dts_globalmean = []
        
for k in range(len(data1)):
    
    print(k + 1)
    
    infile = open(data_dir2 + '/TS/dts_monthly' + np.str(k + 1) + '_1955', 'rb')
    dts = np.array(pickle.load(infile))
    infile.close()
    dts = dts[:,-32:,:]
    
    dts_globalmean.append(np.nansum(np.nansum(np.multiply(np.nanmean(dts, axis=2), weight),axis=1),axis=0)/np.nansum(weight))
    print(dts_globalmean[k])

outfile = open('dts_arctic_annual_1955', 'wb')
pickle.dump(dts_globalmean, outfile)
outfile.close()

np.savetxt('dts_arctic_annual_1955.csv', dts_globalmean, delimiter = ',')

print('done')
