# -*- coding: utf-8 -*-
"""
Created on Fri May 28 13:04:24 2021

@author: Justin
"""

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import numpy.matlib
import os
import pickle

##### Albedo feedback

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
weight_mon = np.tile(weight.T, [12, 1, 1]).T

#T Data
filenames1b = os.listdir(data_dir2 + "/" + 'T')

filenames1a = []
for i in range(len(filenames1b)):
    if '.nc' in filenames1b[i]:
        filenames1a.append(filenames1b[i])
        
alb_feedback = []
alb_feedback_mon = []
alb_feedback_flux = []
alb_feedback_flux_mon = []
        
for k in range(len(filenames1a)):
    
    print(k + 1)
    
    infile = open(data_dir2 + '/TS/dts_monthly' + np.str(k + 1) + '_1955', 'rb')
    dts = np.array(pickle.load(infile))
    infile.close()
    
    infile = open(data_dir2 + "/FSNS/FSNS_Base_monthly" + np.str(k + 1) + '_1955', 'rb')
    SW_sfc_net_1 = np.array(pickle.load(infile))
    infile.close()
   
    infile = open(data_dir2 + "/FSNS/dFSNS_monthly" + np.str(k + 1) + '_1955', 'rb')
    SW_sfc_net_2 = np.array(pickle.load(infile)) + SW_sfc_net_1
    infile.close() 
    
    if k == 0:
   
        infile = open(data_dir2 + "/FSDS/FSDS_Base_monthly" + np.str(len(filenames1a)) + '_1955', 'rb')
        SW_sfc_down_1 = np.array(pickle.load(infile))
        infile.close()
        
        infile = open(data_dir2 + "/FSDS/dFSDS_monthly" + np.str(len(filenames1a)) + '_1955', 'rb')
        SW_sfc_down_2 = np.array(pickle.load(infile)) + SW_sfc_down_1
        infile.close()
        
    else:
        
        infile = open(data_dir2 + "/FSDS/FSDS_Base_monthly" + np.str(k) + '_1955', 'rb')
        SW_sfc_down_1 = np.array(pickle.load(infile))
        infile.close()
        
        infile = open(data_dir2 + "/FSDS/dFSDS_monthly" + np.str(k) + '_1955', 'rb')
        SW_sfc_down_2 = np.array(pickle.load(infile)) + SW_sfc_down_1
        infile.close()
   
    # Calculate the change in global mean surface temperature
    dts_globalmean = np.nansum(np.nansum(np.multiply(np.nanmean(dts,axis=2), weight),axis=1),axis=0)
    dts_globalmean_mon = np.nansum(np.nansum(np.multiply(dts, weight_mon),axis=1),axis=0)
    
    print(dts_globalmean)

    SW_sfc_net_1[SW_sfc_net_1 == 0] = np.nan
    SW_sfc_down_1[SW_sfc_down_1 == 0] = np.nan       
    alb1 = np.squeeze(1-SW_sfc_net_1/SW_sfc_down_1)
    alb1 = ma.filled(alb1, 0)
    alb1[alb1 == np.nan] = 0 
   
    #print(np.nansum(np.nansum(np.nansum(SW_sfc_net_1, 2), 1), 0))
    #print(np.nansum(np.nansum(np.nansum(SW_sfc_down_1, 2), 1), 0))

    SW_sfc_net_2[SW_sfc_net_2 == 0] = np.nan
    SW_sfc_down_2[SW_sfc_down_2 == 0] = np.nan 
    alb2 = np.squeeze(1-SW_sfc_net_2/SW_sfc_down_2)
    alb2 = ma.filled(alb2, 0)
    alb2[alb2 == np.nan] = 0 

    #print(np.nansum(np.nansum(np.nansum(SW_sfc_net_2, 2), 1), 0))
    #print(np.nansum(np.nansum(np.nansum(SW_sfc_down_2, 2), 1), 0))
    
    dalb = (alb2 - alb1) * 100
    
    # Read TOA albedo kernel
    ds9 = nc.Dataset(data_dir + "/" + 'kernels/alb.kernel.nc')
    alb_kernel = ds9['FSNT'][:].T
    
    dSW_alb = alb_kernel * dalb
    
    dSW_alb_globalmean = np.nansum(np.nansum(np.nanmean(dSW_alb, 2) * weight, 1), 0)
    dSW_alb_globalmean_mon = np.nansum(np.nansum((dSW_alb * weight_mon), 1), 0)

    alb_feedback.append(dSW_alb_globalmean/dts_globalmean)
    alb_feedback_mon.append(dSW_alb_globalmean_mon/dts_globalmean_mon)
    alb_feedback_flux.append(dSW_alb_globalmean)
    alb_feedback_flux_mon.append(dSW_alb_globalmean_mon)
    
    print(['Annual surface albedo feedback: ' + np.str(alb_feedback[k]) + ' W m^-2 K^-1'])
    print(['Monthly surface albedo feedback: ' + np.str(np.round(alb_feedback_mon[k], 3)) + ' W m^-2 K^-1'])
    
outfile = open('Albedo_Feedback_1955', 'wb')
pickle.dump(alb_feedback, outfile)
outfile.close()

outfile = open('Albedo_Feedback_Monthly_1955', 'wb')
pickle.dump(alb_feedback_mon, outfile)
outfile.close()

outfile = open('Albedo_Feedback_Flux_1955', 'wb')
pickle.dump(alb_feedback_flux, outfile)
outfile.close()

outfile = open('Albedo_Feedback_Flux_Monthly_1955', 'wb')
pickle.dump(alb_feedback_flux_mon, outfile)
outfile.close()

np.savetxt('Albedo_Feedback_1955.csv', alb_feedback, delimiter = ',')
np.savetxt('Albedo_Feedback_Monthly_1955.csv', alb_feedback_mon, delimiter = ',')

np.savetxt('Albedo_Feedback_Flux_1955.csv', alb_feedback_flux, delimiter = ',')
np.savetxt('Albedo_Feedback_Flux_Monthly_1955.csv', alb_feedback_flux_mon, delimiter = ',')
