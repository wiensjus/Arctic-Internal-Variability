# -*- coding: utf-8 -*-
"""
Created on Thu May 20 13:55:09 2021

@author: Justin
"""

##       Hybrid sigma-pressure vertical coordinate %%                                                                                                                                 
# This script will show you how to calculate 
# Top-Of-Atmosphere, clear-sky radiative feedbacks
# using the CESM-CAM5 radiative kernels. 
# In addition to the kernels and their accompanying 
# data, you'll need a set of T, q, and albedo changes
# on the CESM 0.9 degree grid for each month of the year. 
##############################################

import netCDF4 as nc
import numpy as np
import numpy.matlib
import math
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
weight = weight[:, -32:]
weight_mon = np.tile(weight.T, [12, 1, 1]).T

#T Data
filenames2b = os.listdir(data_dir2 + "/" + 'T')

filenames2a = []
for i in range(len(filenames2b)):
    if '.nc' in filenames2b[i]:
        filenames2a.append(filenames2b[i])
    
# Read midpoint pressure for each grid cell (lat,lon,level,month), [Pa]
    
ds6 = nc.Dataset(data_dir2 + "/" + "T" + "/" + filenames2a[0])
p = ds6['lev'][:] # [hPa]
p = np.tile(p.T, [12, 1]).T
p = np.tile(p, [len(lat), len(lon), 1, 1])
p = np.transpose(p, (1, 0, 2, 3))
            
# Read midpoint pressure for each grid cell (lat,lon,level,month), [Pa]

# Crude tropopause estimate: 100 hPa in the tropics, lowering with
# cosine to 300 hPa at the poles.
x = np.cos(lat*math.pi/180)
p_tropopause_zonalmean = np.subtract(300, np.multiply(200, x))
A = np.tile(p_tropopause_zonalmean.T, [len(lon), 1])
B = np.transpose(np.expand_dims(A, axis=2), (0,1,2))
C = np.tile(B, [1, 1, len(lev)])
p_tropopause=np.tile(np.transpose(np.expand_dims(C, axis=3), (0,1,2,3)), [1, 1, 1, 12])
p = p[:,-32:,:,:]
p_tropopause = p_tropopause[:,-32:,:,:]

# Read TOA Longwave surface temperature kernel
ds4 = nc.Dataset(data_dir + "/" + 'kernels/ts.kernel.nc')
ts_kernel = ds4['FLNT'][:].T
ts_kernel = ts_kernel[:,-32:,:]

ds7 = nc.Dataset(data_dir + "/" + 'kernels/t.kernel.nc')
ta_kernel = ds7['FLNT'][:].T
ta_kernel = ta_kernel[:,-32:,:,:]

t_feedback = []
planck_feedback = []
lapserate_feedback = []
planck_feedback_flux = []
lapserate_feedback_flux = []

t_feedback_mon = []
planck_feedback_mon = []
lapserate_feedback_mon = []
planck_feedback_flux_mon = []
lapserate_feedback_flux_mon = []

for k in range(len(filenames2a)):
    
    print(k + 1)
    
    # Read surface temperature change
    infile = open(data_dir2 + '/TS/dts_monthly' + np.str(k + 1) + '_1955', 'rb')
    dts = np.array(pickle.load(infile))
    infile.close()
    dts = dts[:,-32:,:]

    # Calculate the change in global mean surface temperature
    dts_globalmean = np.nansum(np.nansum(np.multiply(np.nanmean(dts, axis=2), weight), axis=1), axis=0)/np.nansum(weight)
    dts_globalmean_mon = np.nansum(np.nansum(np.multiply(dts, weight_mon),axis=1),axis=0)
    
    for i in range(len(dts_globalmean_mon)):
        dts_globalmean_mon[i] = dts_globalmean_mon[i]/np.nansum(weight)

    ### Temperature feedback calculation

    # Multiply monthly mean TS change by the TS kernels (function of lat, lon, month) (units W/m2)
    dLW_ts = ts_kernel * dts[k]

     # Read air temperature change [lon,lat,level,month]
    #ds5 = nc.Dataset(data_dir + "/" + changefile)
    #dta = ds5['temp'][:].T
    
    if k == 0:
    
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_1_1955', 'rb')
        dta1 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_2_1955', 'rb')
        dta2 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_3_1955', 'rb')
        dta3 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_4_1955', 'rb')
        dta4 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_5_1955', 'rb')
        dta5 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_6_1955', 'rb')
        dta6 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_7_1955', 'rb')
        dta7 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_8_1955', 'rb')
        dta8 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_9_1955', 'rb')
        dta9 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_10_1955', 'rb')
        dta10 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_11_1955', 'rb')
        dta11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(filenames2a)) + '_month_12_1955', 'rb')
        dta12 = np.array(pickle.load(infile)).T
        infile.close()
        
        dta = [dta1, dta2, dta3, dta4, dta5, dta6, dta7, dta8, dta9, dta10, dta11, dta12]
        dta = np.transpose(dta, (3, 2, 1, 0))
        dta = dta[:,-32:,:,:]
        
    else:
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_1_1955', 'rb')
        dta1 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_2_1955', 'rb')
        dta2 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_3_1955', 'rb')
        dta3 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_4_1955', 'rb')
        dta4 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_5_1955', 'rb')
        dta5 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_6_1955', 'rb')
        dta6 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_7_1955', 'rb')
        dta7 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_8_1955', 'rb')
        dta8 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_9_1955', 'rb')
        dta9 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_10_1955', 'rb')
        dta10 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_11_1955', 'rb')
        dta11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(k) + '_month_12_1955', 'rb')
        dta12 = np.array(pickle.load(infile)).T
        infile.close()
        
        dta = [dta1, dta2, dta3, dta4, dta5, dta6, dta7, dta8, dta9, dta10, dta11, dta12]
        dta = np.transpose(dta, (3, 2, 1, 0))
        dta = dta[:,-32:,:,:]
        
    # Set the temperature change to zero in the stratosphere (mask out stratosphere)
    dta = dta * (p >= p_tropopause)

    # Convolve air temperature kernel with air temperature change
    dLW_ta = np.squeeze(np.sum(ta_kernel * dta,2))

    # Add the surface and air temperature response; Take the annual average and global area average 
    dLW_t_globalmean = np.nansum(np.nansum(np.nanmean(-dLW_ta - dLW_ts,2)*weight,1),0)/np.nansum(weight)
    dLW_t_globalmean_mon = np.nansum(np.nansum(((-dLW_ta - dLW_ts) * weight_mon), 1), 0)

    for i in range(len(dLW_t_globalmean_mon)):
        dLW_t_globalmean_mon[i] = dLW_t_globalmean_mon[i]/np.nansum(weight)

    # Divide by the global annual mean surface warming (units: W/m2/K)
    t_feedback.append(dLW_t_globalmean/dts_globalmean)
    t_feedback_mon.append(dLW_t_globalmean_mon/dts_globalmean_mon)

    #print(['Temperature feedback: ' + np.str(t_feedback) + ' W m^-2 K^-1'])

    ##### Planck feedback: vertically uniform temperature change.                                                                                                                            
    # Project surface temperature change into height 
    dts3d = np.tile(np.transpose(np.expand_dims(dts, axis=3), (0,1,3,2)),[1, 1, 30, 1])
    
    # Mask stratosphere
    dt_planck = dts3d * (p>=p_tropopause)
    
    # Convolve air temperature kernel with 3-d surface air temp change
    dLW_planck = np.squeeze(np.sum(ta_kernel * dt_planck,2))
    
    # Take the annual average and global area average; incorporate the
    # part due to surface temperature change itself 
    dLW_planck_globalmean=np.nansum(np.nansum(np.nanmean(-dLW_planck - dLW_ts, 2)* weight, 1), 0)/np.nansum(weight)
    dLW_planck_globalmean_mon = np.nansum(np.nansum(((-dLW_planck - dLW_ts) * weight_mon), 1), 0)
    
    for i in range(len(dLW_planck_globalmean_mon)):
        dLW_planck_globalmean_mon[i] = dLW_planck_globalmean_mon[i]/np.nansum(weight)
    
    # Divide by the global annual mean surface warming (units: W/m2/K)
    planck_feedback.append(dLW_planck_globalmean/dts_globalmean)
    planck_feedback_mon.append(dLW_planck_globalmean_mon/dts_globalmean_mon)
    planck_feedback_flux.append(dLW_planck_globalmean)
    planck_feedback_flux_mon.append(dLW_planck_globalmean_mon)
     
    #### Lapse rate feedback                                                                                                                                                                 
    # Calculate the departure of temperature change from the surface
    # temperature change
    dt_lapserate = (dta - dt_planck) * (p >= p_tropopause)
    
    # Convolve air temperature kernel with 3-d surface air temp change
    dLW_lapserate = np.squeeze(np.sum(ta_kernel * dt_lapserate, 2))
    
    # Take the annual average and global area average 
    dLW_lapserate_globalmean = np.nansum(np.nansum(np.nanmean(-dLW_lapserate, 2) * weight, 1), 0)/np.nansum(weight)
    dLW_lapserate_globalmean_mon = np.nansum(np.nansum((-dLW_lapserate * weight_mon), 1), 0)
    
    for i in range(len(dLW_lapserate_globalmean_mon)):
        dLW_lapserate_globalmean_mon[i] = dLW_lapserate_globalmean_mon[i]/np.nansum(weight)
    
    # Divide by the global annual mean surface warming (units: W/m2/K)
    lapserate_feedback.append(dLW_lapserate_globalmean/dts_globalmean)
    lapserate_feedback_mon.append(dLW_lapserate_globalmean_mon/dts_globalmean_mon) 
    lapserate_feedback_flux.append(dLW_lapserate_globalmean)
    lapserate_feedback_flux_mon.append(dLW_lapserate_globalmean_mon)

    print('Annual arctic planck feedback: ' + np.str(planck_feedback[k]) + ' W m^-2 K^-1')
    print('Annual arctic lapse rate feedback: ' + np.str(lapserate_feedback[k]) + ' W m^-2 K^-1')
    
    print('Monthly arctic planck feedback: ' + np.str(np.round(planck_feedback_mon[k], 3)) + ' W m^-2 K^-1')
    print('Monthly arctic lapse rate feedback: ' + np.str(np.round(lapserate_feedback_mon[k], 3)) + ' W m^-2 K^-1')
    
outfile = open('Arctic_Planck_Feedback_1955', 'wb')
pickle.dump(planck_feedback, outfile)
outfile.close()

outfile = open('Arctic_Lapse_Rate_Feedback_1955', 'wb')
pickle.dump(lapserate_feedback, outfile)
outfile.close()

outfile = open('Arctic_Planck_Feedback_Monthly_1955', 'wb')
pickle.dump(planck_feedback_mon, outfile)
outfile.close()

outfile = open('Arctic_Lapse_Rate_Feedback_Monthly_1955', 'wb')
pickle.dump(lapserate_feedback_mon, outfile)
outfile.close()

outfile = open('Arctic_Planck_Feedback_Flux_1955', 'wb')
pickle.dump(planck_feedback_flux, outfile)
outfile.close()

outfile = open('Arctic_Lapse_Rate_Feedback_Flux_1955', 'wb')
pickle.dump(lapserate_feedback_flux, outfile)
outfile.close()

outfile = open('Arctic_Planck_Feedback_Flux_Monthly_1955', 'wb')
pickle.dump(planck_feedback_flux_mon, outfile)
outfile.close()

outfile = open('Arctic_Lapse_Rate_Feedback_Flux_Monthly_1955', 'wb')
pickle.dump(lapserate_feedback_flux_mon, outfile)
outfile.close()

np.savetxt('Arctic_Planck_Feedback_1955.csv', planck_feedback, delimiter = ',')
np.savetxt('Arctic_Lapse_Rate_Feedback_1955.csv', lapserate_feedback, delimiter = ',')

np.savetxt('Arctic_Planck_Feedback_Monthly_1955.csv', planck_feedback_mon, delimiter = ',')
np.savetxt('Arctic_Lapse_Rate_Feedback_Monthly_1955.csv', lapserate_feedback_mon, delimiter = ',')

np.savetxt('Arctic_Planck_Feedback_Flux_1955.csv', planck_feedback_flux, delimiter = ',')
np.savetxt('Arctic_Lapse_Rate_Feedback_Flux_1955.csv', lapserate_feedback_flux, delimiter = ',')

np.savetxt('Arctic_Planck_Feedback_Flux_Monthly_1955.csv', planck_feedback_flux_mon, delimiter = ',')
np.savetxt('Arctic_Lapse_Rate_Feedback_Flux_Monthly_1955.csv', lapserate_feedback_flux_mon, delimiter = ',')

print('done')
