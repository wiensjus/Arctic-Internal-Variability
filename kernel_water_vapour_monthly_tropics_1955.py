# -*- coding: utf-8 -*-
"""
Created on Thu May 20 19:06:43 2021

@author: Justin
"""

####### CESM-CAM5 Radiative Kernel Demo: %%%%%%%%%%%%
##       Hybrid sigma-pressure vertical coordinate %%
# This script will show you how to calculate 
# Top-Of-Atmosphere, clear-sky radiative feedbacks
# using the CESM-CAM5 radiative kernels. 
# In addition to the kernels and their accompanying 
# data, you'll need a set of T, q, and albedo changes
# on the CESM 0.9 degree grid for each month of the year. 
##################################################

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
weight = weight[:,63:129]
weight_mon = np.tile(weight.T, [12, 1, 1]).T

#T Data
filenames2b = os.listdir(data_dir2 + "/" + 'T')

filenames2a = []
for i in range(len(filenames2b)):
    if '.nc' in filenames2b[i]:
        filenames2a.append(filenames2b[i])
            
# Read midpoint pressure for each grid cell (lat,lon,level,month), [Pa]
    
ds6 = nc.Dataset(data_dir2 + "/" + "T" + "/" + filenames2b[0])
p = ds6['lev'][:] # [hPa]
p = np.tile(p.T, [12, 1]).T
p = np.tile(p, [len(lat), len(lon), 1, 1])
p = np.transpose(p, (1, 0, 2, 3))

# Crude tropopause estimate: 100 hPa in the tropics, lowering with
# cosine to 300 hPa at the poles.
x = np.cos(lat*math.pi/180)
p_tropopause_zonalmean = np.subtract(300, np.multiply(200, x))
A = np.tile(p_tropopause_zonalmean.T, [len(lon), 1])
B = np.transpose(np.expand_dims(A, axis=2), (0,1,2))
C = np.tile(B, [1, 1, len(lev)])
p_tropopause=np.tile(np.transpose(np.expand_dims(C, axis=3), (0,1,2,3)), [1, 1, 1, 12])
p = p[:,63:129,:,:]
p_tropopause = p_tropopause[:,63:129,:,:]

def calcsatspechum(t, p):
    
    # T is temperature, P is pressure in hPa 
    # Formulae from Buck (1981)
    
    es = (1.0007 + (3.46e-6 * p)) * 6.1121 * np.exp(17.502 * (t - 273.15) / (240.97 + (t - 273.15)))
    wsl = 0.622 * es/(p - es) # saturation mixing ratio wrt liquid water (g/kg)   
    
    es = (1.0003 + (4.18e-6 * p)) * 6.1115 * np.exp(22.452 * (t - 273.15) / (272.55 + (t - 273.15)))
    
    wsi = 0.622 * es/(p - es) # saturation mixing ratio wrt ice (g/kg)
    
    ws = wsl
    ws[t < 273.15] = wsi[t < 273.15]
    
    qs = ws/(1 + ws) # saturation specific humidity, g/kg
    return qs

####### Water vapor feedback : using logarithmic q 

# Calculate the change in moisture per degree warming at constant relative humidity. 
# Run the accompanying NCL script with your input files, or
# implement here. 

# Read kernels
ds5 = nc.Dataset(data_dir + "/" + 'kernels/q.kernel.nc')
q_LW_kernel = ds5['FLNT'][:].T
q_SW_kernel = ds5['FSNT'][:].T
q_LW_kernel = q_LW_kernel[:,63:129,:,:]
q_SW_kernel = q_SW_kernel[:,63:129,:,:]

q_feedback = []
q_feedback_mon = []
q_feedback_flux = []
q_feedback_flux_mon = []

for k in range(len(filenames2a)):
    
    print(k + 1)
    
    infile = open(data_dir2 + '/TS/dts_monthly' + np.str(k + 1) + '_1955', 'rb')
    dts = np.array(pickle.load(infile))
    infile.close()
    dts = dts[:,63:129,:]

    # Calculate the change in global mean surface temperature
    dts_globalmean = np.nansum(np.nansum(np.multiply(np.nanmean(dts, axis=2), weight), axis=1), axis=0)/np.nansum(weight)
    dts_globalmean_mon = np.nansum(np.nansum(np.multiply(dts, weight_mon),axis=1),axis=0)
    
    for i in range(len(dts_globalmean_mon)):
        dts_globalmean_mon[i] = dts_globalmean_mon[i]/np.nansum(weight)

    # Read data
    
    if k == 0:
        
        #Q_Base data
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_1_1955', 'rb')
        q01 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_2_1955', 'rb')
        q02 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_3_1955', 'rb')
        q03 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_4_1955', 'rb')
        q04 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_5_1955', 'rb')
        q05 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_6_1955', 'rb')
        q06 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_7_1955', 'rb')
        q07 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_8_1955', 'rb')
        q08 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_9_1955', 'rb')
        q09 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_10_1955', 'rb')
        q010 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_11_1955', 'rb')
        q011 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(len(filenames2a)) + '_month_12_1955', 'rb')
        q012 = np.array(pickle.load(infile)).T
        infile.close()
            
        #dQ data
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_1_1955', 'rb')
        dq1 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_2_1955', 'rb')
        dq2 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_3_1955', 'rb')
        dq3 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_4_1955', 'rb')
        dq4 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_5_1955', 'rb')
        dq5 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_6_1955', 'rb')
        dq6 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_7_1955', 'rb')
        dq7 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_8_1955', 'rb')
        dq8 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_9_1955', 'rb')
        dq9 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_10_1955', 'rb')
        dq10 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_11_1955', 'rb')
        dq11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(len(filenames2a)) + '_month_12_1955', 'rb')
        dq12 = np.array(pickle.load(infile)).T
        infile.close()
                
        #Tair_Base data
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_1_1955', 'rb')
        t11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_2_1955', 'rb')
        t12 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_3_1955', 'rb')
        t13 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_4_1955', 'rb')
        t14 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_5_1955', 'rb')
        t15 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_6_1955', 'rb')
        t16 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_7_1955', 'rb')
        t17 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_8_1955', 'rb')
        t18 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_9_1955', 'rb')
        t19 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_10_1955', 'rb')
        t110 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_11_1955', 'rb')
        t111 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(len(filenames2a)) + '_month_12_1955', 'rb')
        t112 = np.array(pickle.load(infile)).T
        infile.close()
            
        #dta data
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
                
        q0 = [q01, q02, q03, q04, q05, q06, q07, q08, q09, q010, q011, q012]
        q0 = np.transpose(q0, (1, 2, 3, 0))
        
        dq = [dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12]
        dq = np.transpose(dq, (1, 2, 3, 0))
        
        t1 = [t11, t12, t13, t14, t15, t16, t17, t18, t19, t110, t111, t112]
        t1 = np.transpose(t1, (1, 2, 3, 0))
        
        dta = [dta1, dta2, dta3, dta4, dta5, dta6, dta7, dta8, dta9, dta10, dta11, dta12]
        dta = np.transpose(dta, (3, 2, 1, 0))
        
        q0 = q0[:,63:129,:,:]
        dq = dq[:,63:129,:,:]
        t1 = t1[:,63:129,:,:]
        dta = dta[:,63:129,:,:]
        
    else:
        
        #Q_Base data
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_1_1955', 'rb')
        q01 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_2_1955', 'rb')
        q02 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_3_1955', 'rb')
        q03 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_4_1955', 'rb')
        q04 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_5_1955', 'rb')
        q05 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_6_1955', 'rb')
        q06 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_7_1955', 'rb')
        q07 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_8_1955', 'rb')
        q08 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_9_1955', 'rb')
        q09 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_10_1955', 'rb')
        q010 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_11_1955', 'rb')
        q011 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'Q_Base_monthly' + np.str(k) + '_month_12_1955', 'rb')
        q012 = np.array(pickle.load(infile)).T
        infile.close()
    
        #dQ data
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_1_1955', 'rb')
        dq1 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_2_1955', 'rb')
        dq2 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_3_1955', 'rb')
        dq3 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_4_1955', 'rb')
        dq4 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_5_1955', 'rb')
        dq5 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_6_1955', 'rb')
        dq6 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_7_1955', 'rb')
        dq7 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_8_1955', 'rb')
        dq8 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_9_1955', 'rb')
        dq9 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_10_1955', 'rb')
        dq10 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_11_1955', 'rb')
        dq11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "Q" + "/" + 'dQ_monthly' + np.str(k) + '_month_12_1955', 'rb')
        dq12 = np.array(pickle.load(infile)).T
        infile.close()
        
        #Tair_Base data
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_1_1955', 'rb')
        t11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_2_1955', 'rb')
        t12 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_3_1955', 'rb')
        t13 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_4_1955', 'rb')
        t14 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_5_1955', 'rb')
        t15 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_6_1955', 'rb')
        t16 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_7_1955', 'rb')
        t17 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_8_1955', 'rb')
        t18 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_9_1955', 'rb')
        t19 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_10_1955', 'rb')
        t110 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_11_1955', 'rb')
        t111 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'Tair_Base_monthly' + np.str(k) + '_month_12_1955', 'rb')
        t112 = np.array(pickle.load(infile)).T
        infile.close()
    
        #dta data
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
        
        q0 = [q01, q02, q03, q04, q05, q06, q07, q08, q09, q010, q011, q012]
        q0 = np.transpose(q0, (1, 2, 3, 0))
        
        dq = [dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12]
        dq = np.transpose(dq, (1, 2, 3, 0))
        
        t1 = [t11, t12, t13, t14, t15, t16, t17, t18, t19, t110, t111, t112]
        t1 = np.transpose(t1, (1, 2, 3, 0))
        
        dta = [dta1, dta2, dta3, dta4, dta5, dta6, dta7, dta8, dta9, dta10, dta11, dta12]
        dta = np.transpose(dta, (3, 2, 1, 0))
        
        q0 = q0[:,63:129,:,:]
        dq = dq[:,63:129,:,:]
        t1 = t1[:,63:129,:,:]
        dta = dta[:,63:129,:,:]
    
    qs1 = calcsatspechum(t1, p) # g/kg
    qs2 = calcsatspechum(t1 + dta, p) # g/kg
    
    dta[dta == 0] = np.nan
    dqsdt = (qs2 - qs1)/dta
 
    rh = 1000 * q0/qs1
    dqdt = rh * dqsdt

    dlogqdt = dqdt/(1000 * q0)

    # Normalize kernels by the change in moisture for 1 K warming at
    # constant RH (log-q kernel)
    logq_LW_kernel = q_LW_kernel/dlogqdt
    logq_SW_kernel = q_SW_kernel/dlogqdt

    # Mask out the stratosphere
    dq = dq * (p >= p_tropopause)

    dlogq = dq/q0

    # Convolve moisture kernel with change in moisture
    dLW_logq = np.squeeze(np.nansum(logq_LW_kernel * dlogq, 2))
    dSW_logq = np.squeeze(np.nansum(logq_SW_kernel * dlogq, 2))

    # Add the LW and SW responses. Note the sign convention difference
    # between LW and SW!
    dR_logq_globalmean = np.nansum(np.nansum(np.nanmean(-dLW_logq + dSW_logq, 2) * weight, 1), 0)/np.nansum(weight)
    dR_logq_globalmean_mon = np.nansum(np.nansum(((-dLW_logq + dSW_logq) * weight_mon), 1), 0)
    
    for i in range(len(dR_logq_globalmean_mon)):
        dR_logq_globalmean_mon[i] = dR_logq_globalmean_mon[i]/np.nansum(weight)

    # Divide by the global annual mean surface warming (units: W/m2/K)
    q_feedback.append(dR_logq_globalmean/dts_globalmean)
    q_feedback_mon.append(dR_logq_globalmean_mon/dts_globalmean_mon)
    q_feedback_flux.append(dR_logq_globalmean)
    q_feedback_flux_mon.append(dR_logq_globalmean_mon)

    print('Annual Tropics water vapor feedback: ' + np.str(q_feedback[k]) + ' W m^-2 K^-1')
    print('Monthly Tropics water vapor feedback: ' + np.str(np.round(q_feedback_mon[k], 3)) + ' W m^-2 K^-1')

outfile = open('Tropics_Water_Vapour_Feedback_1955', 'wb')
pickle.dump(q_feedback, outfile)
outfile.close()

outfile = open('Tropics_Water_Vapour_Feedback_Monthly_1955', 'wb')
pickle.dump(q_feedback_mon, outfile)
outfile.close()

outfile = open('Tropics_Water_Vapour_Feedback_Flux_1955', 'wb')
pickle.dump(q_feedback_flux, outfile)
outfile.close()

outfile = open('Tropics_Water_Vapour_Feedback_Flux_Monthly_1955', 'wb')
pickle.dump(q_feedback_flux_mon, outfile)
outfile.close()

np.savetxt('Tropics_Water_Vapour_Feedback_1955.csv', q_feedback, delimiter = ',')
np.savetxt('Tropics_Water_Vapour_Feedback_Monthly_1955.csv', q_feedback_mon, delimiter = ',')

np.savetxt('Tropics_Water_Vapour_Feedback_Flux_1955.csv', q_feedback_flux, delimiter = ',')
np.savetxt('Tropics_Water_Vapour_Feedback_Flux_Monthly_1955.csv', q_feedback_flux_mon, delimiter = ',')

print('done')

