# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:22:21 2021

@author: Justin
"""
##    Cloud feedback. 
# This script will show you how to calculate 
# the Cloud Feedback at the Top-Of-Atmosphere 
# using the CESM-CAM5 radiative kernels. 
#############################################

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import numpy.matlib
import math
import os
import pickle
from scipy.interpolate import RectBivariateSpline

# File with the changes in climate: (ts, temp) (TS,T,Q)
changefile = 'demodata/changefields.nc'

# File with initial surface SW downwelling and net radiative fields (for calculating
# albedo).
basefile = 'demodata/basefields.nc'

###############################################

## STEP 1. Calculate total-sky and clear-sky feedbacks

# Read in coordinate info
data_dir = "/export/data/CESM-LENS/proc/cam5-kernels/"
data_dir2 = "/export/data/CESM-LENS/"
ds1 = nc.Dataset(data_dir + "/" + 'kernels/PS.nc')
ds2 = nc.Dataset(data_dir + "/" + 'kernels/t.kernel.nc') # Gaussian weights for the CESM grid

lat = ds1['lat'][:]
lon = ds1['lon'][:]
gw = ds2['gw'][:]
lev = ds2['lev'][:]

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

#print(filenames2a)
        
data1 = []
for i in range(len(filenames1a)):
    data1.append(nc.Dataset(data_dir2 + "/" + 'FLNT' + "/" + filenames1a[i]))

# Read midpoint pressure for each grid cell (lat,lon,level,month), [Pa]
ds6 = nc.Dataset(data_dir2 + "/" + "T" + "/" + filenames2a[0])
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

# Make an area weighting matrix
weight = np.tile(gw.T, [len(lon), 1]) 
weight = np.divide(weight, np.nansum(weight))
weight = weight[:,63:129]
weight_mon = np.tile(weight.T, [12, 1, 1]).T

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

# Read TOA Longwave surface temperature kernel
ds5 = nc.Dataset(data_dir + "/" + 'kernels/ts.kernel.nc')
ts_kernel = ds5['FLNT'][:].T
ts_kernel_clearsky = ds5['FLNTC'][:].T
ts_kernel = ts_kernel[:,63:129,:]
ts_kernel_clearsky = ts_kernel_clearsky[:,63:129,:]

# Read air temperature kernel
ds6 = nc.Dataset(data_dir + "/" + 'kernels/t.kernel.nc')
ta_kernel = ds6['FLNT'][:].T
ta_kernel_clearsky = ds6['FLNTC'][:].T
ta_kernel = ta_kernel[:,63:129,:,:]
ta_kernel_clearsky = ta_kernel_clearsky[:,63:129,:,:]

# Read TOA albedo kernel
ds8 = nc.Dataset(data_dir + "/" + 'kernels/alb.kernel.nc')
alb_kernel = ds8['FSNT'][:].T
alb_kernel_clearsky = ds8['FSNTC'][:].T
alb_kernel = alb_kernel[:,63:129,:]
alb_kernel_clearsky = alb_kernel_clearsky[:,63:129,:]

# Read q kernels
ds9 = nc.Dataset(data_dir + "/" + 'kernels/q.kernel.nc')
q_LW_kernel = ds9['FLNT'][:].T
q_SW_kernel = ds9['FSNT'][:].T
q_LW_kernel_clearsky = ds9['FLNTC'][:].T
q_SW_kernel_clearsky = ds9['FSNTC'][:].T
q_LW_kernel = q_LW_kernel[:,63:129,:,:]
q_SW_kernel = q_SW_kernel[:,63:129,:,:]
q_LW_kernel_clearsky = q_LW_kernel_clearsky[:,63:129,:,:]
q_SW_kernel_clearsky = q_SW_kernel_clearsky[:,63:129,:,:]

ghgfile1 = 'PORT_RF_1955_2005/cam.allforcing.01.nc'
ghgfile2 = 'PORT_RF_1955_2005/cam.allforcing.02.nc'
ghgfile3 = 'PORT_RF_1955_2005/cam.allforcing.03.nc'
ghgfile4 = 'PORT_RF_1955_2005/cam.allforcing.04.nc'
ghgfile5 = 'PORT_RF_1955_2005/cam.allforcing.05.nc'
ghgfile6 = 'PORT_RF_1955_2005/cam.allforcing.06.nc'
ghgfile7 = 'PORT_RF_1955_2005/cam.allforcing.07.nc'
ghgfile8 = 'PORT_RF_1955_2005/cam.allforcing.08.nc'
ghgfile9 = 'PORT_RF_1955_2005/cam.allforcing.09.nc'
ghgfile10 = 'PORT_RF_1955_2005/cam.allforcing.10.nc'
ghgfile11 = 'PORT_RF_1955_2005/cam.allforcing.11.nc'
ghgfile12 = 'PORT_RF_1955_2005/cam.allforcing.12.nc'

ds10 = nc.Dataset(data_dir2 + "/" + ghgfile1)
ds11 = nc.Dataset(data_dir2 + "/" + ghgfile2)
ds12 = nc.Dataset(data_dir2 + "/" + ghgfile3)
ds13 = nc.Dataset(data_dir2 + "/" + ghgfile4)
ds14 = nc.Dataset(data_dir2 + "/" + ghgfile5)
ds15 = nc.Dataset(data_dir2 + "/" + ghgfile6)
ds16 = nc.Dataset(data_dir2 + "/" + ghgfile7)
ds17 = nc.Dataset(data_dir2 + "/" + ghgfile8)
ds18 = nc.Dataset(data_dir2 + "/" + ghgfile9)
ds19 = nc.Dataset(data_dir2 + "/" + ghgfile10)
ds20 = nc.Dataset(data_dir2 + "/" + ghgfile11)
ds21 = nc.Dataset(data_dir2 + "/" + ghgfile12)

lat_ghg = ds10['lat'][:]
lon_ghg = ds10['lon'][:]
    
sw1 = ds10['FSNT'][:].T
sw_cs1 = ds10['FSNTC'][:].T
lw1 = ds10['FLNT'][:].T
lw_cs1 = ds10['FLNTC'][:].T

sw2 = ds11['FSNT'][:].T
sw_cs2 = ds11['FSNTC'][:].T
lw2 = ds11['FLNT'][:].T
lw_cs2 = ds11['FLNTC'][:].T

sw3 = ds12['FSNT'][:].T
sw_cs3 = ds12['FSNTC'][:].T
lw3 = ds12['FLNT'][:].T
lw_cs3 = ds12['FLNTC'][:].T

sw4 = ds13['FSNT'][:].T
sw_cs4 = ds13['FSNTC'][:].T
lw4 = ds13['FLNT'][:].T
lw_cs4 = ds13['FLNTC'][:].T

sw5 = ds14['FSNT'][:].T
sw_cs5 = ds14['FSNTC'][:].T
lw5 = ds14['FLNT'][:].T
lw_cs5 = ds14['FLNTC'][:].T

sw6 = ds15['FSNT'][:].T
sw_cs6 = ds15['FSNTC'][:].T
lw6 = ds15['FLNT'][:].T
lw_cs6 = ds15['FLNTC'][:].T

sw7 = ds16['FSNT'][:].T
sw_cs7 = ds16['FSNTC'][:].T
lw7 = ds16['FLNT'][:].T
lw_cs7 = ds16['FLNTC'][:].T

sw8 = ds17['FSNT'][:].T
sw_cs8 = ds17['FSNTC'][:].T
lw8 = ds17['FLNT'][:].T
lw_cs8 = ds17['FLNTC'][:].T

sw9 = ds18['FSNT'][:].T
sw_cs9 = ds18['FSNTC'][:].T
lw9 = ds18['FLNT'][:].T
lw_cs9 = ds18['FLNTC'][:].T

sw10 = ds19['FSNT'][:].T
sw_cs10 = ds19['FSNTC'][:].T
lw10 = ds19['FLNT'][:].T
lw_cs10 = ds19['FLNTC'][:].T

sw11 = ds20['FSNT'][:].T
sw_cs11 = ds20['FSNTC'][:].T
lw11 = ds20['FLNT'][:].T
lw_cs11 = ds20['FLNTC'][:].T

sw12 = ds21['FSNT'][:].T
sw_cs12 = ds21['FSNTC'][:].T
lw12 = ds21['FLNT'][:].T
lw_cs12 = ds21['FLNTC'][:].T


cloud_masking_of_forcing_sw1 = sw_cs1 - sw1
cloud_masking_of_forcing_lw1 = lw_cs1 - lw1

cloud_masking_of_forcing_sw2 = sw_cs2 - sw2
cloud_masking_of_forcing_lw2 = lw_cs2 - lw2

cloud_masking_of_forcing_sw3 = sw_cs3 - sw3
cloud_masking_of_forcing_lw3 = lw_cs3 - lw3

cloud_masking_of_forcing_sw4 = sw_cs4 - sw4
cloud_masking_of_forcing_lw4 = lw_cs4 - lw4

cloud_masking_of_forcing_sw5 = sw_cs5 - sw5
cloud_masking_of_forcing_lw5 = lw_cs5 - lw5

cloud_masking_of_forcing_sw6 = sw_cs6 - sw6
cloud_masking_of_forcing_lw6 = lw_cs6 - lw6

cloud_masking_of_forcing_sw7 = sw_cs7 - sw7
cloud_masking_of_forcing_lw7 = lw_cs7 - lw7

cloud_masking_of_forcing_sw8 = sw_cs8 - sw8
cloud_masking_of_forcing_lw8 = lw_cs8 - lw8

cloud_masking_of_forcing_sw9 = sw_cs9 - sw9
cloud_masking_of_forcing_lw9 = lw_cs9 - lw9

cloud_masking_of_forcing_sw10 = sw_cs10 - sw10
cloud_masking_of_forcing_lw10 = lw_cs10 - lw10

cloud_masking_of_forcing_sw11 = sw_cs11 - sw11
cloud_masking_of_forcing_lw11 = lw_cs11 - lw11

cloud_masking_of_forcing_sw12 = sw_cs12 - sw12
cloud_masking_of_forcing_lw12 = lw_cs12 - lw12
    
#print(np.shape(cloud_masking_of_forcing_sw2))
#print(len(lat))
#print(len(lon))
 
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw1.T)
cloud_sw1 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw1.T)
cloud_lw1 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw2.T)
cloud_sw2 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw2.T)
cloud_lw2 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw3.T)
cloud_sw3 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw3.T)
cloud_lw3 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw4.T)
cloud_sw4 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw4.T)
cloud_lw4 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw5.T)
cloud_sw5 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw5.T)
cloud_lw5 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw6.T)
cloud_sw6 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw6.T)
cloud_lw6 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw7.T)
cloud_sw7 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw7.T)
cloud_lw7 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw8.T)
cloud_sw8 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw8.T)
cloud_lw8 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw9.T)
cloud_sw9 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw9.T)
cloud_lw9 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw10.T)
cloud_sw10 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw10.T)
cloud_lw10 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw11.T)
cloud_sw11 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw11.T)
cloud_lw11 = ip(lat, lon)

ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_sw12.T)
cloud_sw12 = ip(lat, lon)   
ip = RectBivariateSpline(lat_ghg, lon_ghg, cloud_masking_of_forcing_lw12.T)
cloud_lw12 = ip(lat, lon)

cloud_masking_of_forcing_sw = np.array([cloud_sw1, cloud_sw2, cloud_sw3, cloud_sw4, cloud_sw5, cloud_sw6, cloud_sw7, cloud_sw8, cloud_sw9, cloud_sw10, cloud_sw11, cloud_sw12]).T
cloud_masking_of_forcing_lw = np.array([cloud_lw1, cloud_lw2, cloud_lw3, cloud_lw4, cloud_lw5, cloud_lw6, cloud_lw7, cloud_lw8, cloud_lw9, cloud_lw10, cloud_lw11, cloud_lw12]).T
cloud_masking_of_forcing_sw = cloud_masking_of_forcing_sw[:,63:129,:]
cloud_masking_of_forcing_lw = cloud_masking_of_forcing_lw[:,63:129,:]

LWc_feedback = []
LWc_feedback_mon = []
LWc_feedback_flux = []
LWc_feedback_flux_mon = []

SWc_feedback = []
SWc_feedback_mon = []
SWc_feedback_flux = []
SWc_feedback_flux_mon = []

for k in range(len(data1)):
    
    print(k + 1)

    ### Temperature feedback 
    
    if k == 0:
        
        infile = open(data_dir2 + '/TS/dts_monthly' + np.str(k + 1), 'rb')
        dts = np.array(pickle.load(infile))
        infile.close()
        dts = dts[:,63:129,:]
    
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_1_1955', 'rb')
        dta1 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_2_1955', 'rb')
        dta2 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_3_1955', 'rb')
        dta3 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_4_1955', 'rb')
        dta4 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_5_1955', 'rb')
        dta5 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_6_1955', 'rb')
        dta6 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_7_1955', 'rb')
        dta7 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_8_1955', 'rb')
        dta8 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_9_1955', 'rb')
        dta9 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_10_1955', 'rb')
        dta10 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_11_1955', 'rb')
        dta11 = np.array(pickle.load(infile)).T
        infile.close()
        
        infile = open(data_dir2 + "/" + "T" + "/" + 'dta_monthly' + np.str(len(data1)) + '_month_12_1955', 'rb')
        dta12 = np.array(pickle.load(infile)).T
        infile.close()
        
        dta = [dta1, dta2, dta3, dta4, dta5, dta6, dta7, dta8, dta9, dta10, dta11, dta12]
        dta = np.transpose(dta, (3, 2, 1, 0))
        dta = dta[:,63:129,:,:]
        
    else:
        
        infile = open(data_dir2 + '/TS/dts_monthly' + np.str(k + 1) + '_1955', 'rb')
        dts = np.array(pickle.load(infile))
        infile.close()
        dts = dts[:,63:129,:]
        
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
        dta = dta[:,63:129,:,:]
    
    # Calculate the change in global mean surface temperature
    dts_globalmean = np.nansum(np.nansum(np.multiply(np.nanmean(dts,axis=2), weight),axis=1),axis=0)/np.nansum(weight)
    dts_globalmean_mon = np.nansum(np.nansum(np.multiply(dts, weight_mon),axis=1),axis=0)
    
    for i in range(len(dts_globalmean_mon)):
        dts_globalmean_mon[i] = dts_globalmean_mon[i]/np.nansum(weight)

    # Multiply monthly mean TS change by the TS kernels (function of lat, lon, month) (units W/m2)
    dLW_ts = ts_kernel * dts  
    dLW_ts_cs = ts_kernel_clearsky * dts

    dta = dta * (p >= p_tropopause)

    # Convolve air temperature kernel with air temperature change
    ta_kernel = ma.filled(ta_kernel, fill_value=np.nan)
    ta_kernel_clearsky = ma.filled(ta_kernel_clearsky, fill_value=np.nan)
    dta = ma.filled(dta, fill_value=np.nan)
    dLW_ta = np.squeeze(np.sum(ta_kernel * dta, 2))
    dLW_ta_cs = np.squeeze(np.sum(ta_kernel_clearsky * dta, 2))

    # Add the surface and air temperature response; Take the annual average and global area average 
    #dLW_t_globalmean = np.nansum(np.nansum(np.nanmean(-dLW_ta - dLW_ts, 2) * weight, 1), 0)
    
    ##### Albedo feedback
    
    infile = open(data_dir2 + "/FSNS/FSNS_Base_monthly" + np.str(k + 1) + '_1955', 'rb')
    SW_sfc_net_1 = np.array(pickle.load(infile))
    infile.close()
   
    infile = open(data_dir2 + "/FSNS/dFSNS_monthly" + np.str(k + 1) + '_1955', 'rb')
    SW_sfc_net_2 = np.array(pickle.load(infile)) + SW_sfc_net_1
    infile.close()     
    SW_sfc_net_1 = SW_sfc_net_1[:,63:129,:]
    SW_sfc_net_2 = SW_sfc_net_2[:,63:129,:]
    
    if k == 0:
   
        infile = open(data_dir2 + "/FSDS/FSDS_Base_monthly" + np.str(len(filenames1a)) + '_1955', 'rb')
        SW_sfc_down_1 = np.array(pickle.load(infile))
        infile.close()
        
        infile = open(data_dir2 + "/FSDS/dFSDS_monthly" + np.str(len(filenames1a)) + '_1955', 'rb')
        SW_sfc_down_2 = np.array(pickle.load(infile)) + SW_sfc_down_1
        infile.close()
        SW_sfc_down_1 = SW_sfc_down_1[:,63:129,:]
        SW_sfc_down_2 = SW_sfc_down_2[:,63:129,:]
        
    else:
        
        infile = open(data_dir2 + "/FSDS/FSDS_Base_monthly" + np.str(k) + '_1955', 'rb')
        SW_sfc_down_1 = np.array(pickle.load(infile))
        infile.close()
        
        infile = open(data_dir2 + "/FSDS/dFSDS_monthly" + np.str(k) + '_1955', 'rb')
        SW_sfc_down_2 = np.array(pickle.load(infile)) + SW_sfc_down_1
        infile.close()
        SW_sfc_down_1 = SW_sfc_down_1[:,63:129,:]
        SW_sfc_down_2 = SW_sfc_down_2[:,63:129,:]
    
    alb1 = np.squeeze(1-SW_sfc_net_1/SW_sfc_down_1)
    alb1 = ma.filled(alb1, 0) 
    
    alb2 = np.squeeze(1-SW_sfc_net_2/SW_sfc_down_2)
    alb2 = ma.filled(alb2, 0) 
    
    dalb = (alb2 - alb1) * 100
    
    dSW_alb = alb_kernel * dalb
    dSW_alb_cs = alb_kernel_clearsky * dalb
    
    #dSW_alb_globalmean = np.nansum(np.nansum(np.nanmean(dSW_alb, 2) * weight, 1), 0)
    
    ####### Water vapor feedback
    
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
                
        q0 = [q01, q02, q03, q04, q05, q06, q07, q08, q09, q010, q011, q012]
        q0 = np.transpose(q0, (1, 2, 3, 0))
        
        dq = [dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12]
        dq = np.transpose(dq, (1, 2, 3, 0))
        
        t1 = [t11, t12, t13, t14, t15, t16, t17, t18, t19, t110, t111, t112]
        t1 = np.transpose(t1, (1, 2, 3, 0))
        
        q0 = q0[:,63:129,:,:]
        dq = dq[:,63:129,:,:]
        t1 = t1[:,63:129,:,:]

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
        
        q0 = [q01, q02, q03, q04, q05, q06, q07, q08, q09, q010, q011, q012]
        q0 = np.transpose(q0, (1, 2, 3, 0))
        
        dq = [dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12]
        dq = np.transpose(dq, (1, 2, 3, 0))
        
        t1 = [t11, t12, t13, t14, t15, t16, t17, t18, t19, t110, t111, t112]
        t1 = np.transpose(t1, (1, 2, 3, 0))
        
        q0 = q0[:,63:129,:,:]
        dq = dq[:,63:129,:,:]
        t1 = t1[:,63:129,:,:]
        
    #qs1 = calcsatspechum(t1,p)
    #qs2 = calcsatspechum(t1+dta,p)
    #dqsdt = (qs2 - qs1)/dta
    #rh = q1/qs1
    #dqdt = rh * dqsdt
    
    qs1 = calcsatspechum(t1, p) # g/kg
    qs2 = calcsatspechum(t1 + dta, p) # g/kg
    
    dta[dta == 0] = np.nan
    dqsdt = (qs2 - qs1)/dta
 
    rh = 1000 * q0/qs1
    dqdt = rh * dqsdt
    
    # Normalize kernels by the change in moisture for 1 K warming at
    # constant RH (linear)
    
    q_LW_kernel_2 = q_LW_kernel/dqdt
    q_SW_kernel_2 = q_SW_kernel/dqdt
    q_LW_kernel_clearsky_2 = q_LW_kernel_clearsky/dqdt
    q_SW_kernel_clearsky_2 = q_SW_kernel_clearsky/dqdt
    
    # Mask out the stratosphere
    dq = dq * (p >= p_tropopause)
    
    # Convolve moisture kernel with change in moisture
    dLW_q = np.squeeze(np.nansum(q_LW_kernel_2 * dq, 2))
    dSW_q = np.squeeze(np.nansum(q_SW_kernel_2 * dq, 2))
    dLW_q_cs = np.squeeze(np.nansum(q_LW_kernel_clearsky_2 * dq, 2))
    dSW_q_cs = np.squeeze(np.nansum(q_SW_kernel_clearsky_2 * dq, 2))
    
    # Add the LW and SW responses. Note the sign convention difference
    # between LW and SW!
    #dR_q_globalmean=np.nansum(np.nansum(np.nanmean(-dLW_q + dSW_q, 2) * weight, 1), 0)
    
    ##### Change in Cloud Radiative Effect (CRE) 
    
    infile = open(data_dir2 + "/FSNT/dFSNT_monthly" + np.str(k + 1) + '_1955', 'rb')
    d_sw = np.array(pickle.load(infile))
    infile.close() 
    
    infile = open(data_dir2 + "/FSNTC/dFSNTC_monthly" + np.str(k + 1) + '_1955', 'rb')
    d_sw_cs = np.array(pickle.load(infile))
    infile.close() 
    
    infile = open(data_dir2 + "/FLNT/dFLNT_monthly" + np.str(k + 1) + '_1955', 'rb')
    d_lw = np.array(pickle.load(infile))
    infile.close() 
    
    infile = open(data_dir2 + "/FLNTC/dFLNTC_monthly" + np.str(k + 1) + '_1955', 'rb')
    d_lw_cs = np.array(pickle.load(infile))
    infile.close() 
    
    d_cre_sw = d_sw_cs - d_sw
    d_cre_lw = d_lw_cs - d_lw
    d_cre_sw = d_cre_sw[:,63:129,:]
    d_cre_lw = d_cre_lw[:,63:129,:]
    
    #### Cloud masking of radiative forcing
    #ghgfile = 'forcing/ghg.forcing.nc'
    #ds10 = nc.Dataset(data_dir + "/" + ghgfile)
    #sw = ds10['FSNT'][:].T
    #sw_cs = ds10['FSNTC'][:].T
    #lw =  ds10['FLNT'][:].T
    #lw_cs = ds10['FLNTC'][:].T
    #ghg_sw = sw_cs - sw
    #ghg_lw = lw_cs - lw
    
    #aerosolfile='forcing/aerosol.forcing.nc'
    #ds11 = nc.Dataset(data_dir + "/" + aerosolfile)
    #sw = ds11['FSNT'][:].T
    #sw_cs = ds11['FSNTC'][:].T
    #lw = ds11['FLNT'][:].T
    #lw_cs = ds11['FLNTC'][:].T
    
    #aerosol_sw = sw_cs - sw
    #aerosol_lw = lw_cs - lw
    
    #cloud_masking_of_forcing_sw = aerosol_sw + ghg_sw
    #cloud_masking_of_forcing_lw = aerosol_lw + ghg_lw
    
    ###### Cloud feedback. 
    ### CRE + cloud masking of radiative forcing + corrections for each feedback
    
    dLW_cloud = -d_cre_lw + cloud_masking_of_forcing_lw + (dLW_q_cs - dLW_q) + (dLW_ta_cs - dLW_ta) + (dLW_ts_cs - dLW_ts)
    dSW_cloud = -d_cre_sw + cloud_masking_of_forcing_sw + (dSW_q_cs - dSW_q) + (dSW_alb_cs - dSW_alb)
    
    #print(np.nansum(np.nansum(np.nanmean(-d_cre_lw, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(cloud_masking_of_forcing_lw, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dLW_q_cs, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dLW_q, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dLW_ta_cs, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dLW_ta, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dLW_ts_cs, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dLW_ts, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(-d_cre_sw, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(cloud_masking_of_forcing_sw, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dSW_q_cs, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dSW_q, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dSW_alb_cs, 2) * weight, 1), 0))
    #print(np.nansum(np.nansum(np.nanmean(dSW_alb, 2) * weight, 1), 0))
    
    # Take global and annual averages
    dLW_cloud_globalmean = np.nansum(np.nansum(np.nanmean(-dLW_cloud, 2) * weight, 1), 0)/np.nansum(weight)
    dLW_cloud_globalmean_mon = np.nansum(np.nansum(-dLW_cloud * weight_mon, 1), 0)
    
    for i in range(len(dLW_cloud_globalmean_mon)):
        dLW_cloud_globalmean_mon[i] = dLW_cloud_globalmean_mon[i]/np.nansum(weight)
    
    dSW_cloud_globalmean = np.nansum(np.nansum(np.nanmean(dSW_cloud, 2) * weight, 1), 0)/np.nansum(weight)
    dSW_cloud_globalmean_mon = np.nansum(np.nansum(dSW_cloud * weight_mon, 1), 0)
    
    for i in range(len(dSW_cloud_globalmean_mon)):
        dSW_cloud_globalmean_mon[i] = dSW_cloud_globalmean_mon[i]/np.nansum(weight)
    
    # Divide by global, annual mean temperature change to get W/m2/K
    LWc_feedback.append(dLW_cloud_globalmean/dts_globalmean)
    LWc_feedback_mon.append(dLW_cloud_globalmean_mon/dts_globalmean_mon)
    LWc_feedback_flux.append(dLW_cloud_globalmean)
    LWc_feedback_flux_mon.append(dLW_cloud_globalmean_mon)
    
    SWc_feedback.append(dSW_cloud_globalmean/dts_globalmean)
    SWc_feedback_mon.append(dSW_cloud_globalmean_mon/dts_globalmean_mon)
    SWc_feedback_flux.append(dSW_cloud_globalmean)
    SWc_feedback_flux_mon.append(dSW_cloud_globalmean_mon)
    
    print('LW Cloud Feedback: ' + np.str(LWc_feedback[k]) + ' W m^-2 K^-1')
    print('SW Cloud Feedback: ' + np.str(SWc_feedback[k]) + ' W m^-2 K^-1')
    
outfile = open('Tropics_LW_Cloud_Feedback_1955', 'wb')
pickle.dump(LWc_feedback, outfile)
outfile.close()

outfile = open('Tropics_LW_Cloud_Feedback_Monthly_1955', 'wb')
pickle.dump(LWc_feedback_mon, outfile)
outfile.close()

outfile = open('Tropics_LW_Cloud_Feedback_Flux_1955', 'wb')
pickle.dump(LWc_feedback_flux, outfile)
outfile.close()

outfile = open('Tropics_LW_Cloud_Feedback_Flux_Monthly_1955', 'wb')
pickle.dump(LWc_feedback_flux_mon, outfile)
outfile.close()

outfile = open('Tropics_SW_Cloud_Feedback_1955', 'wb')
pickle.dump(SWc_feedback, outfile)
outfile.close()

outfile = open('Tropics_SW_Cloud_Feedback_Monthly_1955', 'wb')
pickle.dump(SWc_feedback_mon, outfile)
outfile.close()

outfile = open('Tropics_SW_Cloud_Feedback_Flux_1955', 'wb')
pickle.dump(SWc_feedback_flux, outfile)
outfile.close()

outfile = open('Tropics_SW_Cloud_Feedback_Flux_Monthly_1955', 'wb')
pickle.dump(SWc_feedback_flux_mon, outfile)
outfile.close()

np.savetxt('Tropics_LW_Cloud_Feedback_1955.csv', LWc_feedback, delimiter = ',')
np.savetxt('Tropics_LW_Cloud_Feedback_Monthly_1955.csv', LWc_feedback_mon, delimiter = ',')

np.savetxt('Tropics_LW_Cloud_Feedback_Flux_1955.csv', LWc_feedback_flux, delimiter = ',')
np.savetxt('Tropics_LW_Cloud_Feedback_Flux_Monthly_1955.csv', LWc_feedback_flux_mon, delimiter = ',')

np.savetxt('Tropics_SW_Cloud_Feedback_1955.csv', SWc_feedback, delimiter = ',')
np.savetxt('Tropics_SW_Cloud_Feedback_Monthly_1955.csv', SWc_feedback_mon, delimiter = ',')

np.savetxt('Tropics_SW_Cloud_Feedback_Flux_1955.csv', SWc_feedback_flux, delimiter = ',')
np.savetxt('Tropics_SW_Cloud_Feedback_Flux_Monthly_1955.csv', SWc_feedback_flux_mon, delimiter = ',')

print('done')
    
