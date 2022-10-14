# -*- coding: utf-8 -*-
"""
Created on Fri May 28 13:04:24 2021

@author: Justin
"""

import numpy as np
import numpy.matlib
import os
import pickle

##### Cloud feedback

# File with the changes in climate: (ts, temp) (TS,T,Q)
changefile = 'demodata/changefields.nc'

basefile = 'demodata/basefields.nc'
###############################################

# Read in coordinate info
data_dir = "/export/data/CESM-LENS/proc/cam5-kernels/"
data_dir2 = "/export/data/CESM-LENS/"

#T Data
filenames1b = os.listdir(data_dir2 + "/" + 'T')

filenames1a = []
for i in range(len(filenames1b)):
    if '.nc' in filenames1b[i]:
        filenames1a.append(filenames1b[i])
        
LWc_warming = []
LWc_warming_mon = []
temp1 = []

LWc_warming_arc = []
LWc_warming_arc_mon = []
temp1arc = []

LWc_warming_trop = []
LWc_warming_trop_mon = []
temp1trop = []

SWc_warming = []
SWc_warming_mon = []
temp1b = []

SWc_warming_arc = []
SWc_warming_arc_mon = []
temp1barc = []

SWc_warming_trop = []
SWc_warming_trop_mon = []
temp1btrop = []
        
infile = open('LW_Cloud_Feedback_Flux_1955', 'rb')
LWc_feedback_flux = np.array(pickle.load(infile))
infile.close()

infile = open('LW_Cloud_Feedback_Flux_Monthly_1955', 'rb')
LWc_feedback_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_LW_Cloud_Feedback_Flux_1955', 'rb')
LWc_feedback_arc_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_LW_Cloud_Feedback_Flux_Monthly_1955', 'rb')
LWc_feedback_arc_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_LW_Cloud_Feedback_Flux_1955', 'rb')
LWc_feedback_trop_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_LW_Cloud_Feedback_Flux_Monthly_1955', 'rb')
LWc_feedback_trop_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('SW_Cloud_Feedback_Flux_1955', 'rb')
SWc_feedback_flux = np.array(pickle.load(infile))
infile.close()

infile = open('SW_Cloud_Feedback_Flux_Monthly_1955', 'rb')
SWc_feedback_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_SW_Cloud_Feedback_Flux_1955', 'rb')
SWc_feedback_arc_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_SW_Cloud_Feedback_Flux_Monthly_1955', 'rb')
SWc_feedback_arc_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_SW_Cloud_Feedback_Flux_1955', 'rb')
SWc_feedback_trop_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_SW_Cloud_Feedback_Flux_Monthly_1955', 'rb')
SWc_feedback_trop_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/proc/Temperature_Feedback_1955/Planck_Feedback_1955', 'rb')
planck_feedback =  np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/proc/Temperature_Feedback_1955/Planck_Feedback_Monthly_1955', 'rb')
planck_feedback_mon =  np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/proc/Temperature_Feedback_1955/Arctic_Planck_Feedback_1955', 'rb')
planck_feedback_arc = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/proc/Temperature_Feedback_1955/Arctic_Planck_Feedback_Monthly_1955', 'rb')
planck_feedback_arc_mon = np.array(pickle.load(infile))
infile.close()

for i in range(len(filenames1a)):
    print(i)
    
    LWc_warming.append(LWc_feedback_flux[i]/abs(planck_feedback[i]))
    LWc_warming_arc.append(LWc_feedback_arc_flux[i]/abs(planck_feedback[i]))
    LWc_warming_trop.append(LWc_feedback_trop_flux[i]/abs(planck_feedback[i]))
    
    SWc_warming.append(SWc_feedback_flux[i]/abs(planck_feedback[i]))
    SWc_warming_arc.append(SWc_feedback_arc_flux[i]/abs(planck_feedback[i]))
    SWc_warming_trop.append(SWc_feedback_trop_flux[i]/abs(planck_feedback[i]))
    
    for j in range(12):
        temp1.append(LWc_feedback_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp1arc.append(LWc_feedback_arc_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp1trop.append(LWc_feedback_trop_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        
        temp1b.append(SWc_feedback_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp1barc.append(SWc_feedback_arc_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp1btrop.append(SWc_feedback_trop_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        
    LWc_warming_mon.append(temp1)
    LWc_warming_arc_mon.append(temp1arc)
    LWc_warming_trop_mon.append(temp1arc)
    
    SWc_warming_mon.append(temp1b)
    SWc_warming_arc_mon.append(temp1barc)
    SWc_warming_trop_mon.append(temp1barc)
    
    temp1 = []
    temp1arc = []
    temp1trop = []
    
    temp1b = []
    temp1barc = []
    temp1btrop = []
    
    print(LWc_warming[i])
    
outfile = open('LW_Cloud_Feedback_Warming_1955', 'wb')
pickle.dump(LWc_warming, outfile)
outfile.close()

outfile = open('LW_Cloud_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(LWc_warming_mon, outfile)
outfile.close()

outfile = open('Arctic_LW_Cloud_Feedback_Warming_1955', 'wb')
pickle.dump(LWc_warming_arc, outfile)
outfile.close()

outfile = open('Arctic_LW_Cloud_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(LWc_warming_arc_mon, outfile)
outfile.close()

outfile = open('Tropics_LW_Cloud_Feedback_Warming_1955', 'wb')
pickle.dump(LWc_warming_trop, outfile)
outfile.close()

outfile = open('Tropics_LW_Cloud_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(LWc_warming_trop_mon, outfile)
outfile.close()

outfile = open('SW_Cloud_Feedback_Warming_1955', 'wb')
pickle.dump(SWc_warming, outfile)
outfile.close()

outfile = open('SW_Cloud_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(SWc_warming_mon, outfile)
outfile.close()

outfile = open('Arctic_SW_Cloud_Feedback_Warming_1955', 'wb')
pickle.dump(SWc_warming_arc, outfile)
outfile.close()

outfile = open('Arctic_SW_Cloud_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(SWc_warming_arc_mon, outfile)
outfile.close()

outfile = open('Tropics_SW_Cloud_Feedback_Warming_1955', 'wb')
pickle.dump(SWc_warming_trop, outfile)
outfile.close()

outfile = open('Tropics_SW_Cloud_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(SWc_warming_trop_mon, outfile)
outfile.close()

np.savetxt('LW_Cloud_Feedback_Warming_1955.csv', LWc_warming, delimiter = ',')
np.savetxt('LW_Cloud_Feedback_Warming_Monthly_1955.csv', LWc_warming_mon, delimiter = ',')

np.savetxt('Arctic_LW_Cloud_Feedback_Warming_1955.csv', LWc_warming_arc, delimiter = ',')
np.savetxt('Arctic_LW_Cloud_Feedback_Warming_Monthly_1955.csv', LWc_warming_arc_mon, delimiter = ',')

np.savetxt('Tropics_LW_Cloud_Feedback_Warming_1955.csv', LWc_warming_trop, delimiter = ',')
np.savetxt('Tropics_LW_Cloud_Feedback_Warming_Monthly_1955.csv', LWc_warming_trop_mon, delimiter = ',')

np.savetxt('SW_Cloud_Feedback_Warming_1955.csv', SWc_warming, delimiter = ',')
np.savetxt('SW_Cloud_Feedback_Warming_Monthly_1955.csv', SWc_warming_mon, delimiter = ',')

np.savetxt('Arctic_SW_Cloud_Feedback_Warming_1955.csv', SWc_warming_arc, delimiter = ',')
np.savetxt('Arctic_SW_Cloud_Feedback_Warming_Monthly_1955.csv', SWc_warming_arc_mon, delimiter = ',')

np.savetxt('Tropics_SW_Cloud_Feedback_Warming_1955.csv', SWc_warming_trop, delimiter = ',')
np.savetxt('Tropics_SW_Cloud_Feedback_Warming_Monthly_1955.csv', SWc_warming_trop_mon, delimiter = ',')

print('done')
