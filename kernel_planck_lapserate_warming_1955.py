# -*- coding: utf-8 -*-
"""
Created on Fri May 28 13:04:24 2021

@author: Justin
"""

import numpy as np
import numpy.matlib
import os
import pickle

# Read in coordinate info
data_dir = "/export/data/CESM-LENS/proc/cam5-kernels/"
data_dir2 = "/export/data/CESM-LENS/"

#T Data
filenames1b = os.listdir(data_dir2 + "/" + 'T')

filenames1a = []
for i in range(len(filenames1b)):
    if '.nc' in filenames1b[i]:
        filenames1a.append(filenames1b[i])
        
pla_warming = []
pla_warming_mon = []
lr_warming = []
lr_warming_mon = []
temp1 = []
temp2 = []

pla_warming_arc = []
pla_warming_arc_mon = []
lr_warming_arc = []
lr_warming_arc_mon = []
temp1arc = []
temp2arc = []

pla_warming_trop = []
pla_warming_trop_mon = []
lr_warming_trop = []
lr_warming_trop_mon = []
temp1trop = []
temp2trop = []

infile = open('/export/data/CESM-LENS/TS/dts_globalmean_annual_1955', 'rb')
dts_ann = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/TS/dts_arctic_annual_1955', 'rb')
dts_ann_arc = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/TS/dts_tropics_annual_1955', 'rb')
dts_ann_trop = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/TS/dts_globalmean_monthly_1955', 'rb')
dts_mon = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/TS/dts_arctic_monthly_1955', 'rb')
dts_mon_arc = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/TS/dts_tropics_monthly_1955', 'rb')
dts_mon_trop = np.array(pickle.load(infile))
infile.close()

infile = open('Lapse_Rate_Feedback_Flux_1955', 'rb')
lr_feedback_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Lapse_Rate_Feedback_Flux_Monthly_1955', 'rb')
lr_feedback_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_Lapse_Rate_Feedback_Flux_1955', 'rb')
lr_feedback_arc_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_Lapse_Rate_Feedback_Flux_Monthly_1955', 'rb')
lr_feedback_arc_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_Lapse_Rate_Feedback_Flux_1955', 'rb')
lr_feedback_trop_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_Lapse_Rate_Feedback_Flux_Monthly_1955', 'rb')
lr_feedback_trop_flux_mon = np.array(pickle.load(infile))
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

infile = open('/export/data/CESM-LENS/proc/Temperature_Feedback_1955/Tropics_Planck_Feedback_1955', 'rb')
planck_feedback_trop = np.array(pickle.load(infile))
infile.close()

infile = open('/export/data/CESM-LENS/proc/Temperature_Feedback_1955/Tropics_Planck_Feedback_Monthly_1955', 'rb')
planck_feedback_trop_mon = np.array(pickle.load(infile))
infile.close()

for i in range(len(filenames1a)):
    print(i)
    
    pla_feedback_flux = (planck_feedback[i] - planck_feedback[i]) * dts_ann[i]
    pla_feedback_arc_flux = (planck_feedback_arc[i] - planck_feedback[i]) * dts_ann_arc[i]
    pla_feedback_trop_flux = (planck_feedback_trop[i] - planck_feedback[i]) * dts_ann_trop[i]
    
    pla_warming.append(pla_feedback_flux/abs(planck_feedback[i]))
    pla_warming_arc.append(pla_feedback_arc_flux/abs(planck_feedback[i]))
    pla_warming_trop.append(pla_feedback_trop_flux/abs(planck_feedback[i]))
    
    lr_warming.append(lr_feedback_flux[i]/abs(planck_feedback[i]))
    lr_warming_arc.append(lr_feedback_arc_flux[i]/abs(planck_feedback[i]))
    lr_warming_trop.append(lr_feedback_trop_flux[i]/abs(planck_feedback[i]))
    
    print(lr_feedback_arc_flux[i], planck_feedback[i])

    for j in range(12):
        
        pla_feedback_flux_mon = (planck_feedback_mon[i,j] - planck_feedback_mon[i,j]) * dts_mon[i,j]
        pla_feedback_arc_flux_mon = (planck_feedback_arc_mon[i,j] - planck_feedback_mon[i,j]) * dts_mon_arc[i,j]
        pla_feedback_trop_flux_mon = (planck_feedback_trop_mon[i,j] - planck_feedback_mon[i,j]) * dts_mon_trop[i,j]
        
        temp1.append(pla_feedback_flux_mon/abs(planck_feedback_mon[i,j]))
        temp1arc.append(pla_feedback_arc_flux_mon/abs(planck_feedback_mon[i,j]))
        temp1trop.append(pla_feedback_trop_flux_mon/abs(planck_feedback_mon[i,j]))
        
        temp2.append(lr_feedback_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp2arc.append(lr_feedback_arc_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp2trop.append(lr_feedback_trop_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        
    pla_warming_mon.append(temp1)
    pla_warming_arc_mon.append(temp1arc)
    pla_warming_trop_mon.append(temp1trop)
    
    lr_warming_mon.append(temp2)
    lr_warming_arc_mon.append(temp2arc)
    lr_warming_trop_mon.append(temp2trop)
    
    temp1 = []
    temp1arc = []
    temp1trop = []
    temp2 = []
    temp2arc = []
    temp2trop = []
    #print(pla_warming[i])
    #print(lr_warming[i])
    
outfile = open('Planck_Feedback_Warming_1955', 'wb')
pickle.dump(pla_warming, outfile)
outfile.close()

outfile = open('Planck_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(pla_warming_mon, outfile)
outfile.close()

outfile = open('Arctic_Planck_Feedback_Warming_1955', 'wb')
pickle.dump(pla_warming_arc, outfile)
outfile.close()

outfile = open('Arctic_Planck_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(pla_warming_arc_mon, outfile)
outfile.close()

outfile = open('Tropics_Planck_Feedback_Warming_1955', 'wb')
pickle.dump(pla_warming_trop, outfile)
outfile.close()

outfile = open('Tropics_Planck_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(pla_warming_trop_mon, outfile)
outfile.close()

outfile = open('Lapse_Rate_Feedback_Warming_1955', 'wb')
pickle.dump(pla_warming, outfile)
outfile.close()

outfile = open('Lapse_Rate_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(pla_warming_mon, outfile)
outfile.close()

outfile = open('Arctic_Lapse_Rate_Feedback_Warming_1955', 'wb')
pickle.dump(pla_warming_arc, outfile)
outfile.close()

outfile = open('Arctic_Lapse_Rate_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(pla_warming_arc_mon, outfile)
outfile.close()

outfile = open('Tropics_Lapse_Rate_Feedback_Warming_1955', 'wb')
pickle.dump(pla_warming_trop, outfile)
outfile.close()

outfile = open('Tropics_Lapse_Rate_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(pla_warming_trop_mon, outfile)
outfile.close()

np.savetxt('Planck_Feedback_Warming_1955.csv', pla_warming, delimiter = ',')
np.savetxt('Planck_Feedback_Warming_Monthly_1955.csv', pla_warming_mon, delimiter = ',')

np.savetxt('Arctic_Planck_Feedback_Warming_1955.csv', pla_warming_arc, delimiter = ',')
np.savetxt('Arctic_Planck_Feedback_Warming_Monthly_1955.csv', pla_warming_arc_mon, delimiter = ',')

np.savetxt('Tropics_Planck_Feedback_Warming_1955.csv', pla_warming_trop, delimiter = ',')
np.savetxt('Tropics_Planck_Feedback_Warming_Monthly_1955.csv', pla_warming_trop_mon, delimiter = ',')

np.savetxt('Lapse_Rate_Feedback_Warming_1955.csv', lr_warming, delimiter = ',')
np.savetxt('Lapse_Rate_Feedback_Warming_Monthly_1955.csv', lr_warming_mon, delimiter = ',')

np.savetxt('Arctic_Lapse_Rate_Feedback_Warming_1955.csv', lr_warming_arc, delimiter = ',')
np.savetxt('Arctic_Lapse_Rate_Feedback_Warming_Monthly_1955.csv', lr_warming_arc_mon, delimiter = ',')

np.savetxt('Tropics_Lapse_Rate_Feedback_Warming_1955.csv', lr_warming_trop, delimiter = ',')
np.savetxt('Tropics_Lapse_Rate_Feedback_Warming_Monthly_1955.csv', lr_warming_trop_mon, delimiter = ',')

print('done')
