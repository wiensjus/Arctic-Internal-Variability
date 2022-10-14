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
        
alb_warming = []
alb_warming_mon = []
temp1 = []

alb_warming_arc = []
alb_warming_arc_mon = []
temp1arc = []

alb_warming_trop = []
alb_warming_trop_mon = []
temp1trop = []
        
infile = open('Temperature_Feedback_Flux_1955', 'rb')
alb_feedback_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Temperature_Feedback_Flux_Monthly_1955', 'rb')
alb_feedback_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_Temperature_Feedback_Flux_1955', 'rb')
alb_feedback_arc_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Arctic_Temperature_Feedback_Flux_Monthly_1955', 'rb')
alb_feedback_arc_flux_mon = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_Temperature_Feedback_Flux_1955', 'rb')
alb_feedback_trop_flux = np.array(pickle.load(infile))
infile.close()

infile = open('Tropics_Temperature_Feedback_Flux_Monthly_1955', 'rb')
alb_feedback_trop_flux_mon = np.array(pickle.load(infile))
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
    alb_warming.append(alb_feedback_flux[i]/abs(planck_feedback[i]))
    alb_warming_arc.append(alb_feedback_arc_flux[i]/abs(planck_feedback[i]))
    alb_warming_trop.append(alb_feedback_trop_flux[i]/abs(planck_feedback[i]))
    for j in range(12):
        temp1.append(alb_feedback_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp1arc.append(alb_feedback_arc_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
        temp1trop.append(alb_feedback_trop_flux_mon[i,j]/abs(planck_feedback_mon[i,j]))
    alb_warming_mon.append(temp1)
    alb_warming_arc_mon.append(temp1arc)
    alb_warming_trop_mon.append(temp1trop)
    temp1 = []
    temp1arc = []
    temp1trop = []
    print(alb_warming[i])
    
outfile = open('Temperature_Feedback_Warming_1955', 'wb')
pickle.dump(alb_warming, outfile)
outfile.close()

outfile = open('Temperature_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(alb_warming_mon, outfile)
outfile.close()

outfile = open('Arctic_Temperature_Feedback_Warming_1955', 'wb')
pickle.dump(alb_warming_arc, outfile)
outfile.close()

outfile = open('Arctic_Temperature_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(alb_warming_arc_mon, outfile)
outfile.close()

outfile = open('Tropics_Temperature_Feedback_Warming_1955', 'wb')
pickle.dump(alb_warming_trop, outfile)
outfile.close()

outfile = open('Tropics_Temperature_Feedback_Warming_Monthly_1955', 'wb')
pickle.dump(alb_warming_trop_mon, outfile)
outfile.close()

np.savetxt('Temperature_Feedback_Warming_1955.csv', alb_warming, delimiter = ',')
np.savetxt('Temperature_Feedback_Warming_Monthly_1955.csv', alb_warming_mon, delimiter = ',')

np.savetxt('Arctic_Temperature_Feedback_Warming_1955.csv', alb_warming_arc, delimiter = ',')
np.savetxt('Arctic_Temperature_Feedback_Warming_Monthly_1955.csv', alb_warming_arc_mon, delimiter = ',')

np.savetxt('Tropics_Temperature_Feedback_Warming_1955.csv', alb_warming_trop, delimiter = ',')
np.savetxt('Tropics_Temperature_Feedback_Warming_Monthly_1955.csv', alb_warming_trop_mon, delimiter = ',')

print('done')
