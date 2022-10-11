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

data_dir = "/export/data/CESM-LENS/FLNT/"
filenames2 = os.listdir(data_dir)

filenames = []
for i in range(len(filenames2)):
    if '.nc' in filenames2[i]:
        filenames.append(filenames2[i])

outfile = open('FLNT_monthly_filenames', 'wb')
pickle.dump(filenames, outfile)
outfile.close()

data = []
for i in range(len(filenames)):
    data.append(nc.Dataset(data_dir + "/" + filenames[i]))
    
FLNT = []

for i in range(len(data)):
    FLNT.append(np.array(data[i]['FLNT'][:]))
    
dFLNT = []
temp = []
temp2 = []

for k in range(len(data)):
    print(k + 1)
    for i in range(np.shape(FLNT[0])[2]):
        for j in range(np.shape(FLNT[0])[1]):
            for m in range(12):
                temp.append(np.nanmean(FLNT[k][792+m:1032+m:12,j,i]) - np.nanmean(FLNT[k][420+m:660+m:12,j,i]))
            temp2.append(temp)
            temp = []
        dFLNT.append(temp2)
        temp2 = []
    
    outfile = open('dFLNT_monthly' + np.str(k + 1) + '_1955', 'wb')
    pickle.dump(dFLNT, outfile)
    outfile.close()
    dFLNT = []
    
print('done')
