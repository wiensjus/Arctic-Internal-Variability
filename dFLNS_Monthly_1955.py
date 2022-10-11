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

data_dir = "/export/data/CESM-LENS/FLNS/"
filenames2 = os.listdir(data_dir)

filenames = []
for i in range(len(filenames2)):
    if '.nc' in filenames2[i]:
        filenames.append(filenames2[i])

outfile = open('FLNS_monthly_filenames', 'wb')
pickle.dump(filenames, outfile)
outfile.close()

data = []
for i in range(len(filenames)):
    data.append(nc.Dataset(data_dir + "/" + filenames[i]))
    
FLNS = []

for i in range(len(data)):
    FLNS.append(np.array(data[i]['FLNS'][:]))
    
dFLNS = []
temp = []
temp2 = []

for k in range(len(data)):
    print(k + 1)
    for i in range(np.shape(FLNS[0])[2]):
        for j in range(np.shape(FLNS[0])[1]):
            for m in range(12):
                temp.append(np.nanmean(FLNS[k][792+m:1032+m:12,j,i]) - np.nanmean(FLNS[k][420+m:660+m:12,j,i]))
            temp2.append(temp)
            temp = []
        dFLNS.append(temp2)
        temp2 = []
    
    outfile = open('dFLNS_monthly' + np.str(k + 1) + '_1955', 'wb')
    pickle.dump(dFLNS, outfile)
    outfile.close()
    dFLNS = []
    
print('done')
