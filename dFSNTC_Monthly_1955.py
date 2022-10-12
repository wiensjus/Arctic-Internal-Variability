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

data_dir = "/export/data/CESM-LENS/FSNTC/"
filenames2 = os.listdir(data_dir)

filenames = []
for i in range(len(filenames2)):
    if '.nc' in filenames2[i]:
        filenames.append(filenames2[i])

outfile = open('FSNTC_monthly_filenames', 'wb')
pickle.dump(filenames, outfile)
outfile.close()

data = []
for i in range(len(filenames)):
    data.append(nc.Dataset(data_dir + "/" + filenames[i]))
    
FSNTC = []

for i in range(len(data)):
    FSNTC.append(np.array(data[i]['FSNTC'][:]))
    
dFSNTC = []
temp = []
temp2 = []

for k in range(len(data)):
    print(k + 1)
    for i in range(np.shape(FSNTC[0])[2]):
        for j in range(np.shape(FSNTC[0])[1]):
            for m in range(12):
                if k == (len(data) - 1):
                    temp.append(np.nanmean(FSNTC[k][1632+m:1872+m:12,j,i]) - np.nanmean(FSNTC[k][1260+m:1500+m:12,j,i]))
                else:
                    temp.append(np.nanmean(FSNTC[k][792+m:1032+m:12,j,i]) - np.nanmean(FSNTC[k][420+m:660+m:12,j,i]))
            temp2.append(temp)
            temp = []
        dFSNTC.append(temp2)
        temp2 = []
        
    if k == (len(data) - 1):
        
        outfile = open('dFSNTC_monthly1_1955', 'wb')
        pickle.dump(dFSNTC, outfile)
        outfile.close()
        dFSNTC = []
    
    else:
    
        outfile = open('dFSNTC_monthly' + np.str(k + 2) + '_1955', 'wb')
        pickle.dump(dFSNTC, outfile)
        outfile.close()
        dFSNTC = []

done = [1,2,3]
outfile = open('done', 'wb')
pickle.dump(done, outfile)
outfile.close()

print('done')
