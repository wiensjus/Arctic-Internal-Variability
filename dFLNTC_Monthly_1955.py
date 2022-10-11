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

data_dir = "/export/data/CESM-LENS/FLNTC/"
filenames2 = os.listdir(data_dir)

filenames = []
for i in range(len(filenames2)):
    if '.nc' in filenames2[i]:
        filenames.append(filenames2[i])

outfile = open('FLNTC_monthly_filenames', 'wb')
pickle.dump(filenames, outfile)
outfile.close()

data = []
for i in range(len(filenames)):
    data.append(nc.Dataset(data_dir + "/" + filenames[i]))
    
FLNTC = []

for i in range(len(data)):
    FLNTC.append(np.array(data[i]['FLNTC'][:]))
    
dFLNTC = []
temp = []
temp2 = []

for k in range(len(data)):
    print(k + 1)
    for i in range(np.shape(FLNTC[0])[2]):
        for j in range(np.shape(FLNTC[0])[1]):
            for m in range(12):
                if k == (len(data) - 1):
                    temp.append(np.nanmean(FLNTC[k][1632+m:1872+m:12,j,i]) - np.nanmean(FLNTC[k][1260+m:1500+m:12,j,i]))
                else:
                    temp.append(np.nanmean(FLNTC[k][792+m:1032+m:12,j,i]) - np.nanmean(FLNTC[k][420+m:660+m:12,j,i]))
            temp2.append(temp)
            temp = []
        dFLNTC.append(temp2)
        temp2 = []
    
    if k == (len(data) - 1):
        
        outfile = open('dFLNTC_monthly1_1955', 'wb')
        pickle.dump(dFLNTC, outfile)
        outfile.close()
        dFLNTC = []
        
    else:
        
        outfile = open('dFLNTC_monthly' + np.str(k + 2) + '_1955', 'wb')
        pickle.dump(dFLNTC, outfile)
        outfile.close()
        dFLNTC = []
    
done = [1,2,3]
outfile = open('done', 'wb')
pickle.dump(done, outfile)
outfile.close()
    
print('done')
