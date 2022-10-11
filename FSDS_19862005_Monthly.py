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

data_dir = "/export/data/CESM-LENS/FSDS/"
filenames2 = os.listdir(data_dir)

filenames = []
for i in range(len(filenames2)):
    if '.nc' in filenames2[i]:
        filenames.append(filenames2[i])

outfile = open('FSDS_monthly_filenames', 'wb')
pickle.dump(filenames, outfile)
outfile.close()

data = []
for i in range(len(filenames)):
    data.append(nc.Dataset(data_dir + "/" + filenames[i]))
    
FSDS = []

for i in range(len(data)):
    FSDS.append(np.array(data[i]['FSDS'][:]))
    
FSDS_19812005 = []
temp = []
temp2 = []

for k in range(len(data)):
    print(k + 1)
    for i in range(np.shape(FSDS[0])[2]):
        for j in range(np.shape(FSDS[0])[1]):
            for m in range(12):
                if k == (len(filenames) - 1):
                    temp.append(np.mean(FSDS[k][1632+m:1872+m:12,j,i]))
                else:
                    temp.append(np.mean(FSDS[k][792+m:1032+m:12,j,i]))
            temp2.append(temp)
            temp = []
        FSDS_19812005.append(temp2)
        temp2 = []
    
    if k == (len(filenames) - 1):
        outfile = open('FSDS_19862005_monthly1', 'wb')
        pickle.dump(FSDS_19812005, outfile)
        outfile.close()
        FSDS_19812005 = []
        
    else:
        outfile = open('FSDS_19862005_monthly' + np.str(k + 2), 'wb')
        pickle.dump(FSDS_19812005, outfile)
        outfile.close()
        FSDS_19812005 = []
    

print('done')
