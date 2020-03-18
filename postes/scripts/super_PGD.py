# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 13 févr. 2020

@author: cluzetb
Add data (station number, type of station and massif ID) to the PGD.nc (from spinup), creating a super_PGD.nc
which will be used for the post-processing of data.

'''
import os
import shutil

import netCDF4

import numpy as np
import pandas as pd
from postes.explore_metadata import find_metadata
import xml.etree.ElementTree as ET


massifs = [8, 9, 12, 13, 15, 16]
year = 2013
xpid = 'ref2'
vconf = 'postes_8_9_12_13_15_16_csv'
rootDir = '/home/cluzetb/vortexpath/s2m/' + vconf + '/spinup/'
pgdPath = rootDir + '/pgd/PGD_' + vconf + '.nc'
superpgdPath = rootDir + '/pgd/super_PGD_' + vconf + '.nc'
forcing = '/home/cluzetb/vortexpath/safran/' + vconf + '/' + xpid + '/mb0000/meteo/FORCING_{0}080106_{1}080106.nc'.format(year, year + 1)
obsPath = '/home/cluzetb/assim_postes/OBS_2010080100_2020021023_types1346789.csv'
print(pgdPath)
if os.path.exists(pgdPath):
    if os.path.exists(superpgdPath):
        os.remove(superpgdPath)
    shutil.copyfile(pgdPath, superpgdPath)
else:
    raise Exception
# A. read the stations in the FORCING
f = netCDF4.Dataset(forcing, 'r')
superpgd = netCDF4.Dataset(superpgdPath, 'a')
station = f.variables['station']
print('station', station[:])
# B. read the type  in the observations file
df = pd.read_csv(obsPath, delimiter = ';')

# B.1. PREPARE the Dataframe
df = df.pivot(index = 'Date UTC', columns = 'NUMPOST', values = 'type')

# B.2 read the type
type_df = df.mean()
typestat = np.empty((len(station[:]),))
typestat = np.array([int(type_df[s]) for s in station[:]])

# C. read the massif number in the metadata file

tree = ET.parse('/home/cluzetb/snowtools_git/DATA/METADATA.xml')
root = tree.getroot()
statmassif = []
for stat in station[:]:
    for site in root[1].findall('Site'):
        if str(stat) == site.find('number').text.strip().lstrip('0'):
            statmassif.append(site.find('massif').text.strip())

# D. WRITE
statPgd = superpgd.createVariable('station', 'int32', ('Number_of_points',), fill_value = station._FillValue)
typePgd = superpgd.createVariable('type', 'int32', ('Number_of_points',), fill_value = station._FillValue)
massifPgd = superpgd.createVariable('massif', 'int32', ('Number_of_points',), fill_value = station._FillValue)
statPgd[:] = station[:]
typePgd[:] = typestat
massifPgd[:] = statmassif
superpgd.close()
