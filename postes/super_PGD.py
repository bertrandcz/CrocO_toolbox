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

import netCDF4  # @UnresolvedImport

import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET


def super_PGD(massifs, year, xpidforcing, obsPath= os.environ['CROCOPATH']  + '/s2m/postes/obs_csv/bd-clim_ML/OBS_2009080100_2020080123_types13456789.csv', vconf = None):

    if vconf is None:
        vconf = '_'.join(['postes'] + list(map(str, massifs)) + ['csv'])
    rootDir = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/spinup/'
    pgdPath = rootDir + '/pgd/PGD_' + vconf + '.nc'
    superpgdPath = rootDir + '/pgd/super_PGD_' + vconf + '.nc'

    forcing = os.environ['VORTEXPATH'] + '/safran/' + vconf + '/' + xpidforcing + '/mb0000/meteo/FORCING_{0}080106_{1}080106.nc'.format(year, year + 1)

    if os.path.exists(pgdPath):
        if os.path.exists(superpgdPath):
            os.remove(superpgdPath)
        shutil.copyfile(pgdPath, superpgdPath)
    else:
        print(pgdPath)
        raise Exception
    # A. read the stations in the FORCING
    f = netCDF4.Dataset(forcing, 'r')
    superpgd = netCDF4.Dataset(superpgdPath, 'a')
    station = f.variables['station']
    print('station', station[:])
    # B. read the type of post in the observations file
    df = pd.read_csv(obsPath, delimiter = ';')

    # B.1. PREPARE the Dataframe
    df1 = df.pivot(index = 'Date UTC', columns = 'NUMPOST', values = 'TYPE_NIVO')
    df2 = df.pivot(index = 'Date UTC', columns = 'NUMPOST', values = 'TYPE_POSTE_ACTUEL')
    df3 = df.pivot(index = 'Date UTC', columns = 'NUMPOST', values = 'RESEAU_POSTE_ACTUEL')

    # B.2 read the type
    type_df = df1.mean()
    typestat = np.array([int(type_df[s]) for s in station[:]])

    type_df = df2.mean()
    typestat_poste_actuel = np.array([int(type_df[s]) for s in station[:]])
    type_df = df3.mean()
    typestat_reseau = np.array([int(type_df[s]) for s in station[:]])

    # C. read the massif number in the metadata file

    tree = ET.parse(os.environ['SNOWTOOLS_CEN'] + '/DATA/METADATA.xml')
    root = tree.getroot()
    statmassif = []
    for stat in station[:]:
        for site in root[1].findall('Site'):
            if str(stat) == site.find('number').text.strip().lstrip('0'):
                statmassif.append(site.find('massif').text.strip())

    # D. WRITE
    statPgd = superpgd.createVariable('station', 'int32', ('Number_of_points',), fill_value = station._FillValue)
    type_nivo = superpgd.createVariable('type_nivo', 'int32', ('Number_of_points',), fill_value = station._FillValue)
    type_poste_actuel = superpgd.createVariable('type_poste_actuel', 'int32', ('Number_of_points',), fill_value = station._FillValue)
    reseau_poste_actuel = superpgd.createVariable('reseau_poste_actuel', 'int32', ('Number_of_points',), fill_value = station._FillValue)
    massifPgd = superpgd.createVariable('massif', 'int32', ('Number_of_points',), fill_value = station._FillValue)
    statPgd[:] = station[:]
    type_nivo[:] = typestat
    type_poste_actuel[:] = typestat_poste_actuel
    reseau_poste_actuel[:] = typestat_reseau
    massifPgd[:] = statmassif
    superpgd.close()
    return superpgdPath
