# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 17 janv. 2020
@author: cluzetb
 -Extract posts forcings from the reanalysis and add MOCAGE to it.
 - extract observations from the csv file and write it into the corresponding format
 for the first part, inspired on Dev/article2/new_years_MOCAGE.py
 - extract oper simulation runs from pro at postes
'''

import os
import pickle

import netCDF4

import numpy as np
from postes.explore_metadata import slice_file_listpnb
from postes.explore_obs import find_common_massif_obsbase
from utilcrocO import todates
from utils.prosimu import prosimu

years = [2013,2014,2015,2016]
massifs = [8, 9, 12, 13, 15, 16]
vconf = 'postes_' + '_'.join(list(map(str, massifs))) + '_csv'
xpidobs = 'all'
forcingOutDir = os.environ['VORTEXPATH'] + '/safran/' + vconf + '/ref2/mb0000/meteo/'
obsOutDir = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/obs/' + xpidobs + '/'
if not os.path.exists(forcingOutDir):
    os.makedirs(forcingOutDir)
if not os.path.exists(obsOutDir):
    os.makedirs(obsOutDir)

# auxiliary files
metadata = '/home/cluzetb/snowtools_git/DATA/METADATA.xml'
obscsvPath = '/home/cluzetb/assim_postes/OBS_2010080100_2020021023_types1346789.csv'

forcMOC = '/home/cluzetb/vortexpath/safran/12/FORCING_FRANCOIS_MOCAGE_Lautaret_2013_2018.nc'
forcSAFIMP_ref = '/home/cluzetb/vortexpath/safran/12/FORCING_reference/FORCING_2016080106_2017080106.nc'

# preliminary loadings

df, _, listpnb, listpname = find_common_massif_obsbase(obscsvPath, massifs = massifs)

# read impurities variables from the impurities file.
lvar = ['IMPWET1', 'IMPWET2', 'IMPDRY1', 'IMPDRY2']
dictimp = dict()
dimpur = netCDF4.Dataset(forcMOC, 'r')
for var in lvar:
    dictimp[var] = dimpur.variables[var][:]
imptime, impdates = todates(dimpur)


# read reference forcing (for standard impurities variables attributes and units)
dforc_ref = netCDF4.Dataset(forcSAFIMP_ref, 'r')
print(listpnb)
for year in years:
    # 1. extract the forcing
    if False:
        forcingPathIn = '/era40/vortex/s2m/postes/reanalysis/meteo/FORCING_{0}080106_{1}080106.nc'.format(year, year + 1)
        forcingOut = forcingOutDir + 'FORCING_extract_{0}080106_{1}080106.nc'.format(year, year + 1)

        # Slice on the selected posts
        if True:
            print('reading in', forcingPathIn)
            print('writing into', forcingOut)
            done = slice_file_listpnb(forcingPathIn, forcingOut, listpnb)

        # 2. add impurities
        dforc_saf = netCDF4.Dataset(forcingOut, 'r')
        pathnew = forcingOutDir + 'FORCING_{0}080106_{1}080106.nc'.format(year, year + 1)
        if os.path.exists(pathnew):
            os.remove(pathnew)
        dforc_new = netCDF4.Dataset(pathnew, 'w', format='NETCDF4_CLASSIC')
        # copy it and create the new variables.
        # BC from maskSentinel2.py and SemiDistributed.py
        for dimName, dim in dforc_saf.dimensions.items():
            dforc_new.createDimension(dimName, len(dim) if not dim.isunlimited() else None)
        # reversed to start copying time.
        for varname, var in reversed(dforc_saf.variables.items()):
            print(varname)
            # varname = var[0]
            # var = var[1]
            # copy variables and attributes
            gg = dforc_new.createVariable(varname, var.datatype, var.dimensions)
            # copy attributes
            gg.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
            gg[:] = var[:]
        # create new variables for impurities (only dims and attributes): copy it from an already created (and working) FORCING)
        dictNew = dict()

        # variables should be float in the new format.
        for varname in lvar:
            var = dforc_ref.variables[varname]
            gg = dforc_new.createVariable(varname, 'float32', var.dimensions)
            # copy attributes
            for k in var.ncattrs():
                print(k, var.getncattr(k))
            print(gg.dtype)
            gg.setncatts({k: var.getncattr(k) for k in ['units', 'long_name']})
            print('tut', dforc_saf.variables['CO2air']._FillValue)
            gg.setncatts({'_FillValue': dforc_saf.variables['CO2air']._FillValue})
            print('tot', gg._FillValue)

            dictNew[varname] = gg
        # write the impurities values : match the dates and expand in space
        saftime, safdates = todates(dforc_saf)
        maskTime = [i for i, d in enumerate(impdates) if d in safdates]
        for var in lvar:
            dictNew[var][:] = np.repeat(np.expand_dims(dictimp[var][maskTime], axis = 1), dforc_saf.dimensions['Number_of_points'].size, axis = 1)

        dforc_saf.close()
        dforc_new.close()

    # OPTIONAL : slice also the PRO from s2m reanalysis
    if True:
        pathPro = '/home/cluzetb/vortexpath/s2m/{0}/oper_{1}/mb0000/'.format(vconf, year, )
        if not os.path.exists(pathPro):
            os.makedirs(pathPro)
        proOut = pathPro + 'PRO_{0}080106_{1}080106.nc'.format(year, year + 1)
        ok = slice_file_listpnb('/era40/vortex/s2m/postes/reanalysis/pro/PRO_{0}080106_{1}080106.nc'.format(year, year + 1),
                                proOut, listpnb)
        pro = prosimu(proOut)
        depOper = np.expand_dims(pro.read('DSN_T_ISBA'), axis = 2)
        timeOper = pro.readtime()
        oper = dict()
        oper['DEP'] = depOper
        oper['time'] = timeOper
        pathOper = '/home/cluzetb/vortexpath/s2m/{0}/oper_{1}/crocO/beaufix/'.format(vconf, year)
        if not os.path.exists(pathOper):
            os.makedirs(pathOper)
        with open(pathOper + 'oper.pkl', 'wb') as f:
            pickle.dump(oper, f, protocol = pickle.HIGHEST_PROTOCOL)
dimpur.close()
dforc_ref.close()

# 2. prepare the observations
# columns of the df are already in listpnb order (same for the forcings/Pgd.nc generated from them)
if False:
    for dd in df.index.to_datetime():
        obsOut = obsOutDir + 'obs_' + xpidobs + '_' + vconf + '_'  + dd.strftime('%Y%m%d%H') + '.nc'

        # print(df.loc[dd])
        if os.path.exists(obsOut):
            os.remove(obsOut)
        print('new obs file:', obsOut)
        obs = netCDF4.Dataset(obsOut, 'w', data_model='NETCDF3_CLASSIC')
        obs.createDimension('Number_of_points', len(listpnb))
        dep = obs.createVariable('DEP', 'float', ('Number_of_points', ), fill_value=1e20)
        dep[:] = np.ma.masked_invalid(np.asarray(df.loc[dd]))
        obs.close()
if False:
    # create also a pickle of obs timeseries for evaluation

    obsTs = dict()
    obsTs['time'] = np.array(df.index.to_pydatetime())
    obsTs['DEP'] = np.array(df.values)
    # write it
    obsOutTs = obsOutDir + 'obs_' + xpidobs + '_' + vconf + '_'  + df.index.to_datetime()[0].strftime('%Y%m%d%H') +\
        '_' + df.index.to_datetime()[-1].strftime('%Y%m%d%H') + '.pkl'
    if os.path.exists(obsOutTs):
        os.remove(obsOutTs)
    with open(obsOutTs, 'wb') as f:
        pickle.dump(obsTs, f, protocol=pickle.HIGHEST_PROTOCOL)
