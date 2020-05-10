# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 24 févr. 2020

@author: cluzetb
create a masked version of the postes observations  to select a specific massif only.
'''

import os
import subprocess

from consts import CROCO
from postes.utilpostes import set_conf_everydate
from utilcrocO import Pgd


# # ####### PARAMETERS #####
# vconf = 'postes_12_csv'
# selMassif = [12]
# xpidobs = 'all'
# xpid = 'testmaskobs'
# startY = 2013
# endY = 2017
# assimilate_every = 7
# # in order to mask the obs that are not in the selected massif,
# # we use the classesId argument of crocO.py
def create_obs_massif_exclude_poste(vconf, selMassif, year, assimilate_every, alwaysX = []):
    # the classesId are in Python index (0 for the first class/station.)
    pgdPath = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/spinup/pgd/PGD_' + vconf + '.nc'
    pgd = Pgd(pgdPath)
    os.environ['CROCOPATH'] = os.environ['VORTEXPATH']
    classesId = [i for i, m in enumerate(pgd.massif) if m in selMassif]
    print(pgd.station[classesId])
    print(classesId)
    sensor_in = '_'.join(map(str, selMassif))
    sensors = []
    for i, stat in enumerate(pgd.station[classesId]):
        # prevent from excluding twice
        # classesIdK is the Ids of stations KEPT for assimilation
        classesIdK = list(set([cl for j, cl in enumerate(classesId) if j != i and str(pgd.station[cl]) not in alwaysX]))
        classesIdArg = ','.join(list(map(str, classesIdK)))
        listX = list(sorted(set([str(stat)] + alwaysX)))
        statX = '-'.join(listX)
        sensor = '_'.join(map(str, selMassif)) + '_X' + str(statX)
        sensors.append(sensor)
        xpidobs = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/obs/' + sensor
        pathConf = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/obs/' + sensor + '/conf/s2m_' + vconf + '_' + str(year) + '.ini'
        if not os.path.exists(pathConf):
            # generate a fake conf file with dates every 7days between startY and endY
            _ = set_conf_everydate(year, assimilate_every, pathConf, endY=year + 1)
            args = [

                'python3', CROCO + '/crocO.py',
                '--xpid', 'obs',
                '--todo', 'generobs',
                '--vconf', vconf,
                '-d', 'all',
                '--vars', 'DEP',
                '-o', 'testgener',
                '--classes_id', classesIdArg,
                '--sensor', sensor,
                '--sensor_in', sensor_in,
                '--archive_synth',
                '--pathConf', pathConf
            ]
            ret = subprocess.call(args)
            """
            if ret != 0:
                raise Exception('something went wrong when creating the observationss')
        else:
            shutil.rmtree(xpidobs)
            """
    return list(set(sensors))
