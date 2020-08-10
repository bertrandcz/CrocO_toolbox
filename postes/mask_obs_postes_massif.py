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
import numpy as np
from postes.utilpostes import set_conf_everydate
from utilcrocO import Pgd


def create_obs_massif_exclude_poste(vconf, selMassif, year, assimilate_every, alwaysX = [], return_dict = False, overwrite = False, sensor_in = None):
    """
    this method allows to generate a set of sensors within a given geometry
    by keeping, for each sensor, all observations but one.
    also generates a sensor with all the observations.
    Details:
    -> assimilate only obs from a given selMassif of this geometry
    -> within this selMassif, you can exclude the alwaysX posts from all the observations.
    -> 'all' : all obs, including posts from alwaysX
    -> 'alX' : all obs, excluding posts from alwaysX
    """
    pgdPath = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/spinup/pgd/PGD_' + vconf + '.nc'
    pgd = Pgd(pgdPath)
    os.environ['CROCOPATH'] = os.environ['VORTEXPATH']
    # BC june 2020 : change of functionality: does not make much sense to keep only obs from a specific massif inside the vconf..
    # selMassif is used to operate the exclusion of posts over a set of massifs only
    # classesId = [i for i, m in enumerate(pgd.massif) if m in selMassif]
    classesId = [i for i, m in enumerate(pgd.massif)]
    classesIdXclude = [i for i, m in enumerate(pgd.massif) if m in selMassif]
    if sensor_in is None:
        sensor_in = '_'.join(map(str, selMassif))
    sensors = []

    # BC 10/05 : ensure the all and alX (if necess.) runs are ran.
    # -1 will stand for the 'all' case
    # alX will be necessarily ran since stations in alX are in pgd.station[classesId]
    list_enum = np.append(pgd.station[classesIdXclude], -1)

    if len(alwaysX) < 1:
        alX = False
    else:
        alX = True
    ret_dict = {}
    for i, stat in enumerate(list_enum):
        # prevent from excluding twice
        # classesIdK is the Ids of stations KEPT for assimilation
        classesIdK = list(set([cl for cl in classesId if pgd.station[cl] != stat and str(pgd.station[cl]) not in alwaysX]))
        classesIdArg = ','.join(list(map(str, classesIdK)))
        if alX:
            if stat != -1:
                listX = list(sorted(set([str(stat)] + alwaysX)))
                statX = '-'.join(listX)
                sensor = sensor_in + '_X' + str(statX)
                if str(stat) in alwaysX:
                    ret_dict['alX'] = sensor
                else:
                    ret_dict[str(stat)] = sensor
            else:
                sensor = sensor_in
                ret_dict['all'] = sensor

        else:
            if stat != -1:
                sensor = sensor_in + '_X' + str(stat)
                ret_dict[str(stat)] = sensor
            else:
                sensor = sensor_in
                ret_dict['all'] = sensor
        sensors.append(sensor)

        pathConf = os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/obs/' + sensor + '/conf/s2m_' + vconf + '_' + str(year) + '.ini'
        if overwrite is True or not os.path.exists(pathConf):
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

    if return_dict:
        return ret_dict
    else:
        return list(set(sensors))
