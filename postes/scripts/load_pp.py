# -*- coding: utf-8 -*-
'''
Created on 10 févr. 2020

@author: cluzetb
'''
import datetime

from CramponPp import CramponPp
from crampon import set_options
import matplotlib.pyplot as plt
from postes.explore_metadata import find_name_station

# ########## PARAMS ##########################
run = 'global'  # 'ol' for openloop
nens = 40
neff = 7
assimvars = 'DEP'
year = 2014
assimilate_every = 7
vconf = 'postes_8_9_12_13_15_16_csv'
selMassif = [12]

runs = dict()
for sensor in ['12', '12_X38191400']:
    args = [
        '/home/cluzetb/snowtools_git/assim/crampon.py',
        '--xpid', '{0}_{1}_{2}_{3}_{4}_{5}'.format(year, run, nens, sensor, neff, assimilate_every),
        '--xpidol', '2014_ol_40',
        '--sensor', sensor,
        '--vconf', 'postes_8_9_12_13_15_16_csv',
        '-d', 'all',
        '--vars', 'DEP',
        '--ppvars', 'DEP',
        '-o', 'beaufix',
        # '--readprep',
        '--readobs',
        '--readaux'
    ]
    options, conf = set_options(args)

    runs[sensor] = CramponPp(options, conf)
