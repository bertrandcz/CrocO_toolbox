# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 24 févr. 2020

@author: cluzetb
create a masked version of the postes observations  to select a specific massif only.
'''

import datetime

from CramponPf import Crampon
from crampon import set_options
from utilcrampon import Pgd
# ####### PARAMETERS #####
vconf = 'postes_8_9_12_13_15_16_csv'
selMassif = [12]
xpidobs = 'all'
# in order to mask the obs that are not in the selected massif,
# we use the classesId argument of crampon.py

# the classesId are in Python indexw (0 for the first class/station.)
pgdPath = '/home/cluzetb/vortexpath/s2m/' + vconf + '/spinup/pgd/super_PGD_' + vconf + '.nc'
pgd = Pgd(pgdPath)

classesId = [i for i, m in enumerate(pgd.massif) if m in selMassif]
print(pgd.station[classesId])
print(classesId)
for i, stat in enumerate(pgd.station[classesId]):
    classesIdX = [cl for j, cl in enumerate(classesId) if j != i]
    classesIdArg = ','.join(list(map(str, classesIdX)))
    sensor = '_'.join(map(str, selMassif)) + '_X' + str(stat)
    args = [

        '/home/cluzetb/Code/Dev/crampon.py',
        '--xpidobs', xpidobs,
        '--xpid', 'testmaskobs',
        '--todo', 'generobs',
        '--vconf', vconf,
        '-d', 'all',
        '--vars', 'DEP',
        '-o', 'testgener',
        '--classesId', classesIdArg,
        '--sensor', sensor,
        '--archive_synth'
    ]
    options, conf = set_options(args)
    run = Crampon(options, conf)
    for date in run.conf.assimdates:
        dt = datetime.datetime.strptime(date, '%Y%m%d%H')
        if dt.month > 9 or dt.month < 7:
            run.prepare_obs(date)
