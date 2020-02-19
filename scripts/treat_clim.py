# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 11:14:41 2019
Compute a clim of SWE for grandes-rousses massif.
@author: cluzetb
"""

import calendar
import os
import pickle
from snowtools_git.utils.prosimu import prosimu

import numpy as np
from utilcrocO import dictvarsPro


rootpath = '/era40/vortex/s2m/alp_allslopes/reanalysis/pro/'

startY = 1989
lastY = 2018
clim = dict()
dictVars = dictvarsPro()
pathPkl = '/home/cluzetb/vortexpath/s2m/12/clim/crocO/clim.pkl'
listvar = ['DEP', 'SWE']
for var in listvar:
    clim[var] = np.empty((365, 187, len(range(startY, lastY + 1))))
for iy, year in enumerate(range(startY, lastY + 1)):
    print(year)
    data = prosimu(rootpath + 'PRO_' + str(year) + '080106_' + str(year + 1) + '080106.nc')
    if iy == 0:
        print('leap')
        massif_num = data.read('massif_num')
        mask = np.squeeze(np.array(np.where(massif_num == 12)))
    if calendar.isleap(year + 1):
        print('toto')
        time = data.readtime()
        times = [not(t.month == 2 and t.day == 29) for t in time]
        itimes = [i for i, t in enumerate(times) if t]
    for var in listvar:
        datavar = data.read_var(dictVars[var])[:, mask]
        if calendar.isleap(year + 1):
            datavar = datavar[itimes, :]
        clim[var][:, :, iy] = datavar
    data.close()


if os.path.exists(pathPkl):
    os.remove(pathPkl)
with open(pathPkl, 'wb') as f:
    pickle.dump(clim, f)
