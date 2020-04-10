# -*- coding: utf-8 -*-
'''
Created on 29 oct. 2019
Example of script used to load and manipulate many Crampon experiences.
Purpose: load numerous experiences for exploring purposes.
@author: cluzetb
'''


from CramponPp import CramponPp
from consts import CRAMPON
from crampon import set_options
import time
start_time = time.time()
# ############ PARAMETERS ####################
years = [2013]
dictmembers = {2013: [66, 12, 122, 149],
               2014: [69, 28, 2, 122],
               2015: [92, 97, 14, 141],
               2016: [50, 153, 90, 117]}
nens = 40
assimvars = 'DEP'
ppvars = 'SWE'
runs = ['global', 'rlocal', 'klocal']
readprep = '--readprep'
readobs = '--readobs'
##############################################
suffix = '' if nens == 160 else '_{0}'.format(nens) if 'DEP' not in assimvars else '_{0}_{1}'.format(nens, assimvars)


runs = [r + s for r in runs for s in suffix]


print(runs)
RUN = dict()
for year in years:
    for mbsynth in dictmembers[year]:
        for ik, key in enumerate(sorted(runs)):
            xp = '{0}_{1}_{2}'.format(year, mbsynth, key)
            print('------loading xp ', xp, '------')
            args = [
                CRAMPON + '/crampon.py',
                '--xpid', xp,
                '--xpidol', 'art2_OL_{0}_t1500'.format(year),
                '-d', 'all',
                '--vars', assimvars,
                '--ppvars', ppvars,
                '-o', 'gmd',
                '--classes_e', 'all',
                '--classes_s', 'all',
                '--classes_a', 'all',
                readprep,
                readobs,
            ]
            options = set_options(args)
            RUN[xp] = CramponPp(options)
