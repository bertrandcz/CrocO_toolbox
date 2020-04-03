# -*- coding: utf-8 -*-
'''
Created on 29 oct. 2019
Copied from pp_SWECRPSfrom_pickle_time.py.
Purpose: laod numerous expireiences for exploring purposes.
@author: cluzetb
'''


from CrocOpp import CrocOpp
from crocO import set_options
import sys
import time


start_time = time.time()
# general params
print(len(sys.argv))
if len(sys.argv) > 1:
    pp_kind = sys.argv[1]
else:
    pp_kind = ''
    print('loading kind of xp : ', pp_kind)
# aggregated scores.
aggrClasses = False
aggrTime = True
score = 'CRPS'

# inflation case

if pp_kind == 'inflation_effect':
    years = [2014, ]
    dictmembers = {2014: [122]}
    nens = 40
    Neff = ['1', '10', '7']
    assimvars = 'DEP'
    ppvars = 'SWE'

    runs = ['global_40_DEP',
            'global_40_DEP_1',
            ]
    readprep = ''
    readobs = ''
elif pp_kind == 'plot_bg':
    years = [2014, ]
    dictmembers = {2014: [122]}
    nens = 40
    Neff = ['7']
    assimvars = ['B4,B5', 'DEP', ]
    ppvars = 'SWE'

    runs = [            'klocal5_40',
                        'klocal5_40_DEP',
                        ]
    readprep = '--readprep'
    readobs = ''
elif pp_kind == 'postesdebug':
    years = [2013, ]
    dictmembers = {2013: [66]}
    nens = 40
    Neff = ['7']
    assimvars  = 'DEP'
    ppvars = 'SWE'
    runs = ['postes_rlocal_40_DEP_7']
    readprep = '--readprep'
    readobs = '--readobs'
else:
    years = [2013]
    dictmembers = {2013: [  # 66,
        # 12,
        122,
        # 149
    ],
        2014: [69, 28, 2, 122],
        2015: [92, 97, 14, 141],
        2016: [50, 153, 90, 117]}
    nens = 40
    Neff = ['1', '10', '7']
    assimvars = 'DEP'
    ppvars = 'SWE'
    suffix = '' if nens == 160 else '_{0}'.format(nens) if 'DEP' not in assimvars else '_{0}_{1}'.format(nens, assimvars)

    runs = ['global',
            # 'rlocal', 'klocal'
            ]
    if len(Neff) > 0:
        suffixs = [suffix + '_' + n if n is not '7' else suffix for n in Neff]
        suffix = suffixs
        SUFFIX_NAME = '_40_DEP_neffstudy'

    runs = [r + s for r in runs for s in suffix]
    readprep = '--readprep'
    readobs = '--readobs'

print(runs)
RUN = dict()
for year in years:
    for mbsynth in dictmembers[year]:
        for ik, key in enumerate(sorted(runs)):
            xp = '{0}_{1}_{2}'.format(year, mbsynth, key)
            print('------loading xp ', xp, '------')
            args = [
                '/home/cluzetb/snowtools_git/assim/crocO.py',
                '--xpid', xp,
                # '--xpid', 'art2_OL_{0}_t1500'.format(year),
                '--xpidol', 'art2_OL_{0}_t1500'.format(year),
                '-d', 'all',
                '--vars', assimvars[ik] if pp_kind == 'plot_bg' else assimvars,
                '--ppvars', ppvars,
                '-o', 'compar2to3',
                '--classesE', '1800,2100,2400,2700,3000,3300,3600' if pp_kind != 'postesdebug' else '1200,1500,1800,2100,2400',
                '--classesS', '0,20' if pp_kind != 'postesdebug' else '0',
                '--classesA', 'SW,S,SE,E',
                readprep,
                readobs,
            ]
            options, conf = set_options(args)
            RUN[xp] = CrocOpp(options, conf)
