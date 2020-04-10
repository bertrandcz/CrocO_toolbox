# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 6 nov. 2019

@author: cluzetb

Perform local runs and freely postprocess it
'''
from CramponPf import CramponPf
from CramponPp import CramponPp
from consts import CRAMPON
from crampon import set_options
from plotcrampon import plot_part_cesar_from_run
import time

import matplotlib.pyplot as plt


start_time = time.time()
args = [
    CRAMPON + '/crampon.py',
    '--xpid', 'art2_OL_2014_t1500',
    '-d', '2014110710',
    '--todo', 'pf',
    '--pf', 'rlocal',
    '--synth', '2',
    '--neff', '1',
    '--nloc_pf', '1',
    '--nmembers', '3',
    '--archive_synth',
    '--vars', 'DEP',
    '--ppvars', 'DEP',
    '-o', 'ex_local_run_pf',
    # '--classes_e', '1800,2100,2400,2700,3000,3300,3600',
    '--classes_s', '0,20',
    # '--classes_a', 'all',
    '--readprep',
    '--readobs',
    '--readaux',
    '--notreadpro',

]
options = set_options(args)

run = CramponPf(options)  # full loading, expensive but necessary in exploration mode.
run.run()
pp = CramponPp(options)
for cl in range(85, 87):
    print(pp.ensBg[options.dates[0]].stack['DEP'][85, :])
    print(pp.ensAn[options.dates[0]].stack['DEP'][85, :])
    plot_part_cesar_from_run(pp, cl, 'DEP', options.dates[0], kindobs = 'Arch')
    plt.show()

elapsed_time = time.time() - start_time
print('runtime:', elapsed_time)
