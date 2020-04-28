# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 6 nov. 2019

@author: cluzetb

Perform local runs and freely postprocess it.
example data : inside the data archive of the summitted manuscript to GMD (Cluzet et al., 2020)
'''
import time

from CrocoPf import CrocoPf
from CrocoPp import CrocoPp
from consts import CROCO
from crocO import set_options
import matplotlib.pyplot as plt
from plotcrocO import plot_part_cesar_from_run


start_time = time.time()
args = [
    CROCO + '/crocO.py',
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

run = CrocoPf(options)  # full loading, expensive but necessary in exploration mode.
run.run()
pp = CrocoPp(options)
for cl in range(85, 87):
    print(pp.ensBg[options.dates[0]].stack['DEP'][85, :])
    print(pp.ensAn[options.dates[0]].stack['DEP'][85, :])
    plot_part_cesar_from_run(pp, cl, 'DEP', options.dates[0], kindobs = 'Arch')
    plt.show()

elapsed_time = time.time() - start_time
print('runtime:', elapsed_time)
