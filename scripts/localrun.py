# /usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 6 nov. 2019

@author: cluzetb

Perform local runs and freely postprocess it
'''
from CramponPp import CramponPp
from CramponPf import CramponPf
from crampon import set_options
import matplotlib.pyplot as plt
from plotcrampon import plot_part_cesar_from_run
import time
start_time=time.time()
args = [
    '/home/cluzetb/Code/Dev/crampon.py',
    '--xpid', 'artB31D11_2016@cluzetb',
    '-d', '2016111210',
    '--todo', 'pf',
    '--kind', 'localpp',
    '--pf', 'rlocal',
    '--synth', '12',
    '--neff', '5',
    '--nloc_pf', '1',
    '--nmembers', '35',
    '--archive_synth',
    '--vars', 'DEP',
    '--ppvars', 'DEP',
    '-o', 'testDEPk',
    # '--classesE', '1800,2100,2400,2700,3000,3300,3600',
    '--classesS', '0,20',
    # '--classesA', 'all',
    '--readprep',
    '--readobs',
    '--readaux',
    '--notreadpro',

]
options, conf = set_options(args)

run = CramponPf(options, conf)  # full loading, expensive but necessary in exploration mode.
run.run()
"""
pp = CramponPp(options, conf)
for cl in range(85, 120):
    plt.figure()
    plot_part_cesar_from_run(pp, cl, 'DEP', '2016111210', kindobs = 'Arch')
    plt.show()
"""

elapsed_time = time.time() - start_time
print('runtime:', elapsed_time)