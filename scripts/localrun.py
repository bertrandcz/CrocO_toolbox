'''
Created on 6 nov. 2019

@author: cluzetb

Perform local runs and freely postprocess it
'''
from CrocOrun import CrocOrun
from crocO import set_options
args = [
    '/home/cluzetb/Code/Dev/crocO.py',
    '--xpid', 'artB31D20_2016_150aleeee@cluzetb',
    '-d', '2016111210',
    '--todo', 'pf',
    '--kind', 'localpp',
    '--pf', 'global',
    '--synth', '12',
    '--neff', '20',
    '--nloc_pf', '0',
    '--nmembers', '100',
    '--archive_synth',
    '--vars', 'DEP',
    '--ppvars', 'DEP,SWE',
    '-o', 'testDEPk',
    # '--classesE', '1800,2100,2400,2700,3000,3300,3600',
    '--classesS', '0,20',
    # '--classesA', 'SW,S,SE,E',
]
options, conf = set_options(args)

run = CrocOrun(options, conf)  # full loading, expensive but necessary in exploration mode.
run.run()
