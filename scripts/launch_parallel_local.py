'''
Created on 27 mars 2020

@author: cluzetb

script to test local parallelization of CRAMPON sequence

'''
from CrocOparallel import CrocOparallel
from CrocOpp import CrocOpp
from crocO import set_options
import os
import time
if __name__ == '__main__':
    # #### PARAMETERS #####
    vconf = 'postes_8_9_12_13_15_16_csv'
    sensor = '12'
    pf = 'global'
    neff = 2
    nloc_pf = 35
    assimvars = 'DEP'
    nens = 7
    year = 2014
    assimilate_every = 7
    #######################

    args = [
        '/home/cluzetb/Code/Dev/crocO.py',
        '-d', 'all',
        '-n', os.environ['THESE'] + '/article3/OPTIONS_MOTHER_NO_TARTES.nam',
        '-b', '{0}080106'.format(year),
        '-e', '{0}110106'.format(year),
        '-f', os.environ['VORTEXPATH'] + '/safran/' + vconf + '/forcing_{0}{1}B_31D_11_t1500_160/'.format(year, year + 1),
        '--pf', pf,
        '--vconf', vconf,
        '--sensor', sensor,
        '--neff', str(neff),
        '--nloc_pf', str(nloc_pf),
        '--nmembers', '2',
        '--vars', assimvars,
        '--escroc', 'E1notartes',
        '--xpid', 'test0',
        '--xpidobs', '12',
        '--todo', 'parallel',
        '--spinup', os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/spinup/',
        '--arch', '/home/cluzetb/Documents/these/vortexpath/s2m/postes_8_9_12_13_15_16_csv/arch_test1',
        '-o', 'local2',
        # pp args
        '--ppvars', 'DEP',
        '--readaux', '--readprep',
    ]
    defaultConf = os.environ['THESE'] + '/dev_local_parallel/conf.ini'
    options, conf = set_options(args,
                                pathConf=defaultConf, useVortex=False)
    run = CrocOparallel(options, conf)
    run.run(cleanup=False)
    run.archive()

    # post-processing.
    # by construction, pp is operating on the archive.
    # for this reason, the conf file must be copied to the archive and parser to feed CrocOpp.
    optionspp, confpp = set_options(args, pathConf = options.arch + '/conf/conf.ini', useVortex=False)
    optionspp.xpiddir = options.arch
    start_time = time.time()

    pp = CrocOpp(optionspp, confpp)
    elapsed_time = time.time() - start_time
    print('elapsed time (pp):', elapsed_time)
