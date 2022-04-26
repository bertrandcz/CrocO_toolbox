'''
Created on 27 mars 2020

@author: cluzetb

script to test local parallelization of CROCO sequence

'''
from CrocoParallel import CrocoParallel
from CrocoPp import CrocoPp
from crocO import set_options
import os
import time

if __name__ == '__main__':
    # #### PARAMETERS #####
    vconf = 'postes_8_9_12_13_15_16_csv'
    sensor = '12'
    pf = 'ol'  # openloop mode
    neff = 2
    nloc_pf = 35
    assimvars = 'DEP'
    nens = 7
    year = 2014
    archtest = os.environ['CROCOPATH'] + '/s2m/' + vconf + '/arch_test_ol/'
    #######################

    args = [
        '/home/cluzetb/Documents/these/assim/crocO.py',
        '-d', 'all',
        '-n', os.environ['HOME'] + '/article3/OPTIONS_MOTHER_NO_TARTES.nam',
        '-b', '{0}080106'.format(year),
        '-e', '{0}110106'.format(year),
        '-f', os.environ['VORTEXPATH'] + '/safran/' + vconf + '/forcing_{0}{1}B_31D_11_t1500_160/'.format(year, year + 1),
        '--pf', pf,
        '--vconf', vconf,
        '--sensor', sensor,
        '--neff', str(neff),
        '--nmembers', '2',
        '--vars', assimvars,
        '--escroc', 'E1notartes',
        '--xpid', 'testOl',
        '--todo', 'parallel',
        '--spinup', os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/spinup/',
        '--arch', archtest,
        '-o', 'localpp1304',
        # pp args
        '--ppvars', 'DEP',
        '--readaux', '--readprep',
        '--no_need_masking',
    ]
    defaultConf = os.environ['HOME'] + '/article3/conf.ini'
    options = set_options(args,
                          pathConf=defaultConf, useVortex=False)
    run = CrocoParallel(options)
    run.run(cleanup=False)
    # BC 11/20 deactivate post-processing.
    if 0:
        # post-processing.
        # by construction, pp is operating on the archive.
        # for this reason, the conf file must be copied to the archive and parser to feed CrocoPp.
        args = [
            '/home/cluzetb/Documents/these/assim/crocO.py',
            '-d', 'all',
            '-n', os.environ['HOME'] + '/article3/OPTIONS_MOTHER_NO_TARTES.nam',
            '-b', '{0}080106'.format(year),
            '-e', '{0}110106'.format(year),
            '-f', os.environ['VORTEXPATH'] + '/safran/' + vconf + '/forcing_{0}{1}B_31D_11_t1500_160/'.format(year, year + 1),
            '--pf', pf,
            '--vconf', vconf,
            '--sensor', sensor,
            '--neff', str(neff),
            '--nmembers', '2',
            '--vars', assimvars,
            '--escroc', 'E1notartes',
            '--xpid', archtest,
            '--xpidobs', '12',
            '--todo', 'parallel',
            '--spinup', os.environ['VORTEXPATH'] + '/s2m/' + vconf + '/spinup/',
            '-o', 'localpp1304',
            # pp args
            '--ppvars', 'DEP',
            '--readaux', '--readprep',
            '--no_need_masking',
        ]

        optionspp = set_options(args, pathConf = archtest + '/conf/conf.ini', useVortex=False, mutable = True)
        start_time = time.time()

        pp = CrocoPp(optionspp)
        elapsed_time = time.time() - start_time
        print('elapsed time (pp):', elapsed_time)
