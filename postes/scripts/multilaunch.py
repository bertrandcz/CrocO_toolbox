# -*- coding: utf-8 -*-
'''
Created on 15 oct. 2019

@author: cluzetb
'''
import os
import re

from crocO import set_options
from postes.utilpostes import set_conf_everydate
from tasks.vortex_kitchen import vortex_conf_file
from utilcrocO import Pgd
from utilcrocO import check_namelist_soda
from utilcrocO import get_trailing_number
# ########## PARAMS ##########################
pruns = ['noX_global']  # 'ol' for openloop, prefix with 'noX for assim run with the whole massif (without exclusion of posts)
pnens = 40
pneff = 7
passimvars = 'DEP'
pyears = [2013, 2014, 2015, 2016]
passimilate_every = 7
pvconf = 'postes_8_9_12_13_15_16_csv'
pselMassif = [12]
psensorBase = '_'.join(list(map(str, pselMassif)))
pfact = {r: '1' if 'k' not in r else (get_trailing_number(r)) for r in pruns}


def spawn_crampon(year, run, fact, sensor, nens, neff, assimvars, assimilate_every, vconf):
    if run == 'ol':
        xp = '{0}_{1}_{2}'.format(year, run, nens)
        pf = 'rlocal'  # will not be used because Ol
    else:
        if len(run.split('_')) > 1:  # check if noX
            pf = run.split('_')[1]
        else:
            pf = run
        xp = '{0}_{1}_{2}_{3}_{4}_{5}'.format(year, pf, nens, sensor, neff, assimilate_every)
    args = [
        '/home/cluzetb/Code/Dev/crocO.py',
        '-d', 'all',
        '--pf', pf if pf is not 'ol' else 'rlocal',
        '--vconf', vconf,
        '--sensor', sensor,
        '--neff', str(neff),
        '--nloc_pf', '35' if 'klocal' in run else '1' if ('rlocal' in run or 'ol' in run) else '0',
        '--nmembers', str(nens),
        '--vars', assimvars,
        '--fact', fact[run],
        '--ppvars', 'DEP,SWE',
    ]
    defaultConf = '/home/cluzetb/article3/{0}/s2m_{1}_{2}_{3}.ini'.format(nens, vconf, year, assimilate_every)
    ok1 = set_conf_everydate(year, assimilate_every, defaultConf, nens = nens)
    options, conf = set_options(args,
                                pathConf=defaultConf)
    pathNamRun = '/home/cluzetb/article3/tmp/OPTIONS_{0}.nam'.format(xp)
    check_namelist_soda(options, pathIn = '/home/cluzetb/article3/OPTIONS_MOTHER.nam', pathOut = pathNamRun)
    print('assim dates: ', options.dates)
    if 'all' in options.dates:
        pathConfRun = defaultConf
    else:
        pathConfRun = '/home/cluzetb/article3/tmp/{0}_{1}_{2}.ini'.format(options.vapp, options.vconf, xp)
        if os.path.exists(pathConfRun):
            print('remove vortex conf_file')
            os.remove(pathConfRun)
        confRun = vortex_conf_file(pathConfRun, mode='w')
        confRun.new_class('DEFAULT')
        confRun.write_field('assimdates', ','.join(options.dates))
        confRun.write_field('membersId', ','.join(map(str, conf.membersId)))
        confRun.close()

    # 2. launch the run
    cmd = ['python /home/cluzetb/snowtools_git/tasks/s2m_command.py',
           '-b', '{0}0801'.format(year),
           '-e', '{0}0730'.format(year + 1),
           '-f', 'forcing_{0}{1}B_31D_11_t1500_160'.format(year, year + 1),
           '-m', 'safran',
           '--escroc', 'E1tartes',
           '--crampon', pathConfRun,
           '-n', pathNamRun,
           '-o', xp,
           '-r', vconf,
           '--nforcing', '{0}'.format(nens),
           '--nmembers', '{0}'.format(nens),
           '--nnodes', '{0}'.format(int(nens / 40)),
           '--walltime', '200',
           '--writesx',
           '--pickleit',
           '-x', '2016080106',
           '--sensor ' + sensor if run != 'ol' else '',
           '--openloop' if run == 'ol' else '--real'
           ]
    print('s2m_command :', ' '.join(cmd))
    os.system(' '.join(cmd))
    return 0
# ############################################


pgdPath = '/home/cluzetb/vortexpath/s2m/' + pvconf + '/spinup/pgd/super_PGD_' + pvconf + '.nc'
pgd = Pgd(pgdPath)

# generate runs on beaufix
for year in pyears:
    for run in pruns:
        if run == 'ol':
            sensor = 'osef'
            ok = spawn_crampon(year, run, pfact, sensor, pnens, pneff, passimvars, passimilate_every, pvconf)
        # if required, spawn also the assim without excluding any post
        elif run.startswith('noX'):
            sensor = psensorBase
            ok = spawn_crampon(year, run, pfact, sensor, pnens, pneff, passimvars, passimilate_every, pvconf)
        else:
            classesId = [i for i, m in enumerate(pgd.massif) if m in pselMassif]
            for i, stat in enumerate(pgd.station[classesId]):
                classesIdX = [cl for j, cl in enumerate(classesId) if j != i]
                classesIdArg = ','.join(list(map(str, classesIdX)))
                sensor = psensorBase + '_X' + str(stat)
                # 1.  prepare namelist and conf files.
                ok = spawn_crampon(year, run, pfact, sensor, pnens, pneff, passimvars, passimilate_every, pvconf)
