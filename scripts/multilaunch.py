# -*- coding: utf-8 -*-
'''
Created on 15 oct. 2019

@author: cluzetb
'''
import os
import re

from crocO import set_options
from tasks.vortex_kitchen import vortex_conf_file
from utilcrocO import check_namelist_soda


# ########## PARAMS ##########################
runs = ['global', 'klocal1', 'rlocal']
nens = 40
neff = 7
# kind = 'fromol'
kind = 'postes'
assimvars = 'DEP'
years = [2013, 2014, 2015, 2016]


def get_trailing_number(s):
    m = re.search(r'\d+$', s)
    return str(int(m.group())) if m else None


def get_leading_number(s):
    m = re.search(r'^\d+', s)
    return str(int(m.group())) if m else None


fact = {r: '1' if 'k' not in r else (get_trailing_number(r)) for r in runs}
# ############################################

dictmembers = {2013: [66,
                      12,
                      122,
                      149
                      ],
               2014: [69,
                      28,
                      2,
                      122
                      ],
               2015: [92,
                      97,
                      14,
                      141
                      ],
               2016: [50,
                      153,
                      90,
                      117
                      ]
               }
dictOl = {2013: 'art2_OL_2013_t1500@cluzetb',
          2014: 'art2_OL_2014_t1500',
          2015: 'art2_OL_2015_t1500',
          2016: 'art2_OL_t1500@cluzetb', }
# dates = dict()
# dates['DEP'] = dict()
# dates['DEP'][2013] = '2013112210,2013122910,2014012510,2014022211,2014032810,2014042711'
# dates['DEP'][2016] = '2016111210,2016121410,2017011511,2017021611,2017031310,2017041410,2017051610,2017060811'
# dates['B4,B5'] = dict()
# dates['B4,B5'][2013] = 'all'
# dates['B4,B5'][2016] = 'all'
suffixrun = '' if nens == 160 else '_{0}'.format(nens) if 'DEP' not in assimvars else '_{0}_{1}'.format(nens, assimvars)
suffixrun = '_40_DEP_' + str(neff)
# generate runs on beaufix
for year in years:
    if kind == 'fromol':
        llist = dictmembers[year]
    elif kind == 'postes':
        llist = ['{0}_{1}'.format(mb, kind) for mb in dictmembers[year]]

    elif kind == 'baseline':
        llist = ['baseline']
    for mbsynth in llist:
        if kind == 'postes':
            mbnum = get_leading_number(mbsynth)
        for run in runs:

            # 1.  prepare namelist and conf files.
            xp = '{0}_{1}_{2}{3}'.format(year, mbsynth, run, suffixrun)
            args = [
                '/home/cluzetb/Code/Dev/crocO.py',
                '-d', 'all',
                '--pf', run,
                '--synth', str(mbsynth) if kind != 'postes' else mbnum,
                '--neff', str(neff),
                '--nloc_pf', '35' if 'klocal' in run else '1' if 'rlocal' in run else '0',
                '--nmembers', str(nens),
                '--vars', assimvars,
                '--fact', fact[run],
                '--ppvars', 'DEP,SWE',
            ]
            defaultConf = '/home/cluzetb/article2/{0}/s2m_12_{1}.ini'.format(nens, year)
            options, conf = set_options(args,
                                        pathConf=defaultConf)
            pathNamRun = '/home/cluzetb/article2/tmp/OPTIONS_{0}.nam'.format(xp)
            check_namelist_soda(options, pathIn = '/home/cluzetb/ggg.nam', pathOut = pathNamRun)
            print(options.dates)
            if 'all' in options.dates:
                pathConfRun = defaultConf
            else:
                pathConfRun = '/home/cluzetb/article2/tmp/{0}_{1}_{2}.ini'.format(options.vapp, options.vconf, xp)
                if os.path.exists(pathConfRun):
                    print('remove vortex conf_file')
                    os.remove(pathConfRun)
                confRun = vortex_conf_file(pathConfRun, mode='w')
                confRun.new_class('DEFAULT')
                confRun.write_field('assimdates', ','.join(options.dates))
                confRun.write_field('membersId', ','.join(map(str, conf.membersId)))
                confRun.close()

            # 2.  launch the run
            cmd = ['python /home/cluzetb/snowtools_git/tasks/s2m_command.py',
                   '-b', '{0}0801'.format(year),
                   '-e', '{0}0730'.format(year + 1),
                   '-f', 'forcing_{0}{1}B_31D_11_t1500_160'.format(year, year + 1),
                   '-m', 'safran',
                   '-r', '12',
                   '--escroc', 'E1tartes',
                   '--crampon', pathConfRun,
                   '-n', pathNamRun,
                   '-o', xp,
                   '--nforcing', '{0}'.format(nens),
                   '--nmembers', '{0}'.format(nens),
                   '--nnodes', '{0}'.format(int(nens / 40)),
                   '--walltime', '200',
                   '--writesx',
                   '--pickleit'
                   '-x', '2016080106']
            print('cmd', ' '.join(cmd))
            if kind == 'fromol':
                cmd.append('--synth {0}'.format(mbsynth))
                cmd.append('--sensor ' + 'mb{0:04d}'.format(mbsynth))
            if kind == 'postes':
                cmd.append('--synth {0}'.format(mbnum))
                cmd.append('--sensor ' + 'mb{0:04d}_{1}'.format(int(mbnum), kind))

            else:
                cmd.append('--real')
                cmd.append('--sensor baseline')

            os.system(' '.join(cmd))
