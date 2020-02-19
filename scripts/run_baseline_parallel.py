'''
Created on 24 oct. 2019
Run the baseline experiment locally for the 4 considered years.
@author: cluzetb
'''
import os
import shutil
import subprocess
from threading import Thread

from utilcrocO import read_conf


years = [2013]  # 2014, 2015, 2016]
programs = []
# os.chdir('{0}/s2m/12/baseline/'.format(os.environ['VORTEXPATH']))
os.chdir('/manto/cluzetb/testRafife/baseline_test_rafife')


def baseline(year):
    conf = read_conf('s2m_12_{0}.ini'.format(year))
    begindates = conf.assimdates.copy()
    enddates = conf.assimdates.copy()
    begindates.insert(0, '{0}0801'.format(year))
    enddates.append('{0}0730'.format(year + 1))
    outdir = 'baseline_{0}'.format(year)
    if os.path.exists('baseline_{0}'.format(year)):
        shutil.rmtree(outdir)
    os.makedirs('{0}/prep/'.format(outdir))
    shutil.copyfile('PREP_2016080106.nc', '{0}/prep/PREP_2016080106.nc'.format(outdir))
    for i, (b, e) in enumerate(zip(begindates, enddates)):
        args = ['/home/cluzetb/snowtools_git/tasks/s2m_command.py',
                '-b', b,
                '-e', e,
                '-m', 'safran',
                '-r', '12',
                '-f', '{0}/safran/12/forcing_{1}{2}B_31D_11_t1500_baseline/mb0001/meteo/FORCING_{1}080106_{2}080106.nc'.format(os.environ['VORTEXPATH'], year, year + 1),
                '-n', 'OPTIONS_crampon_baseline.nam'.format(os.environ['VORTEXPATH']),
                '-o', 'baseline_{0}'.format(year),
                # '-s', '/manto/cluzetb/testRafife',
                ]
        if i == 0:
            args.append('-x')
            args.append('2016080106')
        print(args)
        p = subprocess.Popen(args)
        _, _ = p.communicate()
        if p.returncode is not 0:
            raise Exception


for year in years:
    t = Thread(target = baseline, args = (year,))
    t.start()
