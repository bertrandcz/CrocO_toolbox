# -*- coding: utf-8 -*-
'''
Created on 15 oct. 2019

@author: cluzetb
'''

import os

from crampon import set_options
from utilcrampon import ftpconnect, ftp_upload


years = [2013, 2014, 2015, 2016]
dictmembers = {2013: [66, 12, 122, 149],
               2014: [69, 28, 2, 122],
               2015: [92, 97, 14, 141],
               2016: [50, 153, 90, 117]}
dictOl = {2013: 'art2_OL_2013_t1500@cluzetb',
          2014: 'art2_OL_2014_t1500',
          2015: 'art2_OL_2015_t1500',
          2016: 'art2_OL_t1500@cluzetb', }

# kind = 'ol'
kind = 'postes'
# generate observations
for year in years:
    if kind == 'ol' or kind == 'postes':
        llist = dictmembers[year]
    else:
        llist = ['baseline']
    for mbsynth in llist:
        args = ['python3 /home/cluzetb/assim/crampon.py',
                '--xpid', dictOl[year] if (kind == 'ol' or kind == 'postes')
                else 'baseline_{0}'.format(year),
                '--synth', str(mbsynth) if (kind == 'ol' or kind == 'postes') else '1',
                '-d', 'all',
                '--nmembers', '160',
                '--todo', 'generobs',
                '--vars', 'B1,B2,B3,B4,B5,B6,B7,SCF,DEP,SWE',
                '--ppvars', 'B1,B2,B3,B4,B5,B6,B7,SCF,DEP,SWE',
                '-o', 'testpostes',
                '--classesE', '1800,2100,2400,2700,3000,3300,3600' if kind != "postes" else '1200,1500,1800,2100,2400',
                '--classesS', '0,20' if kind != "postes" else '0',
                '--classesA', 'E,SE,S,SW',
                '--noise', '0',
                '--sensor', 'mb{0:04d}'.format(mbsynth) if kind == 'ol'
                else 'mb{0:04d}_postes'.format(mbsynth) if kind == 'postes'
                else 'baseline' if kind == 'baseline' else None,
                '--archive_synth']
        options, conf = set_options(args)
        os.system(' '.join(args))

        # put it on hendrix
        ftp = ftpconnect('hendrix')
        if kind == 'ol':
            os.chdir('{0}crampon/SYNTH/mb{1:04d}/'.format(options.xpiddir, mbsynth))
        elif kind == 'postes':
            os.chdir('{0}crampon/SYNTH/mb{1:04d}_postes/'.format(options.xpiddir, mbsynth))
        elif kind == 'baseline':
            os.chdir('{0}crampon/SYNTH/baseline/'.format(options.xpiddir))
        for file in os.listdir('.'):
            if kind == 'ol':
                remotefile = '/home/cluzetb/vortex/{0}/{1}/obs/mb{2:04d}/{3}'.format(options.vapp, options.vconf, mbsynth, file)
            elif kind == 'postes':
                remotefile = '/home/cluzetb/vortex/{0}/{1}/obs/mb{2:04d}_postes/{3}'.format(options.vapp, options.vconf, mbsynth, file)
            elif kind == 'baseline':
                remotefile = '/home/cluzetb/vortex/{0}/{1}/obs/baseline/{2}'.format(options.vapp, options.vconf, file)
            ftp_upload(file, remotefile, ftp)
        ftp.close()
