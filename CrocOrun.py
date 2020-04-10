# -*- coding: utf-8 -*-
'''
Created on 5 févr. 2019

@author: cluzetb, inspired on SodaXP, test_PF.py and snowtools_git/tasks/runs.py from Lafaysse
'''
from ParticleFilter import ParticleFilter
from SemiDistributed import Synthetic, Real
import os
import random
import shutil
import subprocess
from utilcrocO import convertdate, area, check_namelist_soda

import matplotlib.pyplot as plt


# import numpy as np
class CrocO(object):
    '''
    Mother class for local soda test and post-processing
    '''

    def __init__(self, options, conf):

        self.options = options
        self.conf = conf

        self.xpiddir = options.xpiddir
        if not os.path.exists(self.xpiddir):
            if self.options.todo != 'parallel':
                raise Exception('experiment ' + options.xpid  + 'does not exist at ' + self.xpiddir)
            else:
                os.makedirs(self.xpiddir)
        if not hasattr(self.conf, 'openloop'):  # BC to clean
            self.conf.openloop = 'off'
        # set the observations dir
        if self.options.sensor is None and str(self.conf.openloop) == 'off':
            if hasattr(self.conf, 'sensor'):
                self.sensor = self.conf.sensor
                self.options.sensor = self.conf.sensor  # necessary for pp
            else:
                try:
                    mb = 'mb{0:04d}'.format(self.options.synth)
                except TypeError:
                    raise Exception('if you dont specify obs, please specify a synth member to draw')
                self.sensor = mb + '_v' + ''.join(self.options.vars)  \
                    + '_E' + ''.join(self.options.classesE) + '_A' + ''.join(self.options.classesA) \
                    + '_S' + ''.join(self.options.classesS) + '_N' + str(self.options.noise)
                self.options.sensor = self.sensor

        elif self.options.sensor is not None:
            if self.conf.openloop == 'on':
                self.sensor = self.options.sensor
            else:
                if hasattr(self.conf, 'sensor'):
                    self.sensor = self.conf.sensor
                    self.options.sensor = self.conf.sensor  # necessary for pp
                else:
                    self.sensor = self.options.sensor
        else:  # openloop on and options.sensor is None
            print('openloop on and options.sensor is None')
        if self.options.nmembers is None:
            self.options.nmembers = int(self.conf.nmembers)

        self.machine = os.uname()[1]
        if type(self.conf.assimdates) is str:
            self.conf.assimdates = [str(self.conf.assimdates)]
        else:
            self.conf.assimdates = list(map(str, self.conf.assimdates))
        if hasattr(self.conf, 'stopdates'):
            if type(self.conf.stopdates) is str:
                self.conf.stopdates = [str(self.conf.stopdates)]
            else:
                self.conf.stopdates = list(map(str, self.conf.stopdates))
        else:  # parallel case
            self.stopdates = self.conf.assimdates + [self.options.datefin]
        if 'sxcen' not in self.machine:
            self.exesurfex = os.environ['EXESURFEX']
        else:
            self.exesurfex = None

    def prepare_namelist(self):
        """
        Prepare and check the namelist (LWRITE_TOPO must be false for SODA)
        """
        print('copying namelist', os.getcwd())
        if not os.path.exists('OPTIONS.nam'):
            nampathnormal = self.xpiddir + 'conf/OPTIONS_base.nam'
            namelist = nampathnormal if os.path.exists(nampathnormal) else self.xpiddir + 'conf/namelist.surfex.foo'
            shutil.copyfile(namelist, 'OPTIONS_base.nam')

    def prepare_obs(self, date):
        """
        """
        if self.options.synth is not None:
            # synthetic obs is generated from mbsynth at time date
            self.obs = Synthetic(self.xpiddir, date, self.options)
            # first, a (noisy) copy of the member is archived without any mask
            # it will be used for post-processing
            self.obs.prepare(archive_synth = self.options.archive_synth)

        else:
            # real obs are obtained in xpidobs
            # BC 24/02 dirty dirty
            gg = self.options.vortexpath + '/s2m/' + self.options.vconf + '/spinup/pgd/super_PGD_' + area(self.options.vconf) + '.nc'
            self.obs = Real(self.options.xpidobsdir, self.xpiddir, date, self.options,
                            pgdPath=gg)
            self.obs.prepare(archive_synth = self.options.archive_synth, need_masking = self.options.need_masking)


class CrocOpf(CrocO):
    '''
    class meant to perform LOCAL runs of the pf.
    '''

    def __init__(self, options, conf, setup = True):

        self.options = options
        # /!\before calling mother init, in case of synthetical assimilation, we need to remove the synthetical member !!!!
        if self.options.synth is not None:
            self.mblist = list(range(1, self.options.nmembers))
            # draw the synthetical member
            if self.options.synth == 0:
                self.options.synth = random.choice(list(range(1, options.nmembers + 1)))
            # reduce nmembers by 1
            self.options.nmembers -= 1

        # then, call mother init
        CrocO.__init__(self, self.options, conf)
        # set dirs. Croco os the root for the pf xps.
        if self.options.todo == 'parallel':
            self.crocodir = self.xpiddir
        else:
            self.crocodir = self.xpiddir + 'crocO/'
        if self.options.synth is None:
            self.mblist = list(range(1, self.options.nmembers + 1))
        # setup all dirs
        if setup is True:
            self.setup()

    def setup(self):
        if not os.path.exists(self.crocodir):
            os.mkdir(self.crocodir)
        for dd in self.options.dates:
            if dd in self.conf.assimdates:
                print('##############################################')
                print(dd)
                print('##############################################')
                if os.path.exists(self.xpiddir + dd + '/workSODA'):
                    shutil.rmtree(self.xpiddir + dd + '/workSODA')
                os.makedirs(self.xpiddir + dd + '/workSODA')
                self.prepare_sodaenv(dd)
            else:
                raise Exception("prescribed observation dates (-d) must be in conf.assimdates")
#         if hasattr(self.options, 'datefin'):
#             if os.path.exists(self.xpiddir + self.options.datefin + '/workSODA'):
#                 shutil.rmtree(self.xpiddir + self.options.datefin + '/workSODA')
#             os.makedirs(self.xpiddir + self.options.datefin + '/workSODA')
#             self.prepare_sodaenv(self.options.datefin)

    def prepare_sodaenv(self, date):
        """
        set soda environment for each date : -PGD, links to preps, namelist, ecoclimap etc.
        """

        os.chdir(self.xpiddir + '/' + date + '/workSODA')
        self.prepare_namelist()

        # in the parallel case, the namelist is checked by CrocOparallel class
        if self.options.todo != 'parallel':
            check_namelist_soda(self.options)
        else:
            os.rename('OPTIONS_base.nam', 'OPTIONS.nam')
        if 'sxcen' not in self.machine:
            # prepare ecoclimap binaries
            if not os.path.exists('ecoclimapI_covers_param.bin'):
                os.symlink(self.exesurfex + '/../MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin', 'ecoclimapI_covers_param.bin')
                os.symlink(self.exesurfex + '/../MY_RUN/ECOCLIMAP/ecoclimapII_eu_covers_param.bin', 'ecoclimapII_eu_covers_param.bin')
                # flanner stuff
                os.symlink(self.exesurfex + '/../MY_RUN//DATA/CROCUS/drdt_bst_fit_60.nc', 'drdt_bst_fit_60.nc')
            if not os.path.exists('soda.exe'):
                os.symlink(self.exesurfex + '/SODA', 'soda.exe')
            if not os.path.exists('PGD.nc'):
                if self.options.spinup is None:  # in parallel, should not be None
                    os.symlink(self.options.vortexpath + '/s2m/' + self.options.vconf + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')
                else:
                    os.symlink(self.xpiddir + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')
        # prepare (get or fake) the obs
        self.prepare_obs(date)

        # prepare the pgd and prep files.
        # in ocal parallel sequence, must be done only once the corresponding escroc run has produced the files
        if self.options.todo != 'parallel':
            self.prepare_preps(date)
        os.chdir('..')

    def prepare_preps(self, date):
        '''
        Prepare the PREP files for SODA.
        '''
        dateAssSoda = convertdate(date).strftime('%y%m%dH%H')
        for mb in self.mblist:
            if self.options.synth is not None:
                if mb < self.options.synth:
                    self.build_link(mb, mb, date, dateAssSoda)
                else:
                    self.build_link(mb + 1, mb, date, dateAssSoda)
            else:
                self.build_link(mb, mb, date, dateAssSoda)
        # a link fro PREP...1.nc to PREP.nc is also necessary for SODA
        if not os.path.exists('PREP.nc'):
            os.symlink('PREP_' + dateAssSoda + '_PF_ENS1.nc', 'PREP.nc')

    def build_link(self, mb, imb, date, dateAssSoda):
        if self.options.todo != 'parallel':
            try:
                os.symlink(self.xpiddir + 'mb{0:04d}'.format(mb)  + '/bg/PREP_' + date + '.nc', 'PREP_' + dateAssSoda + '_PF_ENS' + str(imb) + '.nc')
            except Exception:
                os.remove('PREP_' + date + '_PF_ENS' + str(imb + 1) + '.nc')
                os.symlink(self.xpiddir + 'mb{0:04d}'.format(mb) + '/bg/PREP_' + date + '.nc', 'PREP_' + dateAssSoda + '_PF_ENS' + str(imb) + '.nc')
        else:  # in parallel mode, bg and an reps do not exist.
            try:
                os.symlink(self.xpiddir + date + '/mb{0:04d}'.format(mb)  + '/SURFOUT.nc', 'PREP_' + dateAssSoda + '_PF_ENS' + str(imb) + '.nc')
            except Exception:
                os.remove('PREP_' + date + '_PF_ENS' + str(imb + 1) + '.nc')
                os.symlink(self.xpiddir + date + '/mb{0:04d}'.format(mb) + '/SURFOUT.nc', 'PREP_' + dateAssSoda + '_PF_ENS' + str(imb) + '.nc')

    def run(self):
        """spawn soda in each date directory"""
        os.system('ulimit -s unlimited')
        for dd in self.options.dates:
            os.chdir(dd)
            if self.options.todo != 'pfpython':
                # os.system('./soda.exe')
                os.system('watch -n 3 "echo toto"')
            else:
                plot = True
                if plot:
                    plt.figure()
                self.pf = ParticleFilter(self.xpiddir, self.options, dd)
                _, resample = self.pf.localize(k=self.options.nloc_pf, errobs_factor=self.options.fact, plot = plot)
                _, gresample = self.pf.globalize(errobs_factor=self.options.fact, plot = plot)
                print(('gres', gresample))
                self.pf.pag = self.pf.analyze(gresample)
                self.pf.pal = self.pf.analyze(resample)
                if plot:
                    plt.show()
            os.chdir('..')

    def run_parallel(self, date):
        """
        This method allows to run soda inside a parallelized assimilation sequence.
        /!\ soda is not parallelized (NOMPI mode).
        """
        os.chdir(self.xpiddir + date + '/workSODA')
        self.prepare_preps(date)
        print('launching SODA on ', date)
        with open('soda.out', 'w') as f:
            p = subprocess.call('./soda.exe', stdout=f)
        os.chdir('..')
