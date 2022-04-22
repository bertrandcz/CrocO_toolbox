# -*- coding: utf-8 -*-
'''
Created on 5 févr. 2019

@author: cluzetb, inspired on SodaXP, test_PF.py and snowtools_git/tasks/runs.py from Lafaysse
'''
from ParticleFilter import ParticleFilter
from SemiDistributed import Synthetic, Real
import os
import shutil
import subprocess
from utilcrocO import convertdate, check_namelist_soda, safe_create_link

import matplotlib.pyplot as plt


# import numpy as np
class CrocO(object):
    '''
    Mother class for CROCO pf local tests, local CROCO runs and post-processing
    '''

    def __init__(self, options):

        self.options = options

        self.rootdir = options.crocOpath + '/' + options.vapp + '/' + options.vconf + '/'
        self.xpiddir = options.xpiddir

        if not os.path.exists(self.xpiddir):
            if self.options.todo == 'parallel'or (self.options.todo == 'generobs' and self.options.synth is None):
                os.mkdir(self.xpiddir)
            else:
                raise Exception('experiment ' + options.xpid  + 'does not exist at ' + self.xpiddir)
        # set dirs
        if self.options.todo == 'parallel':
            self.crocOdir = self.xpiddir
        else:
            self.crocOdir = self.xpiddir + 'crocO/'
        self.machine = os.uname()[1]

        if 'sxcen' not in self.machine:
            self.exesurfex = os.environ['EXESURFEX']
        else:
            self.exesurfex = None

    def prepare_namelist(self):
        """
        copy the namelist to the local dir
        """
        if os.path.exists('OPTIONS_base.nam'):
            os.remove('OPTIONS_base.nam')
        if self.options.namelist_is_default is True:
            nampathnormal = self.xpiddir + 'conf/OPTIONS_base.nam'
        elif os.path.exists(self.options.namelist):
            nampathnormal = self.options.namelist
        else:
            print('the prescribed namelist ' + self.options.namelist + ' does not exist')

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
            # either we are dealing with real obs or we are masking it.
            self.obs = Real(self.options.sensordir, date, self.options)
            self.obs.prepare(archive_synth = self.options.archive_synth, no_need_masking = self.options.no_need_masking)


class CrocoPf(CrocO):
    '''
    class meant to perform LOCAL runs of the pf
    '''

    def __init__(self, options, setup = True):

        CrocO.__init__(self, options)
        # for the local soda PF, it is safer if mblist is a continous range (no removal of the synth mb but one mb less.
        # handling of the synth member is done in prepare_soda_env
        self.mblist = list(range(1, options.nmembers + 1))
        # setup all dirs
        if setup is True:
            self.setup()

    def setup(self):
        if not os.path.exists(self.crocOdir):
            os.mkdir(self.crocOdir)
        os.chdir(self.crocOdir)
        if self.options.todo != 'parallel':
            saverep = self.options.saverep
            if not os.path.exists(saverep):
                os.mkdir(saverep)
        else:
            saverep = ''

        os.chdir(self.crocOdir + '/' + saverep)
        for dd in self.options.dates:
            path = self.crocOdir + '/' + saverep + '/' + dd + '/workSODA'

            if dd in self.options.dates:
                if os.path.exists(path):
                    # bc slight change for local tests where it is painful to have the rep deleted each time. (pwd error)
                    for dirpath, _, filenames in os.walk(path):
                        # Remove regular files, ignore directories
                        for filename in filenames:
                            os.remove(os.path.join(dirpath, filename))
                else:
                    os.makedirs(path)
                if self.options.pf != 'ol':
                    self.prepare_sodaenv(path, dd)
            else:
                print(('prescribed date ' + dd + 'does not exist in the experiment, remove it.'))
                self.options.dates.remove(dd)

    def prepare_sodaenv(self, path, date):
        """
        set soda environment for each date (=path): -PGD, links to preps, namelist, ecoclimap etc.
        """
        cwd = os.getcwd()
        os.chdir(path)
        if not os.path.exists('PGD.nc'):
            os.symlink(self.options.pathPgd, 'PGD.nc')
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

        # prepare (get or fake) the obs
        self.prepare_obs(date)

        # prepare the pgd and prep files.
        # in local parallel sequence, must be done only once the corresponding escroc run has produced the files
        if self.options.todo != 'parallel':
            self.prepare_preps(date)
        os.chdir(cwd)

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
        # a link from PREP...1.nc to PREP.nc is also necessary for SODA
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
            os.chdir(dd + '/workSODA/')
            if self.options.todo != 'pfpython':
                os.system('./soda.exe')
            else:
                # BC April 2020
                # nice but not maintaned and deprecated class
                # which allowed to test and develop locally python version of the PF
                # before implementing it into fortran.
                plot = True
                if plot:
                    plt.figure()
                self.pf = ParticleFilter(self.xpiddir, self.options, dd)
                _, resample = self.pf.localize(k=self.options.pgd.npts, errobs_factor=self.options.fact, plot = plot)
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
        # print('launching SODA on ', date)
        with open('soda.out', 'w') as f:
            p = subprocess.call('./soda.exe', stdout=f)
            if p != 0:
                raise RuntimeError('SODA crashed, check ', os.getcwd() + f)
        os.chdir('..')


class CrocoObs(CrocO):

    def __init__(self, options):
        CrocO.__init__(self, options)
        self.setup()

    def setup(self):
        if not os.path.exists(self.crocOdir):
            os.mkdir(self.crocOdir)
        os.chdir(self.crocOdir)
        if not os.path.exists(self.options.saverep):
            os.mkdir(self.options.saverep)
        os.chdir(self.options.saverep)
        for dd in self.options.dates:
            if os.path.exists(dd):
                # pass
                shutil.rmtree(dd)
            os.mkdir(dd)
            os.chdir(dd)
            self.prepare_obs(dd)
            os.chdir('..')
