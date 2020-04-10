# -*- coding: utf-8 -*-
'''
Created on 5 févr. 2019

@author: cluzetb, inspired on SodaXP, test_PF.py and snowtools_git/tasks/runs.py from Lafaysse
'''
from ParticleFilter import ParticleFilter
from SemiDistributed import Synthetic, Real
import os
import shutil
from utilcrampon import convertdate, area, check_namelist_soda

import matplotlib.pyplot as plt


# import numpy as np
class Crampon(object):
    '''
    Mother class for CRAMPON pf local tests, local CRAMPON runs and post-processing
    '''

    def __init__(self, options):

        self.options = options

        self.rootdir = options.cramponpath + '/' + options.vapp + '/' + options.vconf + '/'
        self.xpiddir = options.xpiddir

        if not os.path.exists(self.xpiddir):
            raise Exception('experiment ' + options.xpid  + 'does not exist at ' + self.xpiddir)

        # set dirs
        self.crampondir = self.xpiddir + 'crampon/'
        self.machine = os.uname()[1]

        if 'sxcen' not in self.machine:
            self.exesurfex = os.environ['EXESURFEX']
        else:
            self.exesurfex = None

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
            self.obs = Real(self.options.xpidobsdir, date, self.options)
            self.obs.prepare(archive_synth = self.options.archive_synth, need_masking = self.options.need_masking)


class CramponPf(Crampon):
    '''
    class meant to perform LOCAL runs of the pf
    '''

    def __init__(self, options, setup = True):

        Crampon.__init__(self, options)
        # for the local soda PF, it is safer if mblist is a continous range (no removal of the synth mb but one mb less.
        # handling of the synth member is done in prepare_soda_env
        self.mblist = list(range(1, options.nmembers + 1))
        # setup all dirs
        if setup is True:
            self.setup()

    def setup(self):
        if not os.path.exists(self.crampondir):
            os.mkdir(self.crampondir)
        os.chdir(self.crampondir)
        if not os.path.exists(self.options.saverep):
            os.mkdir(self.options.saverep)
        os.chdir(self.options.saverep)
        for dd in self.options.dates:
            if dd in self.options.dates:
                if os.path.exists(dd):
                    # pass
                    shutil.rmtree(dd)
                os.mkdir(dd)
                self.prepare_sodaenv(dd)
            else:
                print(('prescribed date ' + dd + 'does not exist in the experiment, remove it.'))
                self.options.dates.remove(dd)

    def prepare_sodaenv(self, path):
        """
        set soda environment for each date (=path): -PGD, links to preps, namelist, ecoclimap etc.
        """

        os.chdir(path)
        # Prepare the PGD and PREP for assim
        dateAssSoda = convertdate(path).strftime('%y%m%dH%H')
        for mb in self.mblist:
            if self.options.synth is not None:
                if mb < self.options.synth:
                    self.build_link(mb, mb, path, dateAssSoda)
                else:
                    self.build_link(mb + 1, mb, path, dateAssSoda)
            else:
                self.build_link(mb, mb, path, dateAssSoda)
        if not os.path.exists('PREP.nc'):
            os.symlink('PREP_' + dateAssSoda + '_PF_ENS1.nc', 'PREP.nc')

        if not os.path.exists('PGD.nc'):
            os.symlink(self.options.pathPgd, 'PGD.nc')

        # Prepare and check the namelist (LWRITE_TOPO must be false for SODA)
        if not os.path.exists('OPTIONS.nam'):
            shutil.copyfile(self.xpiddir + 'conf/namelist.surfex.foo', 'OPTIONS_base.nam')
        check_namelist_soda(self.options)

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
        self.prepare_obs(path)
        os.chdir('..')

    def build_link(self, mb, imb, path, dateAssSoda):
        try:
            os.symlink(self.xpiddir + 'mb{0:04d}'.format(mb)  + '/bg/PREP_' + path + '.nc', 'PREP_' + dateAssSoda + '_PF_ENS' + str(imb) + '.nc')
        except Exception:
            os.remove('PREP_' + path + '_PF_ENS' + str(imb + 1) + '.nc')
            os.symlink(self.xpiddir + 'mb{0:04d}'.format(mb) + '/bg/PREP_' + path + '.nc', 'PREP_' + dateAssSoda + '_PF_ENS' + str(imb) + '.nc')

    def run(self):
        """spawn soda in each date repertory"""
        os.system('ulimit -s unlimited')
        for dd in self.options.dates:
            os.chdir(dd)
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
                _, resample = self.pf.localize(k=self.options.nloc_pf, errobs_factor=self.options.fact, plot = plot)
                _, gresample = self.pf.globalize(errobs_factor=self.options.fact, plot = plot)
                print(('gres', gresample))
                self.pf.pag = self.pf.analyze(gresample)
                self.pf.pal = self.pf.analyze(resample)
                if plot:
                    plt.show()
            os.chdir('..')


class CramponObs(Crampon):

    def __init__(self, options):
        Crampon.__init__(self, options)
        print('ooooobs',)
        self.setup()

    def setup(self):
        if not os.path.exists(self.crampondir):
            os.mkdir(self.crampondir)
        os.chdir(self.crampondir)
        if not os.path.exists(self.options.saverep):
            os.mkdir(self.options.saverep)
        os.chdir(self.options.saverep)
        for dd in self.options.dates:
            if dd in self.options.dates:
                if os.path.exists(dd):
                    # pass
                    shutil.rmtree(dd)
                os.mkdir(dd)
                os.chdir(dd)
                self.prepare_obs(dd)
                os.chdir('..')
