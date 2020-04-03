# -*- coding: utf-8 -*-
'''
Created on 5 févr. 2019

@author: cluzetb, inspired on SodaXP, test_PF.py and snowtools_git/tasks/runs.py from Lafaysse
'''
from ParticleFilter import ParticleFilter
from PostCroco import PostCroco
from SemiDistributed import Synthetic, Real
import os
import random
import shutil
from utilcrocO import convertdate, area, check_namelist_soda

import matplotlib.pyplot as plt


# import numpy as np
class CrocO(object):
    '''
    Class for local soda test
    '''

    def __init__(self, options, conf):

        self.options = options
        self.conf = conf

        self.rootdir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/'
        self.xpiddir = options.xpiddir
        if not os.path.exists(self.xpiddir):
            raise Exception('experiment ' + options.xpid  + 'does not exist at ' + self.xpiddir)
        if not hasattr(self.conf, 'openloop'):
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

        # set dirs
        self.crocodir = self.xpiddir + 'crocO/'
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
            # BC 24/02 dirty dirty
            gg = self.options.vortexpath + '/s2m/' + self.options.vconf + '/spinup/pgd/super_PGD_' + area(self.options.vconf) + '.nc'
            self.obs = Real(self.options.xpidobsdir, date, self.options,
                            pgdPath=gg)
            self.obs.prepare(archive_synth = self.options.archive_synth, need_masking = self.options.need_masking)


class CrocOrun(CrocO):
    '''
    class meant to perform LOCAL runs of the pf and post-process it
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

        if self.options.synth is None:
            self.mblist = list(range(1, self.options.nmembers + 1))
        # setup all dirs
        if setup is True:
            self.setup()

    def setup(self):
        if not os.path.exists(self.crocodir):
            os.mkdir(self.crocodir)
        os.chdir(self.crocodir)
        if not os.path.exists(self.options.saverep):
            os.mkdir(self.options.saverep)
        os.chdir(self.options.saverep)
        for dd in self.options.dates:
            if dd in self.conf.assimdates:
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
            os.symlink(self.options.vortexpath + '/s2m/' + self.options.vconf + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')

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
#         try:
#             os.symlink(self.xpiddir + 'mb{0:04d}'.format(mb)  + '/an/PREP_' + path + '.nc', 'SURFOUT' + str(imb) + '.nc')
#         except:
#             if os.path.exists('SURFOUT' + str(imb + 1) + '.nc'):
#                 os.remove('SURFOUT' + str(imb + 1) + '.nc')
#             else:
#                 pass

    def run(self):
        """spawn soda in each date repertory"""
        os.system('ulimit -s unlimited')
        for dd in self.options.dates:
            os.chdir(dd)
            if self.options.todo != 'pfpython':
                os.system('./soda.exe')
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

    def post_proc(self, options):
        if self.options.todo == 'pfpython':
            postp = PostCroco(self.xpiddir, self.xpiddir, options, pf = self.pf)
            pb = postp.run()
            return pb
        else:
            postp = PostCroco(self.xpiddir, self.xpiddir, options)
            postp.run()
