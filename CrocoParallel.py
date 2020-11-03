# -*- coding: utf-8 -*-
'''
Created on 26 mars 2020

@author: cluzetb\

The aim of CrocoParallel class is to provide an alternative to vortex to perform CROCO parallalized runs.
Many of its features are similar to snowtools_git/tasks/crocO* files
'''
import datetime
import glob
import multiprocessing
import os
import shutil
import subprocess
import time
from tools.change_prep import prep_tomodify
from tools.update_namelist import update_surfex_namelist_file
from utils.ESCROCsubensembles import ESCROC_subensembles
from utils.dates import get_list_dates_files

# from tqdm import tqdm

from CrocoPf import CrocO, CrocoPf
from utilcrocO import area, check_namelist_soda, dump_conf


class CrocoParallel(CrocO):
    '''
    Setup of the data assimilation sequence and run
    '''

    def __init__(self, options):
        '''
        Constructor
        '''
        self.start_time = time.time()
        # call mother init
        CrocO.__init__(self, options)
        self.dump_options()
        # duplicate with offline pool... but cannot move to CrocO class yet (due tu the synth case removal)
        self.mblist = list(range(1, self.options.nmembers + 1))
        self.mbdirs = ['mb{0:04d}'.format(mb) + '/' for mb in self.mblist]

        self.get_spinup()

        self.setup()

    def dump_options(self):
        """
        rough dump of :
        - options into conf file in conf dir
        - namelist copy into conf dir
        """
        _ = dump_conf(self.xpiddir + '/conf/s2m_' + self.options.vconf + '.ini', self.options)

        # the namelist is roughly copied and receives
        # common necessary modifications (updateloc ?)
        # a check for SODA
        shutil.copyfile(self.options.namelist, self.xpiddir + '/conf/OPTIONS_base.nam')
        cwd = os.getcwd()
        os.chdir(self.xpiddir + '/conf/')
        update_surfex_namelist_file(
            datetime.datetime.strptime(self.options.datedeb, '%Y%m%d%H'),
            dateend = datetime.datetime.strptime(self.options.datefin, '%Y%m%d%H'),
            namelistfile="OPTIONS_base.nam",
            # the two following arguments enable to save lot of setup time
            # by assuming that the location has been previously set (e.g. by the spinup experiment)
            # and that the simulation ends before the forcing lat time step.
            updateloc = False,
            no_caution = True,
            cselect = self.options.provars
        )
        os.chdir(cwd)
        check_namelist_soda(self.options, self.xpiddir + '/conf/OPTIONS_base.nam', self.xpiddir + '/conf/OPTIONS_base.nam')

    def get_spinup(self):

        if os.path.exists(self.xpiddir + '/spinup'):
            shutil.rmtree(self.xpiddir + '/spinup')
        shutil.copytree(self.options.spinup, self.xpiddir + '/spinup')
        # modify the PREP (datemust correspond
        spinup_prepfile = glob.glob(self.xpiddir + '/spinup/prep/PREP*.nc')[0]
        prep = prep_tomodify(spinup_prepfile)
        prep.change_date(datetime.datetime.strptime(self.options.datedeb, '%Y%m%d%H'))
        prep.close()
        shutil.move(spinup_prepfile, self.xpiddir + '/spinup/prep/PREP.nc')

    def setup(self):

        # setup all the simulation dirs (all dates at once !)
        # (same simulation architecture as on beaufix (a little bit simpler maybe)
        # feed them with the constants
        # prepare the observations and put it into the dirs
        # => CrocoPf is used to prepare all that (it can be fed with a list of dates !!)
        # however, the links to prep must be done just before the soda run itself
        # hence, they are encapsulated inside run_parallel() class method
        self.sodas = CrocoPf(self.options)  # first because prepare the directories.
        self.soda_time = time.time() - self.start_time
        self.escrocs = OfflinePools(self.options)
        self.escroc_time = time.time() - self.start_time - self.soda_time
        self.setup_time = time.time() - self.start_time
        print('setup duration:', self.setup_time)
        print('|         soda:', self.soda_time)
        print('|       escroc:', self.escroc_time)

    def run(self, cleanup=False):
        '''
        run
        '''
        # progress_bar
        # pbar = tqdm(self.options.assimdates)
        # for dd in pbar:
        for dd in self.options.assimdates:
            # pbar.set_description('Propagating until: ' + dd)
            #  - spawn offline
            self.escrocs.run(dd)

            # - spawn escroc (if not openloop)
            if self.options.pf and self.options.pf != 'ol':
                self.sodas.run_parallel(dd)
            # report on the time spent
        # last propagation
        self.escrocs.run(self.options.datefin)
        elapsed_time = time.time() - self.start_time
        print('elapsed time(setup and simu) :', elapsed_time)

        # archive
        self.archive()
        # post-process the outputs.
        # cleanup
        if cleanup is True:
            self.cleanup()

    def archive(self):
        """
        archive simulation outputs. could be parallelized.
        """
        start_time = time.time()
        if self.options.arch is None:
            print("putting the archive in xpid. Not recommended")
            self.options.arch = self.xpid
        print('archiving the outputs to ', self.options.arch)
        if not os.path.exists(self.options.arch):
            os.makedirs(self.options.arch)
        os.chdir(self.options.arch)
        if not os.path.exists('mb0001/bg/'):
            for mbdir in self.mbdirs:
                os.makedirs(mbdir + '/' + 'bg')
                os.makedirs(mbdir + '/' + 'an')
                os.makedirs(mbdir + '/' + 'pro')
            os.mkdir('workSODA/')

        # deal with first date
        date = self.options.stopdates[0]
        [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/SURFOUT.nc',
                         mbdir + '/bg/PREP_' + date + '.nc') for mbdir in self.mbdirs]
        [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/ISBA_PROGNOSTIC.OUT.nc',
                         mbdir + '/pro/PRO_' + self.options.datedeb + '_' + date + '.nc') for mbdir in self.mbdirs]
        # following dates
        for idd, date in enumerate(self.options.stopdates[1:]):
            # trap here !! idate starts at 0 instaed of 1.
            # @TODO :integration tes should check the PREP dates correspond to their names.
            idate = idd + 1
            [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/SURFOUT.nc',
                             mbdir + '/bg/PREP_' + date + '.nc') for mbdir in self.mbdirs]
            if self.options.pf and self.options.pf != 'ol':
                [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/PREP.nc',
                                 mbdir + '/an/PREP_' + self.options.stopdates[idate - 1] + '.nc') for mbdir in self.mbdirs]
            [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/ISBA_PROGNOSTIC.OUT.nc',
                             mbdir + '/pro/PRO_' + self.options.stopdates[idate - 1] + '_' + date + '.nc') for mbdir in self.mbdirs]
        if self.options.pf and self.options.pf != 'ol':
            for date in self.options.stopdates[0:-1]:
                shutil.copyfile('/'.join([self.xpiddir, date, 'workSODA']) + '/PART', 'workSODA/PART_' + date + '.txt')
                shutil.copyfile('/'.join([self.xpiddir, date, 'workSODA']) + '/ALPHA', 'workSODA/ALPHA_' + date + '.txt')
                # bg_corr and IMASK only exist with the klocal pf.
                if os.path.exists('/'.join([self.xpiddir, date, 'workSODA']) + '/BG_CORR'):
                    shutil.copyfile('/'.join([self.xpiddir, date, 'workSODA']) + '/BG_CORR', 'workSODA/BG_CORR_' + date + '.txt')
                    shutil.copyfile('/'.join([self.xpiddir, date, 'workSODA']) + '/IMASK', 'workSODA/IMASK_' + date + '.txt')
        if self.options.arch != self.xpiddir:
            if os.path.exists('conf/'):
                shutil.rmtree('conf/')
            shutil.copytree(self.xpiddir + 'conf/', 'conf/')
        # @TODO: also archive the observations.

        elapsed_time = time.time() - start_time
        print('elapsed time(archiving only)', elapsed_time)

    def cleanup(self):
        os.chdir(self.xpiddir)
        for dd in self.options.stopdates:
            shutil.rmtree(dd)


class OfflinePools(CrocO):
    '''
    Class meant to prepare the escroc part of the CROCO sequence
    and parallelize the members between each assimilation date
    '''

    def __init__(self, options):
        CrocO.__init__(self, options)

        # B 26/03/20 CAREFUL WITH SYNTHETIC RUNS
        self.mblist = list(range(1, self.options.nmembers + 1))
        self.mbdirs = ['mb{0:04d}'.format(mb) + '/' for mb in self.mblist]

        # setup
        self.setup()
        # setup the workers (i.e 1 "worker" function that will be mapped on a list of list of arguments.

    def setup(self):

        # prepare escroc configurations
        # BC bugfix 03/06/20:
        # self.mblist instead of self.options.members_id
        # (was running with the 40 first members of E1_notartes....)
        self.escroc_confs = ESCROC_subensembles(self.options.escroc, self.options.members_id)
        p = multiprocessing.Pool(min(multiprocessing.cpu_count(), self.options.nmembers * len(self.options.stopdates)))
        p.map(self.mb_prepare, [[date, idate, mbdir, mb]for idate, date in enumerate(self.options.stopdates)
                                for (mb, mbdir) in zip(self.mblist, self.mbdirs)])
        p.close()
        p.join()
        os.chdir(self.xpiddir)

    def prepare_offline_env(self, date, idate, mbdir, mb):
        """
        set offline environment for each date (=path): -PGD, links to preps, namelist, ecoclimap etc.
        """
        if not os.path.exists('/'.join([self.xpiddir, date, mbdir])):
            os.makedirs('/'.join([self.xpiddir, date, mbdir]))
        os.chdir('/'.join([self.xpiddir, date, mbdir]))
        #  get the forcing and prepare it.
        date_begin_forc, date_end_forc, _, _ = get_list_dates_files(
            datetime.datetime.strptime(self.options.datedeb, '%Y%m%d%H'),
            datetime.datetime.strptime(self.options.datefin, '%Y%m%d%H'),
            'yearly')  # return lists with only one item
        date_begin_forc = date_begin_forc[0]  # replace one-item list by item.
        date_end_forc = date_end_forc[0]
        if os.path.exists('FORCING.nc') or os.path.islink('FORCING.nc'):
            os.remove('FORCING.nc')
        os.symlink(self.options.forcing + '/' + mbdir + '/meteo/FORCING_' + date_begin_forc.strftime('%Y%m%d%H') + '_' + date_end_forc.strftime('%Y%m%d%H') + '.nc',
                   'FORCING.nc')
        # prepare ecoclimap binaries
        if not os.path.exists('ecoclimapI_covers_param.bin'):
            os.symlink(self.exesurfex + '/../MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin', 'ecoclimapI_covers_param.bin')
            os.symlink(self.exesurfex + '/../MY_RUN/ECOCLIMAP/ecoclimapII_eu_covers_param.bin', 'ecoclimapII_eu_covers_param.bin')
            # flanner stuff
            os.symlink(self.exesurfex + '/../MY_RUN//DATA/CROCUS/drdt_bst_fit_60.nc', 'drdt_bst_fit_60.nc')
        if not os.path.exists('offline.exe'):
            os.symlink(self.exesurfex + '/OFFLINE', 'offline.exe')
        if not os.path.exists('PGD.nc'):
            os.symlink(self.xpiddir + 'spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')

        # prepare the namelist with the right escroc options
        self.prepare_namelist_offline(date, idate, mb)

    def prepare_namelist_offline(self, date, idate, mb):
        """
        prepare the namelist (begin/end date, escroc configuration and check DA settings)
        """
        self.prepare_namelist()

        # Prepare the escroc namelist
        # BC copied from vortex.cen.algo.ensemble.py
        physical_options = self.escroc_confs.physical_options[mb - 1]  # to check : I think mb starts at 1
        snow_parameters = self.escroc_confs.snow_parameters[mb - 1]
        # SOOOOO SLOW....
        update_surfex_namelist_file(
            datetime.datetime.strptime(
                self.options.stopdates[idate - 1] if idate - 1 > 0 else self.options.datedeb, '%Y%m%d%H'),
            dateend         = datetime.datetime.strptime(date, '%Y%m%d%H'),
            physicaloptions = physical_options,
            snowparameters  = snow_parameters,
            namelistfile="OPTIONS_base.nam",
            updateloc = False,
            no_caution = True,
        )
        shutil.copyfile("OPTIONS_base.nam", "OPTIONS.nam")
        # check is useless (done uin the mother namelist)
        # check_namelist_soda(self.options)

    def prepare_prep(self, date, dateprev, mb):
        '''
        Prepare the PREP file for OFFLINE from SODA (or OL)
        '''

        if date == self.options.stopdates[0]:
            try:
                os.symlink(self.xpiddir + '/spinup/prep/PREP.nc', 'PREP.nc')
            except Exception:
                os.remove('PREP.nc')
                os.symlink(self.xpiddir + '/spinup/prep/PREP.nc', 'PREP.nc')

        else:
            self.link_build(mb, dateprev)

    def link_build(self, mb, dateprev):
        if self.options.pf and self.options.pf != 'ol':
            try:
                os.symlink(
                    self.xpiddir + '/' + dateprev + '/workSODA/SURFOUT' + str(mb) + '.nc',
                    'PREP.nc')
            except Exception:
                os.remove('PREP.nc')
                os.symlink(
                    self.xpiddir + '/' + dateprev + '/workSODA/SURFOUT' + str(mb) + '.nc',
                    'PREP.nc')
        # ol case
        else:
            try:
                os.symlink(
                    self.xpiddir + '/' + dateprev + '/mb{0:04d}'.format(mb) + '/SURFOUT.nc',
                    'PREP.nc')
            except Exception:
                os.remove('PREP.nc')
                os.symlink(
                    self.xpiddir + '/' + dateprev + '/mb{0:04d}'.format(mb) + '/SURFOUT.nc',
                    'PREP.nc')

    def run(self, date):
        # print('launching escroc until ', date)
        p = multiprocessing.Pool(min(multiprocessing.cpu_count(), self.options.nmembers))
        p.map(self.mb_run, [['/'.join([self.xpiddir, date, mbdir]), mbdir] for mbdir in self.mbdirs])
        p.close()
        p.join()
        # print('escroc step is done.')

    def mb_prepare(self, largs):
        "'' parallelized preparation since it takes sooo much time to edit namelists..."""
        date = largs[0]
        idate = largs[1]
        mbdir = largs[2]
        mb = largs[3]
        self.prepare_offline_env(date, idate, mbdir, mb)
        # prepare the pgd and prep files.
        # in Local parallel sequence, must be done only once the corresponding SODA run has produced the files
        self.prepare_prep(date, self.options.stopdates[idate - 1], mb)

    def mb_run(self, largs):
        path = largs[0]
        mbdir = largs[1]
        os.chdir(path)
        with open('offline.out', 'w') as f:
            subprocess.call('./offline.exe', stdout=f)
