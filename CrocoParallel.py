# -*- coding: utf-8 -*-
'''
Created on 26 mars 2020

@author: cluzetb

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
from snowtools.tools.change_prep import prep_tomodify
from snowtools.tools.update_namelist import update_surfex_namelist_file
from snowtools.utils.ESCROCsubensembles import ESCROC_subensembles
from snowtools.utils.dates import get_list_dates_files

from tqdm import tqdm

from CrocoPf import CrocO, CrocoPf
from utilcrocO import area, check_namelist_soda, dump_conf, safe_create_link


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
        if self.options.pf is None:
            raise Exception('please specify a pf variant, or openloop')
        self.dump_options()
        # duplicate with offline pool... but cannot move to CrocO class yet (due tu the synth case removal)
        self.mblist = list(range(1, self.options.nmembers + 1))
        self.mbdirs = ['mb{0:04d}'.format(mb) + '/' for mb in self.mblist]

        self.get_spinup()

        self.setup()

    def dump_options(self):
        '''
        rough dump of :
        - options into conf file in conf dir
        - namelist copy into conf dir
        '''
        _ = dump_conf(self.xpiddir + '/conf/s2m_' + self.options.vconf + '.ini', self.options)

        # the namelist is roughly copied and receives
        # common necessary modifications (updateloc ?)
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
        # Even though this namelist is not used by soda, it must be compatible with the subsequent SODA run.
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
        '''
        setup all the simulation dirs (all dates at once !)
        (same simulation architecture as on beaufix (a little bit simpler maybe)
        feed them with the constants
        prepare the observations and put it into the dirs
        => CrocoPf is used to prepare all that (it can be fed with a list of dates !!)
        however, the links to prep must be done just before the soda run itself
        hence, they are encapsulated inside run_parallel() class method
        '''
        if self.options.pf != 'ol':
            self.sodas = CrocoPf(self.options)  # first because prepare the directories.
        self.soda_time = time.time() - self.start_time
        self.escrocs = OfflinePools(self.options)
        self.escroc_time = time.time() - self.start_time - self.soda_time
        self.setup_time = time.time() - self.start_time
        print('setup duration:', self.setup_time)
        if self.options.pf != 'ol':
            print('|         soda:', self.soda_time)
        print('|       escroc:', self.escroc_time)

    def run(self, cleanup=False):
        '''
        run, archive, and (optionally) cleanup.
        '''
        # progress_bar
        pbar = tqdm(self.options.stopdates)
        for dd in pbar:
            pbar.set_description('ESCROC simulation until ' + dd)
            #  - spawn offline
            self.escrocs.run(dd)

            # - spawn soda (if not openloop and not on the last timestep.
            if self.options.pf != 'ol' and dd != self.options.datefin:
                self.sodas.run_parallel(dd)
        # report on the time spent
        elapsed_time = time.time() - self.start_time
        print('elapsed time(setup and simu) :', elapsed_time)

        # archive
        self.archive()
        # post-process the outputs.
        # cleanup
        if cleanup is True:
            self.cleanup()

    def archive(self):
        '''
        archive simulation outputs. could be parallelized.
        /!\Erases previous archives on the same path
        '''
        start_time = time.time()
        if self.options.arch is None:
            print("putting the archive in xpid. Not recommended")
            self.options.arch = self.xpid

        print('archiving the outputs to ', self.options.arch)
        if os.path.exists(self.options.arch):
            shutil.rmtree(self.options.arch)
        os.makedirs(self.options.arch)
        os.chdir(self.options.arch)

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
            # trap here !! idate starts at 0 instead of 1.
            # @TODO :integration test should check the PREP dates correspond to their names.
            idate = idd + 1
            [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/SURFOUT.nc',
                             mbdir + '/bg/PREP_' + date + '.nc') for mbdir in self.mbdirs]
            if self.options.pf != 'ol':
                [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/PREP.nc',
                                 mbdir + '/an/PREP_' + self.options.stopdates[idate - 1] + '.nc') for mbdir in self.mbdirs]
            [shutil.copyfile('/'.join([self.xpiddir, date, mbdir]) + '/ISBA_PROGNOSTIC.OUT.nc',
                             mbdir + '/pro/PRO_' + self.options.stopdates[idate - 1] + '_' + date + '.nc') for mbdir in self.mbdirs]
        if self.options.pf != 'ol':
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

        self.mblist = list(range(1, self.options.nmembers + 1))
        self.mbdirs = ['mb{0:04d}'.format(mb) + '/' for mb in self.mblist]
        # some times we may want to use less forcings than members
        # (e.g.: 1 2 3 1 2 3 1 2 3, nforcing=3, nmembers = 9).
        self.mbdirs_forc = ['mb{0:04d}'.format(((mb - 1) % self.options.nforcing) + 1) + '/' for mb in self.mblist]

        # setup
        self.setup()

    def setup(self):
        '''
        Setup the execution directory for every member, every date.
        Implies time-consuming editing of the namelists for OFFLINE.
        For that reason, this process is parallelised.
        '''
        # prepare escroc configurations
        self.escroc_confs = ESCROC_subensembles(self.options.escroc, self.options.members_id)
        size_setpool = min(multiprocessing.cpu_count(), self.options.nmembers * len(self.options.stopdates))
        p = multiprocessing.Pool(size_setpool)
        p.map(self.mb_prepare, [[date, idate, mbdir, mb, mbdirs_forc] for idate, date in enumerate(self.options.stopdates)
                                for (mb, mbdir, mbdirs_forc) in zip(self.mblist, self.mbdirs, self.mbdirs_forc)])
        p.close()
        p.join()
        os.chdir(self.xpiddir)

    def prepare_offline_env(self, date, idate, mbdir, mb, mbdirs_forc):
        '''
        set offline environment for each date (=path): -PGD, links to preps, namelist, ecoclimap etc.
        '''
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
        try:
            safe_create_link(self.options.forcing + '/' + mbdirs_forc + '/meteo/FORCING_' + date_begin_forc.strftime('%Y%m%d%H') + '_' + date_end_forc.strftime('%Y%m%d%H') + '.nc',
                             'FORCING.nc')
        except IOError:
            # BC dirty fix to load forcings with custom begin/end dates.
            safe_create_link(glob.glob(self.options.forcing + '/' + mbdirs_forc + '/meteo/*')[0],
                             'FORCING.nc')

        # prepare ecoclimap binaries

        safe_create_link(self.exesurfex + '/../MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin', 'ecoclimapI_covers_param.bin')
        safe_create_link(self.exesurfex + '/../MY_RUN/ECOCLIMAP/ecoclimapII_eu_covers_param.bin', 'ecoclimapII_eu_covers_param.bin')
        # flanner stuff
        safe_create_link(self.exesurfex + '/../MY_RUN//DATA/CROCUS/drdt_bst_fit_60.nc', 'drdt_bst_fit_60.nc')
        safe_create_link(self.exesurfex + '/OFFLINE', 'offline.exe')
        safe_create_link(self.xpiddir + 'spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')

        # prepare the namelist with the right escroc options
        self.prepare_namelist_offline(date, idate, mb)

    def prepare_namelist_offline(self, date, idate, mb):
        '''
        prepare the namelist (begin/end date, escroc configuration and check DA settings)
        '''
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

    def prepare_prep(self, date, dateprev, mb):
        '''
        Prepare the PREP file for OFFLINE from SODA (or OL)
        '''

        if date == self.options.stopdates[0]:
            safe_create_link(self.xpiddir + '/spinup/prep/PREP.nc', 'PREP.nc')

        else:
            self.link_build(mb, dateprev)

    def link_build(self, mb, dateprev):
        '''
        SODA: The SURFOUT{1..NENS}.nc from the previous run are used as init for the next timestep
        openloop case: take initial states from the corresponding previous mb run.
        '''
        # this links are broken on creation, but exist once SURFOUT have been created.
        if self.options.pf != 'ol':
            safe_create_link(
                self.xpiddir + '/' + dateprev + '/workSODA/SURFOUT' + str(mb) + '.nc',
                'PREP.nc', exc_broken = False)
        # ol case
        else:
            safe_create_link(
                self.xpiddir + '/' + dateprev + '/mb{0:04d}'.format(mb) + '/SURFOUT.nc',
                'PREP.nc', exc_broken = False)

    def run(self, date):
        '''
        Proper parallelized OFFLINE run for a given date.
        '''
        size_process_pool = min(multiprocessing.cpu_count(), self.options.nmembers)
        print('number of cores used:', size_process_pool)
        p = multiprocessing.Pool(size_process_pool)
        p.map(self.mb_run, [['/'.join([self.xpiddir, date, mbdir])] for mbdir in self.mbdirs])
        p.close()
        p.join()

    def mb_prepare(self, largs):
        '''
        parallelized preparation since it takes sooo much time to edit namelists...
        '''
        date = largs[0]
        idate = largs[1]
        mbdir = largs[2]
        mb = largs[3]
        mbdirs_forc = largs[4]
        self.prepare_offline_env(date, idate, mbdir, mb, mbdirs_forc)
        # prepare the pgd and prep files.
        # in Local parallel sequence, must be done only once the corresponding SODA run has produced the files
        self.prepare_prep(date, self.options.stopdates[idate - 1], mb)

    def mb_run(self, largs):
        path = largs[0]
        os.chdir(path)
        with open('offline.out', 'w') as f:
            gg = subprocess.call('./offline.exe', stdout=f, stderr=f)
        if gg != 0:
            print('an ESCROC member crashed, check:', largs[0], 'offline.out')
