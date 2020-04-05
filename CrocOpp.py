# -*- coding: utf-8 -*-
'''
Created on 12 juin 2019

@author: cluzetb
'''

from CrocOrun import CrocO
from Ensemble import PrepEnsBg, PrepEnsAn
from Ensemble import PrepEnsOl
from SemiDistributed import FromXp, Real
from bronx.datagrip.namelist import NamelistParser
import datetime
import os
from utilcrocO import Pgd, setlistvars_obs, dictvarsPro
from utilcrocO import ftpconnect, area
from utilpp import read_alpha, read_part, read_mask, read_BG, load_pickle2
from utils.dates import check_and_convert_date
from utils.prosimu import prosimu

import numpy as np
import pickle as pickle


class CrocOpp(CrocO):
    '''
    class meant to post-process beaufix cycled assimilation runs (similar to deprecated SodaXP)
    '''

    def __init__(self, options):
        self.options = options
        CrocO.__init__(self, options)
        if hasattr(self.options, 'synth'):
            self.mbsynth = int(self.options.synth) - 1  # be careful to id offset

        self.initial_context = os.getcwd()
        # setup + loading
        self.setup()
        self.read(readpro = not self.options.notreadpro, readprep=self.options.readprep, readaux = self.options.readaux, readobs = self.options.readobs,
                  readoper=self.options.readoper)

        if 'notebooks' in self.initial_context:
            os.chdir(self.initial_context)

    def setup(self):
        if not os.path.exists(self.crocodir + '/' + self.options.saverep):
            os.makedirs(self.crocodir + '/' + self.options.saverepself.options.saverep)
        os.chdir(self.crocodir + '/' + self.options.saverep)

        # load the PGD (useful for pp)
        if not os.path.exists('PGD.nc'):
            os.symlink(self.options.vortexpath + '/s2m/' + self.options.vconf + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')
        self.pgd = Pgd('PGD.nc')

        # set subdirs
        self.subdirs = [self.xpiddir + 'mb{0:04d}'.format(mb) + '/' for mb in self.mblist]
        self.mbdirs = ['mb{0:04d}/'.format(mb) for mb in self.mblist]
        # set dates
        self.setDates()

        # important : check whether the XP is openloop or not.
        if str(self.options.openloop) == 'on':
            self.isOl = True
            self.readOl = True
        else:
            self.isOl = False
            self.read_opts_in_namelist()
            if hasattr(self.options, 'xpidoldir'):
                self.readOl = True
                self.xpidoldir = self.options.xpidoldir
                if not os.path.exists(self.xpidoldir):
                    raise Exception('the prescribed OL associated experiment ' + self.xpidoldir + 'does not exist')
                else:
                    if not os.path.exists(self.xpidoldir + '/PGD.nc'):
                        os.symlink(self.options.vortexpath + '/s2m/' + self.options.vconf + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc',
                                   self.xpidoldir + 'PGD.nc')
            else:
                self.readOl = False

    def setDates(self):
        self.datedeb = datetime.datetime.strptime(str(self.options.datedeb), "%Y-%m-%d %H:%M:%S")
        self.datefin = datetime.datetime.strptime(str(self.options.datefin), "%Y-%m-%d %H:%M:%S")
        self.stopdates = list(map(check_and_convert_date, self.options.stopdates)) if hasattr(self.options, 'stopdates') else list(map(check_and_convert_date, self.options.assimdates + [self.options.datefin.strftime('%Y%m%d%H')]))
        self.assimdates = list(map(check_and_convert_date, self.options.assimdates))

        self.begprodates = [self.datedeb] + self.assimdates
        self.endprodates = self.assimdates + [self.datefin]

    def read(self, readpro = True, readprep = False, readaux = False, readobs = False, readoper = False):

        # set the list of vars
        self.listvar = setlistvars_obs(self.options.ppvars)
        self.dictvarsPro = dictvarsPro()
        if readpro is True:
            self.readEnsPro(readOl = self.readOl, readClim = self.options.clim)
        if readprep is True:
            self.readEns()
        # read auxiliary files
        if readaux is True:
            self.readAux()
        if readobs is True:
            self.readObs()
            if hasattr(self.options, 'synth') and not self.isOl and self.readOl:
                self.readTruth()
        if readoper is True:
            self.readOper()

    def readEns(self):
        print('initializing ens')
        if not self.isOl:
            self.ensBg = self.loadEnsPrep('bg')
            self.ensAn = self.loadEnsPrep('an')
        if self.isOl or self.readOl:
            self.ensOl = self.loadEnsPrep('ol', isOl = self.isOl)

    def readAux(self):
        if not self.isOl:
            self.alpha = read_alpha(self.options)
            self.part = read_part(self.options)
            if self.options.pf_crocus == 'klocal':
                self.mask = read_mask(self.options)
                self.BG = read_BG(self.options)

    def readEnsPro(self, catPro = False, readOl = False, readClim = False):
        if not self.isOl:
            self.ensProAn = self.loadEnsPro('An', catPro = catPro, isOl=self.isOl)
        if self.isOl or readOl:
            self.ensProOl = self.loadEnsPro('Ol', catPro = catPro, isOl=self.isOl)
        if readClim:
            self.ensProClim = self.loadEnsPro('Cl', catPro = catPro, isOl=self.isOl)

    def readTruth(self):
        """
        BC 24/10/19 : function to read the truth in the OL from the corresponding run.
        - Sometimes, if the OL is a subset of a bigger OL, and the truth is in this bigger OL, need to read in it.
        k is the variable.
        - find the baseline experiment if necessary.
        """

        def load_from_dict(ddict, mbsynth):
            out = dict()
            for var in ddict.keys():
                if 'time' not in var:
                    out[var] = ddict[var][:, :, mbsynth]
                else:
                    out['time'] = ddict['time']
            return out
        try:
            self.truth = load_from_dict(self.ensProOl, self.mbsynth)
        except IndexError:
            print('\n\nWARNING : there is no corresponding mbsynth in this openloop not enough members\n looking in the bigger OL xp.\n\n')
            print('loading ' + self.xpidoldir[0:-4] + '/crocO/' + self.options.saverep + '/ensProOl.pkl')
            pathPklBigOl = self.xpidoldir[0:-4] + '/crocO/' + self.options.saverep + '/ensProOl.pkl'
            gg = load_pickle2(pathPklBigOl)
            self.truth = load_from_dict(gg, self.mbsynth)

    def loadEnsPrep(self, kind, isOl = False):

        # if local pp (e.g. pp of a run on local machine)
        if self.options.kind == 'localpp':
            directFromXp = False
        else:
            directFromXp = True

        if kind == 'bg':
            locEns = {dd: PrepEnsBg(self.options, dd, directFromXp=directFromXp) for dd in self.options.dates}
        elif kind == 'an':
            locEns = {dd: PrepEnsAn(self.options, dd, directFromXp=directFromXp) for dd in self.options.dates}
        else:
            locEns = {dd: PrepEnsOl(self.options, dd, isOl=isOl, directFromXp=directFromXp) for dd in self.options.dates}

        for dd in self.options.dates:
            if (kind == 'ol' and isOl is False):
                pathPkl = self.xpidoldir + '/crocO/' + self.options.saverep + '/' +\
                    kind + '_' + dd + '.pkl'
            else:
                pathPkl = kind + '_' + dd + '.pkl'
            pathpklbeauf = pathPkl + '.foo'
            if not os.path.islink(pathPkl) or not os.path.exists(pathPkl):
                if os.path.exists(pathpklbeauf):
                    try:
                        os.symlink(pathpklbeauf, pathPkl)
                    except FileExistsError:
                        pass
                else:
                    print(('loading ' + kind + ' ens for date : ', dd))
                    locEns[dd].stackit()
                    with open(pathPkl, 'wb') as f:
                        print(('saaving ' + kind + ' to pickle !'))
                        pickle.dump(locEns[dd].stack, f, protocol=pickle.HIGHEST_PROTOCOL)
            else:
                locEns[dd].stack = load_pickle2(pathPkl)
                locEns[dd].isstacked = True
        return locEns

    def loadEnsPro(self, kind, catPro = False, isOl = False):
        if kind == 'Cl':
            pathPkl = self.xpiddir + '../clim/crocO/clim.pkl'
            print(pathPkl)
        elif (kind == 'Ol' and isOl is False):
            if not os.path.exists(self.xpidoldir + '/crocO/' + self.options.saverep + '/'):
                os.makedirs(self.xpidoldir + '/crocO/' + self.options.saverep + '/')
            pathPkl = self.xpidoldir + '/crocO/' + self.options.saverep + '/' + 'ensPro' + kind + '.pkl'
        else:
            pathPkl = 'ensPro' + kind + '.pkl'
        # with pickling on beaufix (february 2020 on), pickle files are in xpdir/crocO and with .foo extension added.
        pathpklbeauf = pathPkl + '.foo'
        if not os.path.islink(pathPkl):
            if os.path.exists(pathpklbeauf):
                try:
                    os.symlink(pathpklbeauf, pathPkl)
                except FileExistsError:
                    pass
        if not os.path.exists(pathPkl):
            print(('preparing the PRO ' + kind + ' pickle file'))
            locPro = dict()
            # DOWNLOADING, concatenating, reading.
            if kind == 'An' or isOl:
                memberdirs = self.subdirs
            else:
                memberdirs = [self.xpidoldir + '/mb{0:04d}'.format(mb) + '/' for mb in self.mblist]

            if not memberdirs[0] + 'pro/':
                self.getProFiles()

            # then check for concatenated file, if not, create it and delete the sequence of files.
            for iS, s in enumerate(memberdirs):
                print(('reading member ' + str(iS + 1)))

                procat = s + 'pro/' + 'PRO_' + self.datedeb.strftime("%Y%m%d%H") + '_' + self.datefin.strftime("%Y%m%d%H") + '.nc'
                if catPro is True:

                    if not os.path.exists(procat):

                        cmd = ' '.join([s + 'pro/' + 'PRO_' + dd.strftime("%Y%m%d%H") + '_' + df.strftime("%Y%m%d%H") + '.nc'
                                        for (dd, df) in zip(self.begprodates, self.endprodates)])
                        print(('concatenating PRO in member ' + str(iS + 1)))
                        os.system('ncrcat ' + cmd + ' ' + procat)
                        print(('removing PRO sequence in member ' + str(iS + 1)))
                        for (dd, df) in zip(self.begprodates, self.endprodates):
                            gg = s + 'pro/' + 'PRO_' + dd.strftime("%Y%m%d%H") + '_' + df.strftime("%Y%m%d%H") + '.nc'
                            if os.path.exists(gg):
                                os.remove(gg)
                    print(('stackingHTN ', str(iS + 1)))
                    f = s + '/pro/' + 'PRO_' + self.datedeb.strftime("%Y%m%d%H") + '_' + self.datefin.strftime("%Y%m%d%H") + '.nc'

                else:

                    # if provided a list of files, prosimu automatically aggregates it along time !!
                    if not os.path.exists(procat):
                        f = s + 'pro/' + 'PRO*.nc'
                    else:
                        f = s + 'pro/' + 'PRO_' + self.datedeb.strftime("%Y%m%d%H") + '_' + self.datefin.strftime("%Y%m%d%H") + '.nc'

                # read
                print('fffileeaf', f)
                datahtn = prosimu(str(f))
                if iS == 0:
                    print('llistvars', self.listvar)
                    for var in self.listvar:
                        if 'B' not in var:
                            locPro[var] = np.expand_dims(datahtn.read(self.dictvarsPro[var]), axis=2)
                        else:
                            locPro[var] = np.expand_dims(datahtn.read('SPECMOD')[:, int(var[-1]) - 1, :], axis=2)

                    locPro['time'] = datahtn.readtime()
                    print(('lol', list(locPro.keys())))
                else:
                    for var in self.listvar:
                        if 'B' not in var:
                            locPro[var] = np.concatenate( (locPro[var], np.expand_dims(datahtn.read(self.dictvarsPro[var]), axis=2)), axis=2)
                        else:
                            locPro[var] = np.concatenate((locPro[var], np.expand_dims(datahtn.read('SPECMOD')[:, int(var[-1]) - 1, :], axis=2)), axis = 2)

                datahtn.close()
            # SAVING
            with open(pathPkl, 'wb') as f:
                print(('saaving ensPro' + kind + ' to pickle !'))
                pickle.dump(locPro, f, protocol=pickle.HIGHEST_PROTOCOL)

        else:
            locPro = load_pickle2(pathPkl)

        return locPro

    def getProFiles(self):
        """
        Get pro files from hendrix and concatenate it
        """
        print('getting pro files on hendrix')
        print((os.getcwd()))
        gg = os.getcwd()
        os.chdir(self.xpiddir)
        ftpObject = ftpconnect('hendrix')
        for s in self.mbdirs:
            print((s, 'on hendrix'))
            print(('current', os.getcwd()))
            os.chdir(s)
            if not os.path.exists('pro'):
                os.mkdir('pro')
            os.chdir('pro')

            for (dd, df) in zip(self.begprodates, self.endprodates):
                ftpObject.retrbinary('RETR /home/cluzetb/vortex/' + self.options.vapp + '/' + self.options.vconf +
                                     '/' + self.options.xpid[0:self.xpid.find('@')] + '/' + s + 'pro/' + 'PRO_' + dd.strftime("%Y%m%d%H") +
                                     '_' + df.strftime("%Y%m%d%H") + '.nc', open('PRO_' + dd.strftime("%Y%m%d%H") +
                                                                                 '_' + df.strftime("%Y%m%d%H") + '.nc', 'wb').write)
            os.chdir('../..')
        os.chdir(gg)

    def readObs(self):
        # paths settings (ol or not, real obs or not).
        # if performing localpp, it is after a locla run and obs are properly archived.
        if self.options.kind == 'localpp':
            if self.options.synth is not None:
                self.pathArch = self.xpiddir + '/crocO/ARCH/' + self.options.sensor + '/'
                self.pathSynth = self.xpiddir + '/crocO/SYNTH/' + self.options.sensor + '/'
            else:
                self.pathReal = self.options.vortexpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'

        else:
            if str(self.options.openloop) == 'on':
                if hasattr(self.options, 'synth'):
                    self.pathArch = self.xpidoldir + '/crocO/ARCH/' + self.options.sensor + '/'
                    self.pathSynth = self.xpidoldir + '/crocO/SYNTH/' + self.options.sensor + '/'
                else:
                    self.pathReal = self.options.vortexpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'
            else:
                if hasattr(self.options, 'xpidoldir'):
                    if hasattr(self.options, 'synth'):
                        self.pathArch = self.xpidoldir + '/crocO/ARCH/' + self.options.sensor + '/'
                        self.pathSynth = self.xpidoldir + '/crocO/SYNTH/' + self.options.sensor + '/'
                    else:
                        self.pathReal = self.options.vortexpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'
                else:
                    # if no xpidoldir has been prescribed, obs must be real.
                    self.pathReal = self.options.vortexpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'

        # proper loading (synth or real)
        if self.options.kind == 'localpp':
            if self.options.synth is not None:
                print('reading synthetic assimilated obs. in ' + self.pathArch)
                self.obsArch = {dd: FromXp(self.pathArch, dd, self.options) for dd in self.options.dates}
                self.obsSynth = {dd: FromXp(self.pathSynth, dd, self.options) for dd in self.options.dates}
                for dd in self.options.dates:
                    self.obsArch[dd].load()
                    self.obsSynth[dd].load()
            else:
                print('reading real assimilated obs. in ' + self.pathReal)
                self.obsReal = {dd: Real(self.pathReal, dd, self.options) for dd in self.options.dates}
                for dd in self.options.dates:
                    self.obsReal[dd].load()

        else:
            if hasattr(self.options, 'synth'):
                print('reading synthetic assimilated obs. in ' + self.pathArch)
                self.obsArch = {dd: FromXp(self.pathArch, dd, self.options) for dd in self.options.dates}
                self.obsSynth = {dd: FromXp(self.pathSynth, dd, self.options) for dd in self.options.dates}
                for dd in self.options.dates:
                    self.obsArch[dd].load()
                    self.obsSynth[dd].load()
            else:
                print('reading real assimilated obs. in ' + self.pathReal)
                self.obsReal = {dd: Real(self.pathReal, dd, self.options) for dd in self.options.dates}
                for dd in self.options.dates:
                    self.obsReal[dd].load()
                print('reading observation timeseries (pickle from csv file)')
                pathObsTs = self.options.vortexpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor +\
                    '/obs_{0}_{1}_2013010106_2018123106.pkl'.format(self.options.sensor, self.options.vconf)
                if not os.path.exists(pathObsTs):
                    os.symlink(self.options.vortexpath + '/s2m/' + self.options.vconf +
                               '/obs/all/obs_all_{0}_2013010106_2018123106.pkl'.format(self.options.vconf), pathObsTs)
                self.obsTs = load_pickle2(pathObsTs)

    def readOper(self):
        # add an optio one day...
        if self.options.kind != 'localpp':
            pathOper = self.options.vortexpath + '/s2m/' + self.options.vconf + '/oper_{0}/crocO/beaufix/oper.pkl'.format(self.begprodates[0].strftime('%Y'))
            self.oper = load_pickle2(pathOper)

    def read_opts_in_namelist(self):
        n = NamelistParser()
        if self.options.kind != 'localpp':
            if os.path.exists(self.options.xpiddir + 'conf/namelist.surfex.foo'):
                N = n.parse(self.options.xpiddir + 'conf/namelist.surfex.foo')
            else:
                N = n.parse(self.options.xpiddir + '/conf/OPTIONS.nam')
        else:
            print('cwdededed', os.getcwd())
            N = n.parse(self.options.dates[0] + '/OPTIONS.nam')
        if hasattr(N['NAM_ASSIM'], 'NLOC_PF'):
            self.options.nloc_pf = N['NAM_ASSIM'].NLOC_PF
            self.options.neff_pf = N['NAM_ASSIM'].NEFF_PF
            self.options.pf_crocus = N['NAM_ASSIM'].CPF_CROCUS
