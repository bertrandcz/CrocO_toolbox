# -*- coding: utf-8 -*-
'''
Created on 12 juin 2019

@author: cluzetb
'''

import datetime
import multiprocessing
import os
from snowtools.utils.dates import check_and_convert_date
from snowtools.utils.prosimu import prosimu

from CrocoPf import CrocO
from Ensemble import PrepEnsBg, PrepEnsAn
from Ensemble import PrepEnsOl
from SemiDistributed import FromXp, Real
import numpy as np
import pickle as pickle
from utilcrocO import Pgd, setlistvars_obs, dictvarsPro
from utilcrocO import ftpconnect, area
from utilpp import read_alpha, read_part, read_mask, read_BG, load_pickle2


def loadEnsPrepDate_parallel(largs):
    """
    BC 24/06/20 externalized this funtion for python2 compatibility (if class method, throws this exception):
    https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map/7309686#7309686
    """
    pp = largs[0]
    locEns  = largs[1]
    dd = largs[2]
    kind = largs[3]
    isOl = largs[4]
    fromArch = largs[5]
    if kind == 'bg':
        local = PrepEnsBg(pp.options, dd, fromArch=fromArch)
    elif kind == 'an':
        local = PrepEnsAn(pp.options, dd, fromArch=fromArch)
    else:
        local = PrepEnsOl(pp.options, dd, isOl=isOl, fromArch=fromArch)
    if (kind == 'ol' and isOl is False):
        pathPkl = pp.xpidoldir + '/crocO/' + pp.options.saverep + '/' +\
            kind + '_' + dd + '.pkl'
    else:
        pathPkl = kind + '_' + dd + '.pkl'
    pathpklbeauf = pathPkl + '.foo'
    if not os.path.islink(pathPkl):
        if os.path.exists(pathpklbeauf):
            try:
                os.symlink(pathpklbeauf, pathPkl)
            except FileExistsError:
                pass
    if not os.path.exists(pathPkl):
        print(' unsuccesfull loading', pathPkl, )
        print(('loading ' + kind + ' ens for date : ', dd))
        local.stackit()
        with open(pathPkl, 'wb') as f:
            print(('saaving ' + kind + ' to pickle !'))
            pickle.dump(local.stack, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        local.stack = load_pickle2(pathPkl)
        if not all(l in local.stack.keys() for l in pp.listvar):
            raise Exception('Some of the ppvars you mentioned are not present in ', os.getcwd() + '/' + pathPkl,
                            '\n listvars:', pp.listvar, '\npickle file keys:', local.stack.keys())
        local.isstacked = True
    locEns[dd] = local


class CrocoPp(CrocO):
    '''
    class meant to post-process beaufix cycled assimilation runs (similar to deprecated SodaXP)
    '''

    def __init__(self, options):
        self.options = options
        CrocO.__init__(self, options)
        if self.options.synth:
            self.mbsynth = int(self.options.synth) - 1  # be careful to id offset

        # overwrite with the archive path if it exists
        if hasattr(self.options, 'arch'):
            self.xpiddir = self.options.arch
            self.crocOdir = self.options.arch
        self.mblist = self.options.mblist
        self.initial_context = os.getcwd()

        # setup + loading
        self.setup()
        self.read(readpro = not self.options.notreadpro, readprep=self.options.readprep,
                  readaux = self.options.readaux, readobs = self.options.readobs,
                  readtruth=options.readtruth,
                  readoper = self.options.readoper,)

        if 'notebooks' in self.initial_context:
            os.chdir(self.initial_context)

    def setup(self):
        if not os.path.exists(self.crocOdir + '/' + self.options.saverep):
            os.makedirs(self.crocOdir + '/' + self.options.saverep)
        os.chdir(self.crocOdir + '/' + self.options.saverep)

        # load the PGD (useful for pp)
        if not os.path.exists('PGD.nc') and not os.path.islink('PGD.nc'):
            os.symlink(self.options.crocOpath + '/s2m/' + self.options.vconf + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc', 'PGD.nc')
        self.pgd = Pgd('PGD.nc')

        # set subdirs
        self.subdirs = [self.xpiddir + '/mb{0:04d}'.format(mb) + '/' for mb in self.mblist]
        self.mbdirs = ['mb{0:04d}/'.format(mb) for mb in self.mblist]
        # set dates
        self.setDates()

        # important : check whether the XP is openloop or not.
        if str(self.options.openloop) == 'on':
            if self.options.todo == 'pf':
                self.isOl = False
                self.readOl = False

            else:
                self.isOl = True
                self.readOl = True

            self.xpidoldir = self.xpiddir
        else:
            self.isOl = False
            if hasattr(self.options, 'xpidoldir'):
                self.readOl = True
                self.xpidoldir = self.options.xpidoldir
                if not os.path.exists(self.xpidoldir):
                    raise Exception('the prescribed OL associated experiment ' + self.xpidoldir + 'does not exist')
                else:
                    if not os.path.exists(self.xpidoldir + '/PGD.nc'):
                        os.symlink(self.options.crocOpath + '/s2m/' + self.options.vconf + '/spinup/pgd/PGD_' + area(self.options.vconf) + '.nc',
                                   self.xpidoldir + 'PGD.nc')
            else:
                self.readOl = False

    def setDates(self):
        try:
            self.datedeb = datetime.datetime.strptime(str(self.options.datedeb), "%Y%m%d%H")
            self.datefin = datetime.datetime.strptime(str(self.options.datefin), "%Y%m%d%H")
        except ValueError:
            self.datedeb = datetime.datetime.strptime(str(self.options.datedeb), "%Y-%m-%d %H:%M:%S")
            self.datefin = datetime.datetime.strptime(str(self.options.datefin), "%Y-%m-%d %H:%M:%S")
        self.stopdates = list(map(check_and_convert_date, self.options.stopdates)) if hasattr(self.options, 'stopdates') else list(map(check_and_convert_date, self.options.assimdates + [self.options.datefin.strftime('%Y%m%d%H')]))
        self.assimdates = list(map(check_and_convert_date, self.options.assimdates))

        self.begprodates = [self.datedeb] + self.assimdates
        self.endprodates = self.assimdates + [self.datefin]

    def read(self, readpro = True, readprep = False, readaux = False, readobs = False, readtruth = False, readoper = False, ):

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
        if readtruth is True:
            if not self.isOl and self.readOl:
                try:
                    self.readTruth()
                except AttributeError:
                    raise Exception('please prescribe a synthetic member in either options or conf file.')

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
            if self.options.pf in ['rlocal', 'klocal']:
                self.mask = read_mask(self.options)
                if self.options.pf == 'klocal':
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
        if self.options.todo == 'pfpp' or self.options.todo == 'pf':
            fromArch = False
        else:
            fromArch = True
        # on beaufix, cannot straightforwardly parallelize the reading of preps
        # error: The "default" UnstructuredGeometry object does not exist yet"
        # BC 30/09/20: don get what the hell happens on sxcen.
        if 'beaufix' in os.uname()[1] or 'sxcen' in os.uname()[1] or 'belenos' in os.uname()[1]:
            if kind == 'bg':
                locEns = {dd: PrepEnsBg(self.options, dd, fromArch=fromArch) for dd in self.options.dates}
            elif kind == 'an':
                locEns = {dd: PrepEnsAn(self.options, dd, fromArch=fromArch) for dd in self.options.dates if dd != self.datefin.strftime('%Y%m%d%H')}
            else:
                locEns = {dd: PrepEnsOl(self.options, dd, isOl=isOl, fromArch=fromArch) for dd in self.options.dates}

            for dd in locEns.keys():
                if (kind == 'ol' and isOl is False):
                    pathPkl = self.xpidoldir + '/crocO/' + self.options.saverep + '/' +\
                        kind + '_' + dd + '.pkl'
                else:
                    pathPkl = kind + '_' + dd + '.pkl'
                pathpklbeauf = pathPkl + '.foo'
                if not os.path.islink(pathPkl):
                    if os.path.exists(pathpklbeauf):
                        try:
                            os.symlink(pathpklbeauf, pathPkl)
                        except FileExistsError:
                            pass
                if not os.path.exists(pathPkl):
                    print(' unsuccesfull loading', pathPkl, )
                    print(('loading ' + kind + ' ens for date : ', dd))
                    locEns[dd].stackit()
                    with open(pathPkl, 'wb') as f:
                        print(('saaving ' + kind + ' to pickle !'))
                        pickle.dump(locEns[dd].stack, f, protocol=pickle.HIGHEST_PROTOCOL)
                else:
                    locEns[dd].stack = load_pickle2(pathPkl)
                    if not all(l in locEns[dd].stack.keys() for l in self.listvar):
                        raise Exception('Some of the ppvars you mentioned are not present in ', pathPkl,
                                        '\n listvars:', self.listvar, '\npickle file keys:', locEns[dd].stack.keys())
                    locEns[dd].isstacked = True
            return locEns

        else:
            # BC june 2020: manager is used to share the dict between nodes
            try:
                manager = multiprocessing.Manager()
                locEns = manager.dict()
                p = multiprocessing.Pool(min(max(multiprocessing.cpu_count() - 4, 1), len(self.options.dates)))
                p.map(loadEnsPrepDate_parallel, [[self, locEns, dd, kind, isOl, fromArch] for dd in self.options.dates])
            except Exception as _:
                p.close()
            p.close()
            p.join()
            return dict(locEns)

    def loadEnsPro(self, kind, catPro = False, isOl = False):
        if kind == 'Cl':
            pathPkl = self.xpiddir + '../clim/crocO/clim.pkl'
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

        # convert to pkl or load.
        if not os.path.exists(pathPkl):
            locPro = self.nc_to_pkl(pathPkl, kind, catPro = catPro, isOl = isOl)
        else:
            locPro = load_pickle2(pathPkl)
            if not all(l in locPro.keys() for l in self.listvar):
                print('Some of the ppvars you mentioned are not present in ', pathPkl,
                      '\n listvars:', self.listvar, '\npickle file keys:', [locPro.keys()], '\n reloading')

                # mutualize variables because the pkl will be overwritten.
                lkeys = list(locPro.keys())
                lkeys.remove('time')
                self.listvar = set(lkeys + self.listvar)
                locPro = self.nc_to_pkl(pathPkl, kind, catPro = catPro, isOl = isOl)
        return locPro

    def getProFiles(self):
        """
        Get pro files from hendrix and concatenate it
        """
        print('getting pro files on hendrix')
        gg = os.getcwd()
        os.chdir(self.xpiddir)
        ftpObject = ftpconnect('hendrix')
        for s in self.mbdirs:
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
        # if performing pfpp, it is after a locla run and obs are properly archived.
        if self.options.todo == 'pfpp':
            if self.options.synth is not None:
                self.pathArch = self.xpiddir + '/crocO/ARCH/' + self.options.sensor + '/'
                self.pathSynth = self.xpiddir + '/crocO/SYNTH/' + self.options.sensor + '/'
            else:
                self.pathReal = self.options.crocOpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'

        else:
            if str(self.options.openloop) == 'on':
                if self.options.synth is not None:
                    self.pathArch = self.xpidoldir + '/crocO/ARCH/' + self.options.sensor + '/'
                    self.pathSynth = self.xpidoldir + '/crocO/SYNTH/' + self.options.sensor + '/'
                else:
                    self.pathReal = self.options.crocOpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'
            else:
                if hasattr(self.options, 'xpidoldir'):
                    if self.options.synth is not None:
                        self.pathArch = self.xpidoldir + '/crocO/ARCH/' + self.options.sensor + '/'
                        self.pathSynth = self.xpidoldir + '/crocO/SYNTH/' + self.options.sensor + '/'
                    else:
                        self.pathReal = self.options.crocOpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'
                else:
                    # if no xpidoldir has been prescribed, obs must be real.
                    self.pathReal = self.options.crocOpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor + '/'

        # proper loading (synth or real)
        if self.options.todo == 'pfpp':
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
            if self.options.synth is not None:
                print('sssynth', self.options.synth)
                print('reading synthetic truth obs. in ' + self.pathArch)
                print('reading synthetic assimilated obs. in ' + self.pathSynth)

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

                pathObsTs = self.options.crocOpath + '/s2m/' + self.options.vconf + '/obs/' + self.options.sensor +\
                    '/obs_{0}_{1}_{2}100106_{3}063006.pkl'.format(self.options.sensor,
                                                                  self.options.vconf,
                                                                  self.datedeb.strftime('%Y'),
                                                                  int(self.datedeb.strftime('%Y')) + 1)
                print('reading observation timeseries (pickle from csv file) in ')
                print(pathObsTs)
                if not os.path.exists(pathObsTs):
                    if '_X' in self.options.sensor:
                        sensorBase = self.options.sensor.split('X')[0][0:-1]
                        os.symlink(self.options.crocOpath + '/s2m/' + self.options.vconf +
                                   '/obs/{0}/obs_{0}_{1}_{2}100106_{3}063006.pkl'.format(sensorBase,
                                                                                         self.options.vconf,
                                                                                         self.datedeb.strftime('%Y'),
                                                                                         int(self.datedeb.strftime('%Y')) + 1),
                                   pathObsTs)
                    else:
                        # BC 30/06/20: some day replace with a link to a yearly file/
                        try:
                            os.symlink(self.options.crocOpath + '/s2m/' + self.options.vconf +
                                       '/obs/{0}/obs_all_{0}_2013010106_2018123106.pkl'.format(self.options.vconf,
                                                                                               self.datedeb.strftime('%Y'),
                                                                                               int(self.datedeb.strftime('%Y')) + 1
                                                                                               ),
                                       pathObsTs)
                            self.obsTs = load_pickle2(pathObsTs)
                        except FileExistsError:
                            try:
                                self.obsTs = load_pickle2(pathObsTs)
                            except FileNotFoundError:  # no obsTs file (makes some sense too.)
                                pass
                        except FileNotFoundError:  # no obsTs file (makes some sense too.)
                            pass
                else:
                    self.obsTs = load_pickle2(pathObsTs)

    def readOper(self):
        """
        (postes geom only) read pickle version of the operational PRO files.
        (pickles generated by extract_postes_add_mocage.py)
        """
        pathOper = os.environ['CROCOPATH'] + '/s2m/' + self.options.vconf + '/oper_{0}/crocO/beaufix/oper.pkl'.format(self.begprodates[0].strftime('%Y'))
        self.oper = load_pickle2(pathOper)

    def nc_to_pkl(self, pathPkl, kind, catPro = False, isOl = False):
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
        # @TODO: parallelize the reading of PRO files
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
            datahtn = prosimu(str(f))
            if iS == 0:
                for var in self.listvar:
                    if 'B' not in var:
                        locPro[var] = np.expand_dims(datahtn.read(self.dictvarsPro[var]), axis=2)
                    else:
                        locPro[var] = np.expand_dims(datahtn.read('SPECMOD')[:, int(var[-1]) - 1, :], axis=2)

                locPro['time'] = datahtn.readtime()
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
            print(pathPkl)
            pickle.dump(locPro, f, protocol=pickle.HIGHEST_PROTOCOL)
        return locPro
