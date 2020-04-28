# -*- coding: utf-8 -*-
'''
Created on 5 févr. 2019

@author: cluzetb
Classes manipulation Semi-distributed objects :
- preparing/faking/ synthetic observations observations within crocO framework
- loading prep files/obervations
- ...

'''
import os
from utilcrocO import Pgd, convertdate
from utilcrocO import area
from utilcrocO import setlistvars_obs, setlistvars_var, setSubsetclasses,\
    dictvarsPrep
from utils.prosimu import prosimu

import netCDF4

import numpy as np


class SemiDistributed(object):
    '''
    class for semi-distributed files (obs or PREP) bound to a geometry described by a pgd in current crocOdir
    '''
    _abstract = True

    def __init__(self, pgdPath = 'PGD.nc'):
        self.pgd = Pgd(pgdPath)
        self.isloaded = False


class Obs(SemiDistributed):
    '''
    class attached to an observation file, either synthetical or real
    '''
    _abstract = True

    def __init__(self, date, options):
        self.pgd = options.pgd
        self.isloaded = False
        self.options = options
        self.sodaName = 'OBSERVATIONS_' + convertdate(date).strftime('%y%m%dH%H') + '.nc'
        self.vortexname = 'obs_' + options.sensor + '_' + area(options.vconf) + '_' + date + '.nc'
        # set list of vars
        self.listvar = setlistvars_obs(set(options.ppvars + options.vars))
#         if type(options.vars) is not list:
#             self.listvar = [options.vars]
#         else:
#             self.listvar = options.vars
        # self.listvars_write = setlistvars_obs(options.vars)

    def load(self):
        if not self.isloaded:
            dd = prosimu(self.sodaName)
            self.data = dict()
            for var in self.listvar:
                if 'R' not in var:
                    self.data[var] = dd.read(self.loadDict[var])
                else:
                    if var in list(dd.variables.keys()):
                        self.data[var] = dd.read(self.loadDict[var])
                    else:
                        self.data[var] = dd.read(self.loadDict['B' + var[1]]) / dd.read(self.loadDict['B' + var[2]])
            dd.close()
            self.isloaded = True

    def prepare(self, archive_synth = False, no_need_masking = False):
        if not no_need_masking:

            self.load_raw()
            self.check_format()
            # first,  if it doesn't already exist, create an archive version (/!\ without classesMask).
            # BC 24/02 archive versio is not necessary for real obs.
            if isinstance(self, Real):
                self.prepare_realmask()
                self.create_new(self.Raw, self.maskpath, self.options, fromArch=True, maskit = True)
            else:
                archOk = self.prepare_archive(archive_synth=archive_synth)
                if not archOk:
                    Arch = self.create_new(self.Raw, self.archpath, self.options, maskit = False)
                    # add noise to it and close.
                    if self.options.noise is not None:
                        self.add_noise(Arch, self.options.noise)
                    else:
                        raise Exception('please use --noise option to specify a gaussian noise to add to synthetical obs (eventually 0).')
                    Arch.close()
                # second create a masked copy of the archive into the sodadir :
                # loading archive and masking:
                self.load_arch()
                print('masking of archive from', self.archpath)
                maskArch = self.create_new(self.Arch, self.sodaName, self.options, fromArch =True, maskit=True)
                maskArch.close()

                # third, create a synth archive if necessary:
                if archive_synth:
                    synth = self.create_new(self.Arch, self.synthpath, self.options, fromArch =True, maskit=True)
                    synth.close()
                self.close()
        else:
            # if no masking/modif is needed, just link the file:
            if not os.path.exists or not os.path.islink(self.sodaName):
                os.symlink(self.path, self.sodaName)

    def prepare_realmask(self):
        self.maskpath = self.options.xpidobsdir + '/../' + self.options.sensor + '/' + self.vortexname
        if not os.path.exists(self.options.xpidobsdir + '/../' + self.options.sensor + '/'):
            os.mkdir(self.options.xpidobsdir + '/../' + self.options.sensor + '/')

    def prepare_archive(self, archive_synth=False):

        if not os.path.exists('../../ARCH'):
            os.mkdir('../../ARCH')
        newRep = '../../ARCH/' + self.options.sensor + '/'

        if not os.path.exists(newRep):
            os.mkdir(newRep)
        self.archpath = newRep + self.vortexname

        # if archive_synth is True, archive also the masked version (for synthetic experiments)
        if archive_synth:
            if not os.path.exists('../../SYNTH'):
                os.mkdir('../../SYNTH')
            newRepSynth = '../../SYNTH/' + self.options.sensor + '/'
            if not os.path.exists(newRepSynth):
                os.mkdir(newRepSynth)
            self.synthpath = newRepSynth + self.vortexname
        else:
            self.archpath = self.path
        gg = os.path.exists(self.archpath)
        if gg:
            print('loading archive from', self.archpath)
        else:
            print('archiving to', self.archpath)
        return gg

    def load_raw(self):
        self.Raw = netCDF4.Dataset(self.path, 'r')

    def load_arch(self):
        self.Arch = netCDF4.Dataset(self.archpath, 'r')

    def check_format(self):
        pass

    def create_new(self, Raw, newFile, options, fromArch = False, maskit = False):
        New = netCDF4.Dataset(newFile, 'w')
        if maskit:
            print('masking into', newFile)
            _, mask = self.subset_classes(self.pgd, options)
            self.copydimsvars(Raw, New, self.listvar, fromArch = fromArch, mask=mask)
            self.computeratio(Raw, New, self.listvar, mask=mask)
        else:
            self.copydimsvars(Raw, New, self.listvar,)
        return New

    def subset_classes(self, pgd, options):
        """
        in the SODA copy of the observation file, remove classes where the assim shouldn't be performed
        """

        if options.classes_id is None:  # user can specify a list of classes instead of a selection by Elev,A,S
            subset, mask = setSubsetclasses(pgd, options.classes_e, options.classes_a, options.classes_s)
        else:
            subset = options.classes_id
            mask = np.array([True if i in list(map(int, options.classes_id)) else False for i in range(self.pgd.npts)])
        return subset, mask

    def copydimsvars(self, Raw, New, nameVars, fromArch = False, mask = None):
        '''
        BC 5/02/19 from maskSentinel2.py
        copy nameVars from Raw to New netCDF Datasets (and copy dims before)
        /!\ only for vars that already exist in Raw
        also :
        - apply a class mask if necessary.
        - apply a ceil on SWE values (20mm) if necessary.
        '''
        for dimName, dim in list(Raw.dimensions.items()):
            New.createDimension(dimName, len(dim) if not dim.isunlimited() else None)
        for name in nameVars:  # name in arg format (b*, r**...)
            if not fromArch:
                readName = self.dictVarsRead[name]
            else:
                readName = self.dictVarsWrite[name]
            writeName = self.dictVarsWrite[name]
            if readName in list(Raw.variables.keys()):  # if var exists in Raw (and if it is not a ratio)
                var = Raw.variables[readName]
                tmpvar = var[:]
            elif readName == 'DEP':
                var = Raw.variables['WSN_VEG1']
                tmpvar = np.nansum([(Raw.variables['WSN_VEG' + str(i + 1)][:] / Raw.variables['RSN_VEG' + str(i + 1)][:]) /
                                    np.cos(np.arctan(self.pgd.slope[:]))
                                    for i in range(0, 50)], axis = 0)
            elif readName == 'SWE':
                var = Raw.variables['WSN_VEG1']  # for the reading of attributes
                tmpvar = np.nansum([Raw.variables['WSN_VEG' + str(i + 1)] /
                                    np.cos(np.arctan(self.pgd.slope[:]))
                                    for i in range(0, 50)], axis = 0)
            if readName in list(Raw.variables.keys()) or readName in ['DEP', 'SWE']:
                tmp = New.createVariable(writeName, 'float', var.dimensions)
                # copy attributes
                tmp.setncatts({k: var.getncattr(k) for k in var.ncattrs()})

                if mask is None:
                    tmp[:] = tmpvar[:]
                else:
                    if len(var.shape) == 2.:  # synthetic synthetic (patch dimension)
                        tmp[0, mask] = tmpvar[0, mask]
                        tmp[0, np.invert(np.squeeze(mask))] = var.getncattr('_FillValue')  # _FillValue is 1e+20 in the PREP files, MUST be equal XUNDEF.
                    else:
                        tmp[mask] = tmpvar[mask]
                        tmp[np.invert(mask)] = var.getncattr('_FillValue')
            else:
                raise Exception('the var ' + readName + ' does not exist in the ref file')

    def computeratio(self, Raw, New, nameVars, mask = None):
        '''
        @TODO : potentially deprecated (27/05/19)
        '''

        listrat = [r for r in nameVars if 'R' in r]
        for r in listrat:
            varnum = Raw.variables['B' + r[1]]
            varden = Raw.variables['B' + r[2]]
            tmp = New.createVariable(r, 'float', varnum.dimensions)
            if mask is None:
                tmp[:] = varnum[:] / varden[:]
            else:
                tmp[mask] = varnum[mask] / varden[mask]
                tmp[np.invert(mask)] = np.nan

    def close(self):
        self.Raw.close()
        self.Arch.close()


class Synthetic(Obs):
    """
    Class describing synthetic obs
    """

    def __init__(self, xpiddir, date, options, nmembers = 35, ):
        '''
        Constructor
        '''
        Obs.__init__(self, date, options)
        self.date = date
        if options.synth > 0:
            self.path = xpiddir + 'mb{0:04d}'.format(options.synth) + '/bg/PREP_' + date + '.nc'
            self.ptinom = 'synth' + 'mb{0:04d}'.format(options.synth)
            self.member = options.synth
        self.dictVarsRead = dictvarsPrep()
        self.dictVarsWrite = {name: name for name in self.listvar}
        self.loadDict = self.dictVarsWrite
        self.isloaded = False

    def add_noise(self, New, noise):
        if noise > 0:
            for var in self.options.vars:  # do NOT iterate on options.ppvars
                New.variables[var][:] = New.variables[var][:] * (1 + np.random.normal(0, noise, size=np.shape(New.variables[var][:])))


class Real(Obs):
    """
    Class describing real obs
    """

    def __init__(self, xpidobsdir, date, options):
        '''
        Constructor
        '''
        Obs.__init__(self, date, options)
        if self.options.todo == 'generobs':
            self.sodaName = xpidobsdir + 'obs_' + options.xpidobs + '_' + area(options.vconf) + '_' + date + '.nc'
        else:
            # BC 30/03/20 change that could cause bugs
            # self.sodaName = xpidobsdir + self.vortexname
            self.sodaName = options.xpiddir + '/' + date + '/workSODA/' + self.sodaName
            self.vortexname = xpidobsdir + self.vortexname
        # BC30/03/20 idem
        # self.path = self.sodaName
        self.path = self.vortexname
        self.dictVarsRead = {name: name for name in options.vars}
        self.dictVarsWrite = self.dictVarsRead
        self.loadDict = self.dictVarsRead
        self.isloaded = False
        self.ptinom = xpidobsdir


class Archived(Obs):

    def __init__(self, path, date, options, ptinom = 'archive'):
        Obs.__init__(self, date, options)
        self.path = path
        self.sodaName = self.path + 'crocO/ARCH/' + options.sensor + '/' + self.vortexname
        self.dictVarsRead = {name: name for name in self.listvar}
        self.dictVarsWrite = self.dictVarsRead
        self.loadDict = self.dictVarsWrite
        self.isloaded = False
        self.ptinom = ptinom


class FromXp(Obs):

    def __init__(self, path, date, options, ptinom = 'fromXp'):
        Obs.__init__(self, date, options)
        self.path = path
        self.sodaName = self.path + self.vortexname
        self.dictVarsRead = {name: name for name in self.listvar}
        self.dictVarsWrite = self.dictVarsRead
        self.loadDict = self.dictVarsWrite
        self.isloaded = False
        self.ptinom = ptinom


class Prep(SemiDistributed):
    """
    class describing Prep (background) and/or analysis
    """
    _abstract = True

    def __init__(self, options):
        self.pgd = options.pgd
        self.isloaded = False
        self.options = options
        self.listvar = set(options.vars + options.ppvars)
        # if 'DEP' not in self.listvar:
        #    self.listvar.append('DEP')
        self.listvar_soda = setlistvars_var(options.vars)
        self.dictVarsRead = dictvarsPrep()
        self.dictVarsWrite = {name: name for name in options.vars}
        self.loadDict = dictvarsPrep()

    def load(self):
        if not self.isloaded:
            dd = prosimu(self.sodaName)
            self.data = dict()
            for var in self.listvar:
                if var == 'DEP':
                    self.data[var] = np.nansum([dd.read('WSN_VEG' + str(i + 1)) /
                                                dd.read('RSN_VEG' + str(i + 1)) /
                                                np.cos(np.arctan(self.pgd.slope))
                                                for i in range(0, 50)], axis = 0)
                    # print(self.data[var])
                elif var == 'SWE':
                    self.data[var] = np.nansum([dd.read('WSN_VEG' + str(i + 1)) /
                                                np.cos(np.arctan(self.pgd.slope))
                                                for i in range(0, 50)], axis = 0)
                    # print(self.data[var])
                elif 'R' in var:
                    self.data[var] = dd.read(self.loadDict['B' + var[1]]) / dd.read(self.loadDict['B' + var[2]])
                else:
                    self.data[var] = dd.read(self.loadDict[var])
            dd.close()
            self.isloaded = True


class PrepBg(Prep):
    def __init__(self, date, mbid, options, directFromXp = True):
        Prep.__init__(self, options)
        self.date = date
        self.ptinom = 'bg' + str(mbid)
        if not directFromXp:
            self.sodaName = date + '/PREP_' + convertdate(date).strftime('%y%m%dH%H') + '_PF_ENS' + str(mbid) + '.nc'
        else:
            # self.sodaName = '../../../test_160_OL@cluzetb/' + 'mb{0:04d}'.format(mbid) + '/bg/PREP_' + date + '.nc'
            self.sodaName = options.xpiddir + '/mb{0:04d}'.format(mbid) + '/bg/PREP_' + date + '.nc'


class PrepOl(Prep):
    def __init__(self, date, mbid, options, directFromXp = True, isOl = False):
        Prep.__init__(self, options)
        self.date = date
        self.ptinom = 'ol' + str(mbid)
        if not directFromXp:
            self.sodaName = date + '/PREP_' + convertdate(date).strftime('%y%m%dH%H') + '_PF_ENS' + str(mbid) + '.nc'
        else:
            if isOl:
                self.sodaName = options.xpiddir + '/mb{0:04d}'.format(mbid) + '/bg/PREP_' + date + '.nc'
            else:
                self.sodaName = options.xpidoldir + '/mb{0:04d}'.format(mbid) + '/bg/PREP_' + date + '.nc'


class PrepAn(Prep):
    def __init__(self, date, mbid, options, directFromXp = True):
        Prep.__init__(self, options)
        self.date = date
        self.ptinom = 'an' + str(mbid)

        if not directFromXp:
            self.sodaName = date + '/SURFOUT' + str(mbid) + '.nc'
        else:
            self.sodaName = options.xpiddir + '/mb{0:04d}'.format(mbid) + '/an/PREP_' + date + '.nc'


class PrepAbs(Prep):
    '''
    class describing abstract prep-like objects ( computed from prep, not bound to a specific file)
    e.g. medians, covariance etc.
    '''

    def __init__(self, date, options, ptinom):
        Prep.__init__(self, options)
        self.date = date
        self.ptinom = ptinom
        self.data = dict()
