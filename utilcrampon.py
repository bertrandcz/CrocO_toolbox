# -*- coding: utf-8 -*-
'''
Created on 6 févr. 2019

@author: cluzetb

utils suited for crampon interface only
'''

from bronx.datagrip.namelist import NamelistParser
import datetime
from ftplib import FTP
from netrc import netrc
from optparse import Values
import os
import re
import shutil

import netCDF4

import numpy as np


def dictsAspect():
    '''
    returns the dict for aspect and its reverse.
    '''

    gg1 = {'N': 0, 'NE': 45, 'E': 90, 'SE': 135, 'S': 180, 'SW': 225, 'W': 270, 'NW': 315, 'flat': -1}

    return gg1, {v: k for k, v in list(gg1.items())}


def parse_classes(options):
    """
    BC, april 2020
    This function parses option arguments classes_e, classes_a, and classes_s into explicit conditions matching the options.pgd features
    classes_e = 'all' -> classes_e = ['600',...'3600']
    classes_a = 'all' -> classes_a = ['N','NW',...]
    """

    # check for conf file (not necessarily defined in old conf files
    for attr in ['classes_e', 'classes_a', 'classes_s']:
        if not hasattr(options, attr):
            options.__setattr__(attr, None)

    if options.classes_e == 'all' or options.classes_e is None:
        options.classes_e = sorted(np.unique(list(map(str, map(int, options.pgd.elev)))))
    if options.classes_a == 'all' or options.classes_a is None:
        options.classes_a = sorted(list(dictsAspect()[0].keys()))
    if options.classes_s == 'all' or options.classes_s is None:
        options.classes_s = ['0', '20', '40']
    return options


def setSubsetclasses(pgd, selE, selA, selS):
    """
    BC 5/02/19
    Returns a list and a mask corresponding to points whose topographic conditions match any combination of the selection (selE,selA, selS)
    params:
    - selE : string or list of string for elevations (['1800','2100', ...]. 'all' for all elevation bands
    - selA : string or list of strings for aspects (['N', 'NW', ...]), 'all' for all.
    - selS : string or list of strings for slopes (degrees) (['0','20',]), 'all' for all.
       """
    subsetClass = []
    dictElev = {'all': np.unique(pgd.elev)}
    dictAsp, _ = dictsAspect()
    dictAsp['all'] = np.unique(pgd.aspect)
    dictSlope = {'all': ['0', '20', '40']}

    if isinstance(selE, np.ndarray):
        selE = np.ndarray.tolist(selE)
    elif isinstance(selE, str):
        selE = [selE]
    if 'all' not in selE:
        classes_e = list(map(int, selE))
    else:
        classes_e = dictElev['all']

    classes_a = []
    if 'all' not in selA:
        for cl, asp in enumerate(selA):
            classes_a.append(dictAsp[asp])
    else:  # avoid having a list of list
        classes_a = dictAsp['all']

    classes_s = []
    if isinstance(selS, str):
        selS = [selS]
    if 'all' not in selS:
        classes_s = selS
    else:  # avoid having a list of list
        classes_s = dictSlope['all']
    mask = []
    for cl in range(pgd.npts):
        if pgd.elev[cl] in classes_e and (
                (str(int(np.arctan(pgd.slope[cl]) * 180. / np.pi)) in classes_s and pgd.aspect[cl] in classes_a) or
                (pgd.slope[cl] < 0.01 and ('0' in classes_s))):

            subsetClass.append(cl)
            mask.append(True)
        else:
            mask.append(False)

    return subsetClass, np.array(mask)


def dictvarsPrep():
    return {'B1': 'SPM_VEG1', 'B2': 'SPM_VEG2', 'B3': 'SPM_VEG3',
            'B4': 'SPM_VEG4', 'B5': 'SPM_VEG5', 'B6': 'SPM_VEG6',
            'B7': 'SPM_VEG7',
            'SCF': 'WSN_VEG1',  # computing of SCF requires to read SWE_tot
            'R53': 'R53', 'R52': 'R52', 'R51': 'R51', 'R54': 'R54', 'R21': 'R21', 'R23': 'R23', 'R24': 'R24',
            'DEP': 'DEP', 'SWE': 'SWE'}  # caution : DEP_TOT is the depth of the previous timestep !!


def dictvarsPro():
    return {'B1': 'SPM_VEG1', 'B2': 'SPM_VEG2', 'B3': 'SPM_VEG3',
            'B4': 'SPM_VEG4', 'B5': 'SPM_VEG5', 'B6': 'SPM_VEG6',
            'B7': 'SPM_VEG7',
            'SCF': 'WSN_T_ISBA',  # SCF is not in the pro, fake it reading WSN_T_ISBA
            'R53': 'R53', 'R52': 'R52', 'R51': 'R51', 'R54': 'R54', 'R21': 'R21', 'R23': 'R23', 'R24': 'R24',
            'DEP': 'DSN_T_ISBA', 'SWE': 'WSN_T_ISBA'}


class Opt(Values):
    """
    Opt is an object representing options of a command.
    This class is initialized with a dict.
    """

    def __init__(self, **kwargs):
        for name, val in kwargs.items():
            # pydev editor show an error for Values.__setattr__
            super(Opt, self).__setattr__(name, val)


class ImmutableOpt(Opt):
    """
    This class is initialized with a dict.
    Once instanciated, the object is immutable
    """

    def __setattr__(self, name, value):
        raise Exception('This class is immutable.')


class Pgd:
    """
    class to read a semi-distributed PGD file
    slope is the TANGENT of the angle of slope.
    """

    def __init__(self, pathPGD):

        pgd = netCDF4.Dataset(pathPGD)
        bugfix = pgd.variables['BUG']
        if bugfix == 0:
            raise Exception(' Version of your PGD is deprecated (BUG==0) please update it by rerunning a spinup so that BUG>=1')
        self.path = pathPGD
        self.elev = np.squeeze(pgd.variables['MIN_ZS'][:])  # lower altitude
        self.slope = np.squeeze(pgd.variables['SSO_SLOPE'][:])
        self.aspect = np.squeeze(pgd.variables['SSO_DIR'][:])
        self.npts = self.elev.size
        self.lat = np.squeeze(pgd.variables['XY'][:])
        self.lon = np.squeeze(pgd.variables['XX'][:])

        pgd.close()


def convertdate(date):
    '''
    YYYYMMDDHH to datetime.datetime
    '''
    return datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]), int(date[8:10]), 0, 0)


'''
def colorbar(mappable):
    """
    from http://joseph-long.com/writing/colorbars/
    """
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)
'''


def setlistvars_obs(arg):
    """
    BC 6/02/19
    convert a crampon argument options.vars into a list of OBS variables names in soda format
    """

    if arg == 'all':
        listvar = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'DEP']
    else:
        listvar = []
        for var in arg:
            listvar.append(var)
    return listvar


def setlistvars_var(arg):
    """
    BC 6/02/19
    convert a crampon argument options.vars into a list of VAR variables names in soda format
    TODO : same stuff for DEP/ SCF etc.
    """
    if arg == 'all':
        listvar = ['PB1', 'PB2', 'PB3', 'PB4', 'PB5', 'PB6', 'PB7', 'DEP']
    else:
        listvar = []
        for var in arg:
            if (var == 'DEP'):
                listvar.append('DEP')
            elif (var == 'SCF'):
                listvar.append('PSB')
            elif 'R' in var:
                listvar.append(var)
            else:
                listvar.append('P' + var)  # 'b*' -> 'PB*'

    return listvar


def area(geometry):
    if geometry == '13':
        area = 'thabor'
    elif geometry == '12':
        area = 'grandes_rousses'
    else:
        area = geometry

    return area


def read_conf(pathconf):
    '''
    B. Cluzet
    duplicated from evalSODA.util
    '''
    if not os.path.exists(pathconf):
        if os.path.exists(pathconf[0:-4] + '.foo'):
            shutil.copyfile(pathconf[0:-4] + '.foo', pathconf)
        else:
            raise FileNotFoundError('no conf file at this path :', pathconf)
    try:
        from vortex.layout.nodes import ConfigSet
        from vortex.util.config import GenericConfigParser
    except ImportError:
        raise Exception("vortex is necessary to parse conf files.")
    iniparser = GenericConfigParser(pathconf)
    thisconf  = iniparser.as_dict(merged=False)
    updconf = thisconf.get('defaults', dict())
    conf = ConfigSet()
    conf.update(updconf)

    return conf


def dictErrors():
    return {'B1': 0.00071, 'B2': 0.00046, 'B3': 0.00056, 'B4': 0.00056, 'B5': 0.002,
            'B6': 0.0015, 'B7': 0.00078, 'SCF': 0.2, 'DEP': 0.01, 'SWE': 100, 'R52': 0.001, 'R54': 0.001}


def set_errors(argsoda):
    '''
    BC 06/02/19
    set soda canonical errors for the prescribed vars (Wright et al., Charrois et al.)
    argsoda is a list
    '''

    dicterrors = dictErrors()
    ret = []
    for el in argsoda:
        ret.append(dicterrors[el])
    return ret


def set_factors(argsoda, fact):
    '''
    BC 06/02/19
    properly set the error factors for namelist Writing
    fact (list of string is prescribed in namelist :
    - if only one value (default or lazy case), apply it to all variables
    - if list of values (exhaustive): ok
    '''

    if type(fact) is list:
        return list(map(float, fact))
    if len(argsoda) == len([fact]):  # default with 1 var, lazy with 1 var, exhaustive
        if len([fact]) == 1:
            return [fact]  # bc probably bugged
        else:
            return fact  # idem probably bugged
    else:
        if len([fact]) == 1:  # lazy/default
            return list(map(float, [fact])) * len(argsoda)
        else:
            raise Exception('you should either apply the same fact to all vars or specify it for each one')


def read_opts_in_namelist(options):
    """
    read pf parameters in namelist
    """
    n = NamelistParser()
    if hasattr(options, 'kind'):
        if options.todo != 'localpp':
            try:
                N = n.parse(options.xpiddir + 'conf/namelist.surfex.foo')
                print(' read the PF params in namelist:',
                      options.xpiddir + 'conf/namelist.surfex.foo')
            except FileNotFoundError:
                N = n.parse(options.xpiddir + '/conf/OPTIONS.nam')
                print(' read the PF params in namelist:', options.xpiddir + '/conf/OPTIONS.nam')
    else:
        N = n.parse(options.dates[0] + '/OPTIONS.nam')
        print(' read the PF params in namelist:', options.dates[0] + '/OPTIONS.nam')
    try:
        options.nloc_pf = N['NAM_ASSIM'].NLOC_PF
        options.neff_pf = N['NAM_ASSIM'].NEFF_PF
        options.pf_crocus = N['NAM_ASSIM'].CPF_CROCUS
    except AttributeError:
        raise Exception('Some of the PF parameters are not defined in the namelist.')


def check_namelist_soda(options, pathIn= None, pathOut = None):
    """
    - check consistency between the prescribed assim vars and what is written in the namelist
    - check and change LWRITE_TOPO to .FALSE.
    """
    n = NamelistParser()
    if pathIn is None:
        N = n.parse('OPTIONS_base.nam')
    else:
        N = n.parse(pathIn)
    # LWRITE_TOPO must be false if we want ISBA_ANALYSIS.out to write
    print('force LWRITE_TOPO to .FALSE.')
    N['NAM_IO_OFFLINE'].LWRITE_TOPO = False

    # update assim vars in the namelist
    # NAM_ENS
    N['NAM_ENS'].NENS = options.nmembers

    # NAM_OBS
    sodaobs = setlistvars_obs(options.vars)
    if 'SCF' in sodaobs:
        sodaobs.remove('SCF')
    if 'NAM_OBS' not in list(N.keys()):
        N.newblock('NAM_OBS')
    if 'NAM_VAR' not in list(N.keys()):
        N.newblock('NAM_VAR')
    if 'NAM_ASSIM' not in list(N.keys()):
        N.newblock('NAM_ASSIM')
    N['NAM_OBS'].NOBSTYPE = len(sodaobs)
    N['NAM_OBS'].COBS_M = sodaobs
    N['NAM_OBS'].XERROBS_M = set_errors(sodaobs)
    N['NAM_OBS'].XERROBS_FACTOR_M = set_factors(sodaobs, options.fact)
    N['NAM_OBS'].NNCO = [1] * len(sodaobs)
    N['NAM_OBS'].CFILE_FORMAT_OBS = "NC    "

    # NAM_VAR
    sodavar = setlistvars_var(options.vars)
    print(sodavar)
    if 'PSB' in sodavar:
        sodavar.remove('PSB')
    N['NAM_VAR'].NVAR = N['NAM_OBS'].NOBSTYPE
    N['NAM_VAR'].CVAR_M = sodavar
    N['NAM_VAR'].NNCV = N['NAM_OBS'].NNCO

    # NAM_ASSIM
    if hasattr(N['NAM_ASSIM'], 'LSEMIDISTR_CROCUS') or hasattr(N['NAM_ASSIM'], 'LASSIM_CROCUS'):
        print('be careful, old-formatted namelist !', 'LSEMIDISTR_CROCUS' 'LASSIM_CROCUS -> CPF_CROCUS, LCRAMPON')
        N['NAM_ASSIM'].delvar('LSEMIDISTR_CROCUS')
        N['NAM_ASSIM'].delvar('LASSIM_CROCUS')
    N['NAM_ASSIM'].CPF_CROCUS = options.pf.upper()
    N['NAM_ASSIM'].LCRAMPON = True
    N['NAM_ASSIM'].LASSIM = True
    N['NAM_ASSIM'].CASSIM_ISBA = 'PF   '
    N['NAM_ASSIM'].NLOC_PF = options.nloc_pf
    N['NAM_ASSIM'].NEFF_PF = options.neff
    N['NAM_ASSIM'].LEXTRAP_SEA = False
    N['NAM_ASSIM'].LEXTRAP_WATER = False
    N['NAM_ASSIM'].LEXTRAP_NATURE = False
    N['NAM_ASSIM'].LWATERTG2 = True

    namSURFEX = open(pathOut if pathOut is not None else 'OPTIONS.nam', 'w')
    namSURFEX.write(N.dumps())
    namSURFEX.close()


def todates(dataset):
    """
    inspired on dbazin (forum.marine.copernicus.eu)
    returns the netCDF4 variable and the corresponding dates.
    """
    dtime = dataset.variables['time'][:]
    dtime_units = dataset.variables['time'].units
    try:
        dtime_cal = dataset.variables['time'].calendar
    except AttributeError:
        dtime_cal = u'gregorian'
    return dtime, netCDF4.num2date(dtime, units = dtime_units, calendar = dtime_cal)


def ftpconnect(machine):
    """
    copy from Carlo Carmagnola PEARP_get.py
    """

    username, _, password = netrc().authenticators(machine)
    ftp = FTP(machine)
    ftp.set_debuglevel(1)
    ftp.login(username, password)
    return ftp


def get_trailing_number(s):
    m = re.search(r'\d+$', s)
    return str(int(m.group())) if m else None


def get_leading_number(s):
    m = re.search(r'^\d+', s)
    return str(int(m.group())) if m else None


def ftp_upload(localfile, remotefile, ftp):
    '''
    source:
    http://makble.com/upload-new-files-to-ftp-server-with-python
    '''

    fp = open(localfile, 'rb')
    try:
        ftp.storbinary('STOR %s' % remotefile, fp, 1024)
    except Exception:
        print("remotefile not exist error caught" + remotefile)
        path, _ = os.path.split(remotefile)
        print("creating directory: " + remotefile)
        ftp.mkd(path)
        ftp_upload(localfile, remotefile, ftp)
        fp.close()
        return
    fp.close()
    print("after upload " + localfile + " to " + remotefile)
