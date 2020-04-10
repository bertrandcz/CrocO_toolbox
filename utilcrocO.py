# -*- coding: utf-8 -*-
'''
Created on 6 févr. 2019

@author: cluzetb

utils suited for crocO interface only
'''

from bronx.datagrip.namelist import NamelistParser
from configparser import ConfigParser
import datetime
from ftplib import FTP
from netrc import netrc
import os
import re
import shutil
from vortex.layout.nodes import ConfigSet
from vortex.util.config import GenericConfigParser

from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4

import numpy as np
import pandas as pd


def dictsAspect():
    '''
    returns the dict for aspect and its reverse.
    '''

    gg1 = {'N': 0, 'NE': 45, 'E': 90, 'SE': 135, 'S': 180, 'SW': 225, 'W': 270, 'NW': 315, 'flat': -1}

    return gg1, {v: k for k, v in list(gg1.items())}


def setSubsetclasses(pgd, selE, selA, selS):
    """
    BC 5/02/19
    Returns a list of point ids and a mask corresponding to the selection of topographic classes (selE,selA, selS
    params:
    - selE : string or list of string for elevations (['1800,'2100', ...]. 'all' for all elevation bands
    - selA : string or list of strings for aspects (['N, 'NW', ...]), 'all' for all.
    - selS : string or list of strings for slopes (degrees) (['0,20',]), 'all' for all. 
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
        classesE = list(map(int, selE))
    else:
        classesE = dictElev['all']

    classesA = []
    if 'all' not in selA:
        for cl, asp in enumerate(selA):
            classesA.append(dictAsp[asp])
    else:  # avoid having a list of list
        classesA = dictAsp['all']

    classesS = []
    if isinstance(selS, str):
        selS = [selS]
    if 'all' not in selS:
        classesS = selS
    else:  # avoid having a list of list
        classesS = dictSlope['all']
    mask = []
    for cl in range(pgd.npts):
        if pgd.elev[cl] in classesE and (
                (str(int(np.arctan(pgd.slope[cl]) * 180. / np.pi)) in classesS and pgd.aspect[cl] in classesA) or
                (pgd.slope[cl] < 0.01 and ('0' in classesS))):

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
            'DEP': 'DEP_TOT', 'SWE': 'SWE_TOT'}


def dictvarsPro():
    return {'B1': 'SPM_VEG1', 'B2': 'SPM_VEG2', 'B3': 'SPM_VEG3',
            'B4': 'SPM_VEG4', 'B5': 'SPM_VEG5', 'B6': 'SPM_VEG6',
            'B7': 'SPM_VEG7',
            'SCF': 'WSN_T_ISBA',  # SCF is not in the pro, fake it reading WSN_T_ISBA
            'R53': 'R53', 'R52': 'R52', 'R51': 'R51', 'R54': 'R54', 'R21': 'R21', 'R23': 'R23', 'R24': 'R24',
            'DEP': 'DSN_T_ISBA', 'SWE': 'WSN_T_ISBA'}


def niceName(pgd, cl, tolist = False):
    _, revdictAsp = dictsAspect()
    return str(int(pgd.elev[cl])) + '_' + revdictAsp[pgd.aspect[cl]] + '_' + str(int(np.arctan(pgd.slope[cl]) * 180. / np.pi))


def niceLabel(var, score = None, printunits=True):

    if score is None:
        sc = ''
    else:
        sc = score
    units = {'SWE': r'[$\mathrm{\mathsf{kgm^{-2}}}$]',
             'DEP': r'[$\mathrm{\mathsf{m}}$]',
             'B5': '',
             'B4': '',
             }
    if printunits:
        u = units[var]
    else:
        u = ''

    ddict = {'SWE': 'SWE {0} {1}'.format(sc, u),
             'DEP': 'HS {0} {1}'.format(sc, u),
             'B5': 'Band 5 {0} {1}'.format(sc, u),
             'B4': 'Band 4 {0} {1}'.format(sc, u),
             }
    return ddict[var]


def cm2inch(w, h):
    return(0.393701 * w, 0.393701 * h)


class Pgd(object):
    """
    class to read a semi-distributed PGD file
    slope is the TANGENT of the angle of slope.
    """

    def __init__(self, pathPGD):

        pgd = netCDF4.Dataset(pathPGD)
        self.path = pathPGD
        self.elev = np.squeeze(pgd.variables['MIN_ZS'][:])  # lower altitude
        self.slope = np.squeeze(pgd.variables['SSO_SLOPE'][:])
        self.aspect = np.squeeze(pgd.variables['SSO_DIR'][:])
        self.npts = self.elev.size
        self.lat = np.squeeze(pgd.variables['XY'][:])
        self.lon = np.squeeze(pgd.variables['XX'][:])
        # in the case of postes geometries, a "SUPER" PGD (enhanced with numposts ('station') and type of postes) is put at pathPGD
        if 'station' in pgd.variables.keys():
            self.station = np.squeeze(pgd.variables['station'][:])
            self.type = np.squeeze(pgd.variables['type'][:])
            self.massif = np.squeeze(pgd.variables['massif'][:])
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
    convert a crocO argument options.vars into a list of OBS variables names in soda format
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
    convert a crocO argument options.vars into a list of VAR variables names in soda format
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


def unpack_conf(arg):
    """
    BC 01/04/20
    supported types:
    arg = aa,bb,cc : return ['aa','bb','cc']
    arg = 'aa','bb','cc' : return ['aa','bb','cc']
    arg = "'toto'" or '"toto"' : return 'toto'
    arg = '2' : return 2
    arg = '2.3' : return 2.3
    arg = 'toto.nam' :'return toto.nam'
    arg = '2.2,2.3' : return [2.2,2.3]
    arg = None : return None (NoneType)

    /!| tricky case:
    arg='2016080106' (YYYYmmddHH format) : return '2016080106' (and not 2016080106)
    arg=rangex(start: a end: b) : return list(range(a,b+1)) (a vortex type...)
    arg=rangex(start:a end:b): return list(range(a,b+1))
    exception cases :
    - badly formatted str
    - mixed types in lists
    """
    # first deal with the vortex case
    if "rangex" in arg:
        try:
            start = int(re.search(r"start: (\d+)", arg).group(1))
        except AttributeError:
            start = int(re.search(r"start:(\d+)", arg).group(1))
        try:
            end = int(re.search(r"end: (\d+)", arg).group(1))
        except AttributeError:
            end = int(re.search(r"end:(\d+)", arg).group(1))
        ret = list(range(start, end + 1))
    elif arg == 'None':
        ret = None
    else:
        if "," in arg:
            ret = list(map(unpack_conf, arg.split(',')))
        # unpack strings:
        elif arg.startswith("'") or arg.startswith('"'):
            if arg[0] == arg[-1]:
                ret = arg[1:-1]
            else:
                raise ValueError(arg + ': badly formatted string.')
        elif '.' in arg:
            try:
                ret = float(arg)
            except ValueError:  # must be a path ^^
                ret = arg
        else:
            try:
                _ = datetime.datetime.strptime(arg, '%Y%m%d%H')
                ret = arg
            except ValueError:
                try:
                    ret = int(arg)
                except ValueError:
                    ret = arg
    if isinstance(ret, list):
        if all(isinstance(r, type(ret[0])) for r in ret):
            return ret
        else:
            raise TypeError('unconsistent types in list: {0}'.format(ret))
    else:
        return ret


class ConfObj1L(dict):
    def __init__(self, **kwargs):
        self.__dict__.update(self, **kwargs)
        for name, val in kwargs.items():
            super().__setitem__(name, val)

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, val):
        self[name] = val
        super().__setitem__(name, val)  # useful for object description
        self.__dict__[name] = val  # useful for iteration over self.__dict__

    def __str__(self):
        return '<ConfObj1L(' + ', '.join(['{0}:{1}'.format(key, val) for key, val in self.__dict__.items()]) + ')>'

    def __repr__(self):
        return '<ConfObj1L(' + ', '.join(['{0}:{1}'.format(key, val) for key, val in self.__dict__.items()]) + ')>'


def conf2obj(conf):
    '''
    BC 01/04/20
    convert a ConfParser to a one-level dict then a one level dot.dict.
    '''

    dict1 = dict()
    default_sec = conf.default_section
    for k in conf[default_sec]:
        dict1[k] = unpack_conf(conf[default_sec][k])
    for k in conf._sections:
        for kk in conf[k]:
            dict1[kk] = unpack_conf(conf[k][kk])
    return ConfObj1L(**dict1)


def read_conf(pathconf, useVortex=True):
    '''
    B. Cluzet
    duplicated from evalSODA.util
    '''
    if not os.path.exists(pathconf):
        if os.path.exists(pathconf[0:-4] + '.foo'):
            shutil.copyfile(pathconf[0:-4] + '.foo', pathconf)
        else:
            print('no conf file for this experiment :', pathconf)
    if useVortex is True:
        iniparser = GenericConfigParser(pathconf)
        thisconf  = iniparser.as_dict(merged=False)
        updconf = thisconf.get('defaults', dict())
        conf = ConfigSet()
        conf.update(updconf)
    else:
        conf = ConfigParser()
        conf.read(pathconf)
        conf = conf2obj(conf)
    return conf


def colorbar(mappable, ax = None):
    """
    from http://joseph-long.com/writing/colorbars/
    """
    if ax is None:
        ax = mappable.axes
    else:
        ax = ax
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


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


def set_conf_everydate(startY, assimilate_every, filepath, nens=40, endY = None):
    '''
    BC 08/02/20
    on beaufix, generate  conf file with the appropriate assimilation dates (for ex. every 7 days)

    '''
    if endY is None:
        endY = startY + 1

    # assimilate between october and june
    startDate = datetime.datetime(startY, 10, 1, 6, 0, 0)
    endDate = datetime.datetime(endY, 6, 30, 6, 0, 0)

    datelist = pd.date_range(startDate, periods=(endDate - startDate).days).tolist()
    stringdates = ','.join([d.strftime('%Y%m%d%H')for d in datelist[::assimilate_every]])

    # members Id copied from the previous experiments (always the same for reproducibility)
    membersId = [168, 429, 524, 1460, 1568, 698, 1084, 162, 1159, 872, 466, 1646, 1550, 1582, 403, 165, 1, 540,
                 632, 1309, 1665, 1863, 1787, 270, 79, 1188, 385, 1782, 1271, 1380, 864, 880, 971, 54, 609, 289, 1279, 415,
                 1754, 612, 1431, 696, 1215, 1861, 313, 734, 1813, 1710, 567, 1161, 1792, 653, 326, 73, 182, 1428, 453, 628,
                 78, 906, 720, 527, 1212, 1717, 386, 655, 1909, 1352, 1100, 151, 515, 595, 784, 1284, 477, 675, 110, 99, 1052,
                 444, 397, 1697, 284, 1699, 1394, 1585, 715, 193, 848, 1324, 59, 1015, 687, 255, 1491, 7, 1235, 879, 18, 1744,
                 929, 1597, 1183, 384, 1157, 1223, 39, 1641, 1622, 1644, 1648, 1446, 1709, 1653, 888, 1095, 214, 483, 1337,
                 1075, 455, 1009, 1742, 95, 1840, 250, 1407, 189, 1928, 1099, 916, 412, 361, 112, 423, 194, 843, 1857, 439,
                 1294, 71, 1872, 451, 332, 322, 543, 23, 1008, 1856, 1371, 1822, 1884, 1409, 788, 1011, 1277, 1573, 665, 390,
                 77
                 ]
    membersId = list(map(str, membersId))

    # if not os.path.exists(filepath):
    f = open(filepath, 'w')
    f.write('[DEFAULT]')
    f.write('\n')
    f.write('assimdates=' + stringdates)
    f.write('\n')
    f.write('nmembers=' + str(nens))
    f.write('\n')
    f.write('membersId=' + ','.join(membersId[0:nens]))
    f.write('\n')
    f.close()
    return 0
