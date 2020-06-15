#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 8 févr. 2020

@author: cluzetb

on beaufix, generate  conf file with the appropriate assimilation dates (for ex. every 7 days)

'''
import datetime
import os

import netCDF4

import numpy as np
import pandas as pd


def set_conf_everydate(startY, assimilate_every, filepath, nens=40, endY = None):
    '''
    BC 08/02/20
    on beaufix, generate  conf file with the appropriate assimilation dates (for ex. every 7 days)

    '''
    if endY is None:
        endY = startY + 1
    # assimilate between october and june
    for i, year in enumerate(range(startY, endY)):
        startDate = datetime.datetime(year, 10, 1, 6, 0, 0)
        endDate = datetime.datetime(year + 1, 6, 30, 6, 0, 0)
        datelist = pd.date_range(startDate, periods=(endDate - startDate).days).tolist()
        datelist = [d.strftime('%Y%m%d%H') for d in datelist[::assimilate_every]]
        if i == 0:
            stringdates = ','.join(datelist)
        else:
            stringdates = ','.join([stringdates] + datelist)

    # members Id copied from the previous experiments (always the same for reproducibility
    members_id = [168, 429, 524, 1460, 1568, 698, 1084, 162, 1159, 872, 466, 1646, 1550, 1582, 403, 165, 1, 540,
                  632, 1309, 1665, 1863, 1787, 270, 79, 1188, 385, 1782, 1271, 1380, 864, 880, 971, 54, 609, 289, 1279, 415,
                  1754, 612, 1431, 696, 1215, 1861, 313, 734, 1813, 1710, 567, 1161, 1792, 653, 326, 73, 182, 1428, 453, 628,
                  78, 906, 720, 527, 1212, 1717, 386, 655, 1909, 1352, 1100, 151, 515, 595, 784, 1284, 477, 675, 110, 99, 1052,
                  444, 397, 1697, 284, 1699, 1394, 1585, 715, 193, 848, 1324, 59, 1015, 687, 255, 1491, 7, 1235, 879, 18, 1744,
                  929, 1597, 1183, 384, 1157, 1223, 39, 1641, 1622, 1644, 1648, 1446, 1709, 1653, 888, 1095, 214, 483, 1337,
                  1075, 455, 1009, 1742, 95, 1840, 250, 1407, 189, 1928, 1099, 916, 412, 361, 112, 423, 194, 843, 1857, 439,
                  1294, 71, 1872, 451, 332, 322, 543, 23, 1008, 1856, 1371, 1822, 1884, 1409, 788, 1011, 1277, 1573, 665, 390,
                  77
                  ]
    members_id = list(map(str, members_id))

    if not os.path.exists(filepath):
        try:
            os.makedirs(os.path.dirname(filepath))
        except FileExistsError:
            pass
    else:
        os.remove(filepath)
    f = open(filepath, 'w')
    f.write('[DEFAULT]')
    f.write('\n')
    f.write('assimdates=' + stringdates)
    f.write('\n')
    f.write('nmembers=' + str(nens))
    f.write('\n')
    f.write('members_id=' + ','.join(members_id[0:nens]))
    f.write('\n')
    f.close()
    return 0


class Mnt(object):
    """
    class to read a semi-distributed Mnt file
    slope is in degrees
    """

    def __init__(self, pathMnt):
        mnt = netCDF4.Dataset(pathMnt)
        self.elev = np.squeeze(mnt.variables['Elevation'][:])  # hyp a verif : elev est le mn de la classe safran corr
        self.slope = np.squeeze(mnt.variables['Slope'][:])
        self.asp = np.squeeze(mnt.variables['Aspect'][:])
        self.massif = np.squeeze(mnt.variables['Massif'][:])
        self.ridges = np.squeeze(mnt.variables['Ridges'][:])
        self.landcov = np.squeeze(mnt.variables['Landcover'][:]).astype(int)
        self.x = np.squeeze(mnt.variables['x'][:])
        self.y = np.squeeze(mnt.variables['y'][:])
        mnt.close()
        print('mnt loaded')
