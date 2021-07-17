#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 8 févr. 2020

@author: cluzetb

on beaufix, generate  conf file with the appropriate assimilation dates (for ex. every 7 days)

'''
import datetime
import os

import netCDF4  # @UnresolvedImport

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


def set_legend_mountain(pp, id_postes, ax, width, line):
    # set the legend
    # generate proxy scatter for legend entries:
    s = []
    l = []
    for mountain in list(set(mountain_from_massif(pp.pgd.massif[id_postes]))):
        gg = plt.scatter([], [],
                         color=mountain_color([mountain]),
                         )
        s.append(gg)
        l.append(mountain)

    # generate a composite legend entry for line and envelope
    # proxy for envelope:
    p2 = ax.fill(np.NaN, np.NaN, color='C3', alpha=0.5)
    s.append((p2[0], line[0]))
    l.append(
        '{0}-posts-width\nrolling median $\pm 1\sigma$'.format(width))
    ax.legend(s, l, loc='lower_right')
    return 0


def set_itimes_posts(run):
    """
    set masks to put in common an ensemble (an or ol), its obs timeseries and the associated oper run.
    """
    if run.options.openloop == 'on':
        timesEns = run.ensProOl['time']
    else:
        timesEns = run.ensProAn['time']
    year = timesEns[0].year
    mobs = [i for i, t in enumerate(run.obsTs['time'])
            if (t.hour == 6 and ((t.year == year and t.month > 9) or (t.month < 7 and t.year == year + 1)))]
    mens = [i for i, t in enumerate(run.ensProOl['time']) if (
        t.hour == 6 and (t.month > 9 or t.month < 7))]
    moper = [i for i, t in enumerate(run.oper['time']) if (
        t.hour == 6 and (t.month > 9 or t.month < 7))]

    return year, mobs, mens, moper

#
# def mountain_color(mountains):
#     """
#     return a color for each mountain.
#     """
#     cols = []
#     if isinstance(mountains, str):
#         mountains = [mountains]
#     for mountain in mountains:
#         if 'alpes' in mountain:
#             cols.append('C0')
#         elif 'pyrenees' in mountain:
#             cols.append('C1')
#         elif 'haute-ariege_andorre' in mountain:
#             cols.append('C4')
#         elif 'corse' in mountain:
#             cols.append('C2')
#         elif 'other' in mountain:
#             cols.append('C5')
#         else:
#             print(mountain)
#             print(type(mountain))
#             raise Exception(' no color for this mountain at the moment sorry')
#     return cols


def mountain_color(mountains, highlight=['alpes', 'pyrenees']):
    """
    return a color for each mountain.
    """
    dict_cols = {'alpes': 'C0',
                 'pyrenees': 'C1',
                 'corse': 'C2',
                 'beaufortain': 'C7',
                 'haute-ariege_andorre': 'r',
                 'other': 'C5',
                 }
    dcols = {h: dict_cols[h] for h in highlight}
    if isinstance(mountains, str):
        mountains = [mountains]
    cols = [dcols[m] if m in dcols.keys() else 'C0' for m in mountains ]
    return cols


def mountain_from_massif(massifs):
    """
    Determine the mountain (alpes, pyrenees or corse) from the massif number
    in: list
    out:list
    """
    mountains = []
    if isinstance(massifs, np.int32):
        massifs = [massifs]
    for massif in massifs:
        #         if massif in [5]:
        #             mountains.append('beaufortain')
        #         elif massif in [17]:
        #             mountains.append('haute-maurienne')
        if massif <= 23:
            mountains.append('alpes')
        elif massif in [40, 41]:
            mountains.append('corse')
        # elif massif in [70, 71]:
        #    mountains.append('haute-ariege_andorre')
        elif massif > 91:
            print(' No true massifs with num {0} exists'.format(massif))
            mountains.append('other')
        else:
            mountains.append('pyrenees')
    return mountains
