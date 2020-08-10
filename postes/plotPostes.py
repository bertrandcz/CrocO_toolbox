'''
Created on 30 mai 2020

@author: cluzetb
Advances plotting facilities suited to the postes geometry
'''
import os

from pyproj import Proj, transform

import matplotlib.pyplot as plt
import numpy as np
from plotcrocO import Pie
import shapefile as shp  # @UnresolvedImport
from crocO import callunits


class MapPostes(Pie):
    '''
    Inspired on Pie, the principle is to plot maps of a variable (e.g. correlation or score)
    in the postes geopetry: the postes denoted by points are colored according to the variable value.
    The plotting map facility comes from notebooks.postes.first_glance*.ipynb
    '''

    def __init__(self, *args, **kwargs):
        Pie.__init__(self, *args, **kwargs)
        # load shapefile of massifs from the alps and pyrenees into a same shp.
        # (https://gis.stackexchange.com/questions/103033/using-pyshp-to-merge-two-shapefiles)

        self.setMassifShapes()

    def setMassifShapes(self):
        self.massifs = dict()
        for massif in np.unique(self.sdObj.pgd.massif):
            if massif < 91:
                self.massifs[massif] = MassifReal(massif)
            else:
                # BC 30/05/20
                # massifs >91 are special cases of locations outside massifs borders.
                # Treat it as single-point massifs ? group them in a single massif ?
                # @TODO: something smart.
                pass

    def plot(self, tag_postes = False, tag_massifs = True, plotCircle = True, vmin = -1, vmax = 1, cmap = 'RdBu'):
        fig, ax = plt.subplots()

        for mnum, massif in self.massifs.items():
            ax.plot(massif.x, massif.y, color = 'k', zorder = 0)
            if tag_massifs:
                ax.annotate(massif.name, massif.center)
        ax.scatter(self.sdObj.pgd.lon, self.sdObj.pgd.lat, c = self.sdObj.data['DEP'],
                   cmap = cmap, vmin = vmin, vmax = vmax, zorder = 100)
        if tag_postes:
            for ipt in range(self.sdObj.pgd.npts):
                ax.annotate(str(self.sdObj.pgd.station[ipt]), (self.sdObj.pgd.lon[ipt], self.sdObj.pgd.lat[ipt]),)

        if hasattr(self, 'focusCl'):
            ax.scatter(self.sdObj.pgd.lon[self.focusCl], self.sdObj.pgd.lat[self.focusCl],
                       c= 'r', zorder = 1000)
            # stars on the 11 highest values
            #  need to exclude the nans
            nNan = len(self.sdObj.data['DEP']) - self.sdObj.data['DEP'].count()
            print('number of nans', nNan)
            bestargs = np.ma.argsort(self.sdObj.data['DEP'])[(-nNan - 10):(- nNan)][::-1]
            _ = ax.scatter(self.sdObj.pgd.lon[bestargs], self.sdObj.pgd.lat[bestargs],
                           c= 'g', zorder = 1000, marker = '*')
        if 'rlocal' in self.sdObj.ptinom and plotCircle:
            xmin = ax.get_xlim()[0]
            ymax = ax.get_ylim()[1]
            radius = callunits(self.sdObj.ptinom.split(',')[1].split(' ')[0])
            circle = plt.Circle((xmin + radius, ymax - radius), radius, color = 'g', fill = False)
            ax.add_artist(circle)
        ax.set_aspect('equal')
        plt.title(self.sdObj.ptinom)


class Massif(object):
    """
    Clas describing massifs (real from metadata or abstract from scratch)
    """

    def __init__(self, massifnum):
        self.massifnum = massifnum


class MassifReal(Massif):
    """
    Class describing an existing massif and its shape, taken from the shp snowtools_git/DATA/
    """

    def __init__(self, *args):
        Massif.__init__(self, *args)
        self.inProj = Proj(init = 'epsg:27572')  # EPSG code of the native NTF_Lambert_II_Carto
        self.outProj = Proj(init = 'epsg:4326')
        if self.massifnum <= 23:
            sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_alpes.shp')
        elif self.massifnum in [40, 41]:
            sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_corse.shp')
        elif self.massifnum > 91:
            # special case for locations outside massifs borders.
            raise Exception(' No true massifs with num {0} exists'.format(self.massifnum))
        else:
            sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_pyrenees.shp')
        for shape in sf.shapeRecords():
            if shape.record[1] == self.massifnum:
                self.name = shape.record[2]
                x = [i[0] for i in shape.shape.points[:]]
                y = [i[1] for i in shape.shape.points[:]]
                self.x, self.y = transform(self.inProj, self.outProj, x, y)
                self.center = [np.mean(self.x), np.mean(self.y)]
                break


class MassifAbs(Massif):
    """
    Cass describing abstract massifs such as single points outside massif borders
    """

    def __init__(self, massifnum, massifname, x, y, center):
        Massif.__init__(self, massifnum)
        self.name = massifname
        self.x = x
        self.y = y
        self.center = center
