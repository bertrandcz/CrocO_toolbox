'''
Created on 30 mai 2020

@author: cluzetb
Advances plotting facilities suited to the postes geometry
'''
import os

from pyproj import Proj, transform

import numpy as np
import shapefile as shp  # @UnresolvedImport
from postes.utilpostes import mountain_from_massif
from matplotlib.patches import Polygon
from plotcrocO import Pie
import matplotlib.pyplot as plt
from crocO import callunits
from matplotlib.collections import PatchCollection
from utilplot import cm2inch
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap


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
        mass = mountain_from_massif([self.massifnum])[0]
        mmass = mass if 'orlu' not in mass else 'pyrenees'
        sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_' + mmass + '.shp')

        for shape in sf.shapeRecords():
            if shape.record[1] == self.massifnum:
                self.name = shape.record[2]
                x = [i[0] for i in shape.shape.points[:]]
                y = [i[1] for i in shape.shape.points[:]]
                self.x, self.y = transform(self.inProj, self.outProj, x, y)
                self.polygon = Polygon([(xx, yy) for xx, yy in zip(self.x, self.y)])
                self.center = [np.mean(self.x), np.mean(self.y)]
                self.area = shape.record[4]
                break


class MassifAbs(Massif):
    """
    Class describing abstract massifs such as single points outside massif borders
    """

    def __init__(self, massifnum, massifname, x, y, center):
        Massif.__init__(self, massifnum)
        self.name = massifname
        self.x = x
        self.y = y
        self.center = center


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

    def plot(self, ax = None, tag_postes = False, tag_massifs = True,
             plotCircle = True, vmin = -1, vmax = 1, cmap = 'RdBu', colorbar = False, background = False):
        if ax is None:
            fig = plt.figure(figsize = cm2inch(12.5, 6))
            gs = gridspec.GridSpec(1, 2, width_ratios=[0.97, 0.03],
                                   wspace = 0.05,
                                   bottom = 0.08,
                                   top = 0.95,
                                   left = 0.08,
                                   right = 0.9,
                                   )

            ax = plt.subplot(gs[0])

        else:
            fig = plt.gcf()
        for mnum, massif in self.massifs.items():
            ax.plot(massif.x, massif.y, color = 'k', zorder = 2)
            if tag_massifs:
                ax.annotate(massif.name, massif.center)
        gg = ax.scatter(self.sdObj.pgd.lon, self.sdObj.pgd.lat, c = self.sdObj.data['DEP'],
                        cmap = cmap, vmin = vmin, vmax = vmax, zorder = 100, s=6)
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
        xmin = ax.get_xlim()[0]
        ymax = ax.get_ylim()[1]
        ymin = ax.get_ylim()[0]
        xmax = ax.get_xlim()[1]
        # plot a circle
        if 'local' in self.sdObj.ptinom and plotCircle:

            radius = callunits(self.sdObj.ptinom.split(',')[1].split(' ')[0])
            circle = plt.Circle((xmin + 1.1 * radius, ymax - 1.1 * radius), radius, color = 'g', fill = False)
            ax.add_artist(circle)
        ax.set_aspect('equal')

        # shaded background
        if background:

            m = Basemap(projection='cyl', llcrnrlat=ymin, llcrnrlon=xmin,
                        urcrnrlat=ymax, urcrnrlon=xmax,
                        resolution='i')
            m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000)
            m.drawrivers()
            m.drawcountries()
            m.drawparallels(np.arange(42., 47.5, 1.), labels=[True, False, False, True])
            m.drawmeridians(np.arange(-1., 9.5, 1.), labels=[False, True, True, False])

        # colorbar
        if colorbar:
            step = (vmax - vmin) / 10
            cb = plt.colorbar(gg, cax=plt.subplot(gs[1]), ticks=np.arange(vmin, vmax + step, step))
        # title
        # plt.title(self.sdObj.ptinom)

    def plot_massifs_color(self, dictmassif, ax = None, vmin = 0, vmax = 1, cmap = 'viridis'):
        """
        BC, 03/09/20
        plot the polygons of the massifs filled with a value frm dictmassif

        """
        # print(plt.cm.__dict__)
        cmap = plt.cm.__dict__[cmap]
        norm = plt.Normalize(vmin, vmax)

        patches = []
        for massif in self.massifs:
            color = cmap(norm(dictmassif[massif]))
            self.massifs[massif].polygon.set_color(color)
            patches.append(self.massifs[massif].polygon)
        pc = PatchCollection(patches,
                             match_original=True,
                             edgecolor='k',
                             linewidths=1.,
                             zorder=2
                             )
        # print(pc.__dict__)
        ax.add_collection(pc)
        ax.autoscale()
