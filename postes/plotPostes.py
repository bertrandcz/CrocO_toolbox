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
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import matplotlib as mpl
from cartopy.io import srtm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io import LocatedImage
from matplotlib.ticker import FormatStrFormatter
import math


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
        # EPSG code of the native NTF_Lambert_II_Carto
        self.inProj = Proj(init='epsg:27572')
        self.outProj = Proj(init='epsg:4326')
        mass = mountain_from_massif([self.massifnum])[0]
        mmass = mass if 'orlu' not in mass else 'pyrenees'
        sf = shp.Reader(os.environ['SNOWTOOLS_CEN'] + '/DATA/massifs_' +
                        mmass + '.shp')

        for shape in sf.shapeRecords():
            if shape.record[1] == self.massifnum:
                self.name = shape.record[2]
                x = [i[0] for i in shape.shape.points[:]]
                y = [i[1] for i in shape.shape.points[:]]
                self.x, self.y = transform(self.inProj, self.outProj, x, y)
                self.polygon = Polygon([(xx, yy)
                                        for xx, yy in zip(self.x, self.y)])
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
        self.mountains = {0: 'alpes', 1: 'pyrenees'}
        self.axis = {}
        self.xlims = {}
        self.ylims = {}
        self.extent = {}
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

    def prepare_background(self,
                           ax=None,
                           background=True,
                           tag_massifs=False,
                           resolution='high'):
        """
        Bc june 2021 quick hack to submit my paper.
        """
        if ax is not None:
            raise Exception(' not implemented yet, please be patient Bber')
        projection = ccrs.PlateCarree()
        self.fig = plt.figure(figsize=cm2inch(12, 10), constrained_layout=True)
        # height aspect ratio btw alps/pyr is 2.9888... :)
        self.gs = self.fig.add_gridspec(2, 10, height_ratios=[3, 1])
        gs = self.gs
        plt.subplots_adjust(top=.92, bottom=.03)
        for i, mountain in self.mountains.items():
            if i == 0:
                self.axis[mountain] = self.fig.add_subplot(
                    gs[i, 0:7], projection=projection)
            else:
                self.axis[mountain] = self.fig.add_subplot(
                    gs[i, :], projection=projection)
            ax = self.axis[mountain]
            ax.set_anchor('W')
            gl = ax.gridlines(zorder=3, alpha=0.5)
            gl.xlabels_top = False
            gl.xlabels_top = True
            gl.ylabels_left = True
            gl.ylabels_right = False
            gl.xlines = True
            if i == 0:
                gl.xlocator = mticker.FixedLocator(np.arange(5, 8.5, 1.))
                gl.ylocator = mticker.FixedLocator(np.arange(43.5, 47., 0.5))
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
            else:
                gl.xlocator = mticker.FixedLocator(np.arange(-1.5, 3.5, 1.))
                gl.ylocator = mticker.FixedLocator(np.arange(42, 44, 0.5))
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER

            for mnum, massif in self.massifs.items():
                if mountain_from_massif([mnum]) == [mountain]:
                    ax.plot(massif.x, massif.y, color='k', zorder=2)
                # trick to limit to the alps/pyrenees

                if tag_massifs:
                    ax.annotate(massif.name, massif.center)

            self.xlims[mountain] = ax.get_xlim()
            self.ylims[mountain] = ax.get_ylim()

            self.extent[mountain] = self.xlims[mountain] + self.ylims[mountain]
            ax.set_extent(self.extent[mountain])
            if background:
                # Add a background image downloaded installed in .local/.../cartopy and converted into .png (see cartopy doc)
                ax.background_img(extent=self.extent[mountain],
                                  name='ne_shaded_high',
                                  resolution=resolution)

    def plot_scatter_postes(self,
                            tag_postes=False,
                            tag_massifs=True,
                            plotCircle=True,
                            cmap='RdBu',
                            vmin=-0.5,
                            vmax=.5,
                            colorbar=True,
                            cbar_steps=4,
                            cbar_title=''):
        for i, mountain in self.mountains.items():
            ax = self.axis[mountain]
            _ = ax.scatter(self.sdObj.pgd.lon,
                           self.sdObj.pgd.lat,
                           c=self.sdObj.data['DEP'],
                           cmap=cmap,
                           vmin=vmin,
                           vmax=vmax,
                           zorder=100,
                           s=6)
            if tag_postes:
                for ipt in range(self.sdObj.pgd.npts):
                    ax.annotate(
                        str(self.sdObj.pgd.station[ipt]),
                        (self.sdObj.pgd.lon[ipt], self.sdObj.pgd.lat[ipt]),
                    )

            # plot a circle on the alps
            if i == 0:
                if 'local' in self.sdObj.ptinom and plotCircle:

                    radius = callunits(
                        self.sdObj.ptinom.split(',')[1].split(' ')[0])
                    circle = plt.Circle(
                        (self.xlims[mountain][0] + 1.1 * radius,
                         self.ylims[mountain][1] - 1.1 * radius),
                        radius,
                        color='g',
                        fill=False,
                        lw=3)
                    ax.add_artist(circle)
                    angle_rad = 0
                    ax.arrow(self.xlims[mountain][0] + 1.1 * radius,
                             self.ylims[mountain][1] - 1.1 * radius,
                             (radius - 0.1) * math.cos(angle_rad),
                             (radius - 0.1) * math.sin(angle_rad),
                             head_width=0.1,
                             head_length=0.1,
                             fc='g',
                             ec='g')
                    coords = (self.xlims[mountain][0] + 1.2 * radius,
                              self.ylims[mountain][1] - 0.75 * radius)
                    print(coords)
                    ax.annotate('35 km',
                                coords,
                                xycoords='data',
                                ha='center',
                                va='center',
                                fontsize=7, c = 'g')

        if colorbar:
            self.add_colorbar(vmin=vmin,
                              vmax=vmax,
                              cmap=cmap,
                              cbar_steps=cbar_steps,
                              cbar_title=cbar_title)

    def plot_massifs_color(self,
                           dictmassif,
                           ax=None,
                           vmin=0,
                           vmax=1,
                           cmap='viridis',
                           colorbar=False,
                           cbar_steps=4,
                           cbar_title=''):
        """
        BC, 03/09/20
        plot the polygons of the massifs filled with a value from dictmassif

        """
        cmap = plt.cm.__dict__[cmap]
        norm = plt.Normalize(vmin, vmax)
        for _, mountain in self.mountains.items():
            ax = self.axis[mountain]
            patches = []
            for mnum in self.massifs:
                if mountain_from_massif([mnum])[0] == mountain:
                    color = cmap(norm(dictmassif[mnum]))
                    self.massifs[mnum].polygon.set_color(color)
                    patches.append(self.massifs[mnum].polygon)
            pc = PatchCollection(patches,
                                 match_original=True,
                                 edgecolor='k',
                                 linewidths=1.,
                                 zorder=2)
            ax.add_collection(pc)
            ax.autoscale()
        if colorbar:
            self.add_colorbar(vmin=0,
                              vmax=np.max(list(dictmassif.values())),
                              cmap=cmap,
                              cbar_steps=cbar_steps,
                              cbar_title=cbar_title)

    def add_colorbar(self,
                     vmin=-1,
                     vmax=1,
                     cmap='RdBu',
                     cbar_steps=4,
                     cbar_title=''):
        step = (vmax - vmin) / cbar_steps
        #ax = self.fig.add_subplot(self.gs[0, 8:9])
        gssub = self.gs[0, 8:9].subgridspec(2, 1, height_ratios=[1, 8])
        ax = plt.subplot(gssub[0])
        _ = ax.text(0.7,
                    0.,
                    cbar_title,
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)
        ax.axis('off')
        ax = plt.subplot(gssub[1])
        ax.set_anchor('W')
        norm = plt.Normalize(vmin, vmax)
        gggg = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
        gggg.set_array(np.arange(vmin, vmax, 5))  # divide by 10 for obs/km/yr
        cb = plt.colorbar(gggg,
                          cax=ax,
                          ticks=np.arange(vmin, vmax + step, step))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


def shade(located_elevations):
    """
    Given an array of elevations in a LocatedImage, add a relief (shadows) to
    give a realistic 3d appearance.

    """
    new_img = srtm.add_shading(located_elevations.image,
                               azimuth=135,
                               altitude=15)
    return LocatedImage(new_img, located_elevations.extent)
