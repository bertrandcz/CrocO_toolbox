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
import shapefile as shp


class MapPostes(Pie):
    '''
    Inspired on Pie, the principle is to plot maps of a variable (e.g. correlation or score)
    in the postes geometry: the postes denoted by points are colored according to the variable value.
    The plotting map facility comes from notebooks.postes.first_glance*.ipynb
    '''

    def __init__(self, *args, **kwargs):
        MapPostes.__init__(self, *args, **kwargs)
        # load shapefile of massifs from the alps and pyrenees into a same shp.
        # (https://gis.stackexchange.com/questions/103033/using-pyshp-to-merge-two-shapefiles)

        self.setMassifShapes()

    def setMassifsShapes(self):
        self.massifs = dict()
        for massif in np.unique(self.sdObj.pgd.massif):
            self.massifs[massif] = Massif(massif)

    def plot(self):
        fig, ax = plt.subplots()

        for mnum, massif in self.massifs.items():
            ax.plot(massif.x, massif.y, color = 'k', zorder = 0)
            ax.annotate(massif.name, massif.center)


class Massif(object):
    """
    Class describing a massif and its shape, taken from the shp snowtools_git/DATA/
    """

    def __init__(self, massif):
        self.inProj = Proj(init = 'epsg:27572')  # EPSG code of the native NTF_Lambert_II_Carto
        self.outProj = Proj(init = 'epsg:4326')

        if massif <= 23:
            sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_alpes.shp')
        elif massif in [40, 41]:
            sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_corse.shp')

        else:
            sf = shp.Reader(os.environ['SNOWTOOLS_CEN' ] + '/DATA/massifs_pyrenees.shp')
        for shape in sf.shapeRecords():
            if shape.record[1] == massif:
                self.name = shape.record[2]
                x = [i[0] for i in shape.shape.points[:]]
                y = [i[1] for i in shape.shape.points[:]]
                self.x, self.y = transform(self.inProj, self.outProj, x, y)
                self.center = [np.mean(x), np.mean(y)]
                break
