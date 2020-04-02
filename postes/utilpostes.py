#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 8 févr. 2020

@author: cluzetb

on beaufix, generate  conf file with the appropriate assimilation dates (for ex. every 7 days)

'''
import datetime

import netCDF4

import numpy as np


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
