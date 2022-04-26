# /usr/bin/env/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:07:26 2020

@author: cluzetb
"""

import os

import netCDF4  # @UnresolvedImport

import numpy as np
import xml.etree.ElementTree as ET


def find_metadata(filepath = os.environ['SNOWTOOLS_CEN'] + '/DATA/METADATA.xml', massifs = [12], exclude_corsica = True):
    '''
    BC 17/01/20
    Extract all metadata of postes (Sites) within the list of massifs.
    '''
    tree = ET.parse(filepath)
    root = tree.getroot()

    listp = []
    for site in root[1].findall('Site'):
        for massif in massifs:
            if str(massif) == site.find('massif').text.strip():
                listp.append(site)
    listpnb = [int(p.find('number').text) for p in listp]
    listpname = [p.find('name').text for p in listp]

    # reorder lists so that they are sorted by station number
    # build a dict:
    pdict = {nb: [p, name] for nb, p, name in zip(listpnb, listp, listpname)}
    listpnb = [k for k in sorted(pdict.keys())]
    listp = [pdict[k][0] for k in sorted(pdict.keys())]
    listpname = [pdict[k][1] for k in sorted(pdict.keys())]
    print('###')
    print('number of selected posts: ', len(listp))
    return listp, listpnb, listpname


def find_name_station(stations, filepath=os.environ['SNOWTOOLS_CEN'] + '/DATA/METADATA.xml'):
    '''
    BC 17/02/20
    find the string name of a station ID or list of IDS.
    '''

    tree = ET.parse(filepath)
    root = tree.getroot()
    ret = []
    # usually, stations is a numpy array or an int32...
    if isinstance(stations, list):
        pass
    elif np.issubdtype(stations, np.integer):  # @TODO: still a bug here for int...
        stations = [stations]
    for station in stations:
        for site in root[1].findall('Site'):
            if str(station) in site.find('number').text:
                name = site.find('name').text.strip()
        try:
            ret.append(name)
            del name
        except UnboundLocalError:
            print('there is no station ', station, 'in the database ', filepath)

    if len(ret) == 1:
        return ret[0]
    else:
        return ret


def find_duplicates(listp, listpnb, listpname, radius = 1000):
    """
    BC July 2020
    finding the lat/lon duplicates within a 1km radius in a list of posts
    """

    listplatlon = [(float(p.find('lat').text.strip()), float(p.find('lon').text))  for p in listp]
    names = dict()
    numbers = dict()
    dupes = set()
    for _, ll in enumerate(listplatlon):
        dists = distances(ll, listplatlon)
        in_radius = np.where(dists < radius)[0]
        # print(' in_radius', in_radius)
        if len(in_radius) > 1:
            in_radius = list(map(int, in_radius))
            dupes.add(listplatlon[in_radius[0]])
            names[in_radius[0]] = [listpname[ir] for ir in in_radius]
            numbers[in_radius[0]] = [listpnb[ir] for ir in in_radius]
    return names, numbers, dupes


def distances(latlon, listplatlon):
    """distance in m between latlon and each point of listplatlon"""
    def dist_deg(ll1, ll2):
        return np.sqrt((ll1[0] - ll2[0])**2 + (ll1[1] - ll2[1])**2)
    R_t = 6371000  # earth's radius
    dists = np.array([dist_deg(latlon, llp) for llp in listplatlon]) * R_t / 180. * np.pi
    return dists


def slice_file_listpnb(fileIn, fileOut, listpnb):
    '''
    Extract a list of points from a forcing or PRO (provided they have a station variable)
    '''
    forcing = netCDF4.Dataset(fileIn, 'r')
    station = forcing.variables['station'][:]
    statpos = {s: str(i) for i, s in enumerate(station) if s in listpnb}
    forcing.close()
    listfiles = []
    # /!\ CRUCIAL POINT HERE: sorting of the station by number.
    for s in sorted(statpos.keys()):
        ff = '/tmp/out_' + str(s) + '.nc'
        listfiles.append(ff)
        cmd = 'ncks -O -d Number_of_points,' + statpos[s] + ',' + statpos[s] + \
            ' ' + fileIn + ' ' + ff
        os. system(cmd)

    # make number of points record dim
    [os.system('ncpdq -O -a Number_of_points,time ' + f + ' ' + f) for f in listfiles]
    cmd2 = 'ncrcat -O ' + ' '.join(listfiles) + ' ' + fileOut
    print(cmd2)
    os.system(cmd2)
    cmd11 = 'ncpdq -O -a time,Number_of_points ' + fileOut + ' ' + fileOut
    os.system(cmd11)
    [os.remove(f) for f in listfiles]
    return 0
