# /usr/bin/env/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 15:07:26 2020

@author: cluzetb
"""

import os
import time

import netCDF4

import xml.etree.ElementTree as ET


def find_metadata(filepath = '/home/cluzetb/snowtools_git/DATA/METADATA.xml', massifs = [12]):
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


def find_name_station(station, filepath = '/home/cluzetb/snowtools_git/DATA/METADATA.xml'):
    '''
    BC 17/02/20
    find the string name of a station ID
    '''
    tree = ET.parse(filepath)
    root = tree.getroot()
    listp = []
    for site in root[1].findall('Site'):
        if str(station) in site.find('number').text:
            ret = site.find('name').text.strip()
    try:
        return ret
    except UnboundLocalError:
        print('there is no station ', station, 'in the database ', filepath)


def slice_file_listpnb(fileIn, fileOut, listpnb):
    '''
    Extract a list of points from a forcing or PRO (provided they have a station variable)
    '''
    forcing = netCDF4.Dataset(fileIn, 'r')
    station = forcing.variables['station']
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
