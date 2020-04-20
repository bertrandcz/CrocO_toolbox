#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 28 févr. 2019

@author: cluzetb
This script launches series of crampon command to test the PF collapsing conditions.


'''
import glob
import os

import matplotlib.pyplot as plt
xpid = 'artB31D11_2016@cluzetb'
dt = '2017031310'
rootdir = '/home/cluzetb/vortexpath/s2m/12/' + xpid + '/crampon'
os.chdir(rootdir)
plt.figure()
for filename in glob.glob('./coll*/' + dt + '/PART'):
    f = open(filename)
    config = filename[6:-16]

    # print filename
    print(config)
    fact = config.split('_')[0]
    npts = config.split('_')[1]
    for line in f:
        part = line.split()
        unpart = len(set(part))
    f.close()
    # print fact
    # print npts
    print(unpart)

    if float(unpart) > 27.:
        color = 'k'
    elif float(unpart) < 3.:
        color = 'r'
    else:
        color = 'b'
    plt.scatter(float(fact), float(npts), s=unpart * 4, color = color)
    plt.title('cardinal of the PF sample, ' + dt)
    plt.grid(True)
    # Create some sizes and some labels.
    labels = ['Tiny', 'Small', 'Medium', 'Large', 'Huge']
    plt.xlabel('multiplicative factor on observation error')
    plt.xlim([0, 26])
    plt.ylabel('number of classes simultaneously assimilated')

plt.savefig('/home/cluzetb/snowtools_git/assim/scripts/cardinal_' + dt + '.png')
