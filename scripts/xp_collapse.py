#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 28 févr. 2019

@author: cluzetb
This script launches series of crocO command to test the PF collapsing conditions.


'''
import os
import time
start_time = time.time()


# ####### FIXED EXPERIMENT SETTINGS ###############
xpid = 'artB31D11_2016@cluzetb'
dt = '2017010611'
var = 'b4,b5'
# todo = 'None' BC no todo to do.
nmembers = '35'
synth = '0'
list_fact = [1., 1.1, 1.5, 2., 5., 10., 25., 100. ]
# list_fact = [1.]
# #################################################
for fact in list_fact:
    print(("#########factor ", fact))
    for npts in range(1, 45):
        print(('npts ', npts))
        o = 'coll' + str(fact) + '_' + str(npts)

        pts = list(range(140, 140 + npts))
        ipts = ','.join([str(pp) for pp in pts])
        os.system('python /home/cluzetb/snowtools_git/assim/crocO.py ' +
                  '--xpid ' + xpid +
                  ' -d ' + dt +
                  ' --vars ' + var +
                  ' --synth ' + synth +
                  ' -o ' + o +
                  ' --nmembers ' + nmembers +
                  ' --fact ' + str(fact) +
                  ' --classesId ' + ipts
                  )
        elapsed_time = time.time() - start_time
        print(elapsed_time)
