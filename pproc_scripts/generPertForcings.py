# -*- coding: utf-8 -*-
'''
Created on 8 nov. 2018
@author: cluzetb. Base of the script comes from Cesar DB
Launch locally the perturbation of the forcings using snoztpols facilities
If you have access to beaufix, you'd better use snowtools_git/tools/job_pert_forcings.sh (parallelized)
'''
import os
from shutil import copyfile
import sys
machine = os.uname()[1]
print(machine)
if machine == 'sxcen':
    refPath = '/home/cluzetb/Data/safran/12/ref/mb0000/meteo/'
    rootPath = '/cnrm/cen/users/NO_SAVE/cluzetb/vortex/safran/12/'
else:
    refPath = os.environ['CRAMPONPATH'] + '/safran/12/ref/mb0000/meteo/'
    rootPath = os.environ['CRAMPONPATH'] + '/safran/12/'

for year in [2013, 2014, 2015, 2016]:
    forcName = refPath + 'FORCING_{0}080106_{1}080106.nc'.format(year, year + 1)
    print(sys.argv)
    paramPath = sys.argv[1] if len(sys.argv) > 1 else os.environ['SNOWTOOLS_CEN'] + '/tools/param.txt'
    nmembers = 160
    brutalImp = '--brutalImp'
    outDir = rootPath + 'forcing_{0}{1}B_31D_11_t1500_160/'.format(year, year + 1)
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    copyfile(paramPath, outDir + 'param.txt')
    # run the perturbation
    cmd = os.environ['SNOWTOOLS_CEN'] + '/tools/makeForcingEnsemble.py -f ' + forcName + ' -p ' + paramPath\
        + ' -nmembers ' + str(nmembers)  + ' -o ' + outDir + ' ' + brutalImp
    os.system(cmd)

    # visualize
    if not os.path.exists(outDir + "viz"):
        os.mkdir(outDir + 'viz')
    cmd2 = os.environ['SNOWTOOLS_CEN'] + 'tools/VizForEnsemble.py -r ' + outDir + ' -o ' + outDir + 'viz'
    os.system(cmd2)
