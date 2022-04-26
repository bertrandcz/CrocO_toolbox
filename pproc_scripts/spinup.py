#!/usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 8 Feb 2020
copied from evalSODA2_CRST/pproc_scripts/spinup/py

@author: cluzetb
creating spinup files for semi-distributed assimilation
since we only have one year of forcing we loop ten times over it
apply a swe threshold across the summer.
'''

import datetime
import os
import shutil
from tools.change_prep import prep_tomodify
from utils.dates import check_and_convert_date

import tasks.s2m_command as s2m
from utilcrocO import area


# Params
beginDate = "2016080106"
endDate = "2017080106"
spinDate = endDate
sweceil = '100'
extractforcing = ""
vconf = '12'
forcingDir = os.environ['CROCOPATH'] + '/safran/' + vconf + '/ref/mb0000/meteo/'
forcingPath = forcingDir + 'FORCING_' + beginDate + '_' + endDate + '.nc'
rootXpOut = os.environ['CROCOPATH'] + '/s2m/' + vconf + '/spinup/'
gg = datetime.datetime.today().strftime("%Y%m%d-%H-%M-%S")
exesurfex = ' -s ' + os.environ['EXESURFEX']
diroutput = rootXpOut + gg
os.mkdir(diroutput)

nameName = 'OPTIONSV81_spinup.nam'
pathName = rootXpOut + nameName

# groud init : first year with -g arg
shortcommand = "s2m -b " + beginDate + " -e " + endDate + extractforcing + " -f " + forcingPath
commonoptions = " -o " + diroutput  + " -n " + pathName + ' -a ' + sweceil

commandTG = ' -g'
commandInit = shortcommand + commonoptions + commandTG
print(commandInit)
s2m(commandInit.split())

# the rest of the years without -g arg
for year in range(9):
    spincommand = ' -x ' + spinDate + ' -g'
    command = shortcommand + commonoptions + spincommand
    print(command)
    s2m(command.split())


# in the end, change the date and apply threshold on swe (dupli from snowtools_git/tests/test_thres_threshold.)

pathPrep = diroutput + '/prep/PREP_' + endDate + '.nc'
pathPgd = diroutput + '/prep/PGD.nc'
if not os.path.exists(rootXpOut + '/prep/'):
    os.makedirs(rootXpOut + '/prep/')
    os.makedirs(rootXpOut + '/pgd/')
pathNewPrep = rootXpOut + '/prep/' + 'PREP_' + beginDate + '.nc'
pathNewPgd = rootXpOut + '/pgd/' + 'PGD_' + area(vconf) + '.nc'
shutil.copyfile(pathPrep, pathNewPrep)
shutil.copyfile(pathPgd, pathNewPgd)

prepFile = prep_tomodify(pathNewPrep)

pp = prepFile
# apply a threshold on the SWE to prevent accumulation of snow across the summer at high altitude
pp.apply_swe_threshold(100)
pp.change_date(check_and_convert_date(beginDate))
pp.close()
