#!/usr/bin/env/python3
# -*- coding: utf-8 -*-
'''
Created on 8 Feb 2020
copied from evalSODA2_CRST/scripts/spinup/py

@author: cluzetb
creating spinup files for semi-distributed assimilation
since we only have one year of forcing we loop ten times over it
'''

import datetime
import os
import shutil
import sys

import tasks.s2m_command as s2m
from tools.change_prep import prep_tomodify
from utils.dates import check_and_convert_date
# print(sys.path)
# Params
beginDate = "2016080106"
endDate = "2017080106"
massifs = [8, 9, 12, 13, 15, 16]
spinDate = endDate
# massifnumber = "postes"  # 12: gdes rousses, 13 : thabor
sweceil = '100'
extractforcing = ""
vconf = 'postes_' + '_'.join(list(map(str, massifs))) + '_csv'
forcingDir = os.environ['VORTEXPATH'] + '/safran/' + vconf + '/ref/mb0000/meteo/'
forcingPath = forcingDir + 'FORCING_' + beginDate + '_' + endDate + '.nc'
rootXpOut = '/home/cluzetb/vortexpath/s2m/' + vconf + '/spinup/'
gg = datetime.datetime.today().strftime("%Y%m%d-%H-%M-%S")
exesurfex = ' -s /home/cluzetb/SURFEX_V81/cen_release/exe/'
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
s2m.execute(commandInit.split())
# os.remove(diroutput + '/prep/init_TG.nc')

# the rest of the years without -g arg
for year in range(9):
    spincommand = ' -x ' + spinDate + ' -g'
    command = shortcommand + commonoptions + spincommand
    print(command)
    s2m.execute(command.split())


# in the end, change the date and apply threshold on swe (dupli from snowtools_git/tests/test_thres_threshold.)

# pathPrep = '/home/cluzetb/vortexpath/s2m/12/spinup_V81/20181004113826822174/prep/PREP_' + endDate + '.nc'
pathPrep = diroutput + '/prep/PREP_' + endDate + '.nc'
pathPgd = diroutput + '/prep/PGD.nc'
if not os.path.exists(rootXpOut + '/prep/'):
    os.makedirs(rootXpOut + '/prep/')
    os.makedirs(rootXpOut + '/pgd/')
pathNewPrep = rootXpOut + '/prep/' + 'PREP_' + beginDate + '.nc'
pathNewPgd = rootXpOut + '/pgd/' + 'PGD_postes.nc'
shutil.copyfile(pathPrep, pathNewPrep)
shutil.copyfile(pathPgd, pathNewPgd)

prepFile = prep_tomodify(pathNewPrep)

pp = prepFile

pp.apply_swe_threshold(100)
pp.change_date(check_and_convert_date(beginDate))
pp.close()
