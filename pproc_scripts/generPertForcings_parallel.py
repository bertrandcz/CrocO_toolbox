# -*- coding: utf-8 -*-
'''
Created on 10 May 2020
@author: cluzetb.

Launch LOCALLY the perturbation of the forcings using snowtools facilities
(see snowtools.tools.job_gener_pert_forcings.py for a launch on beaufix)
If you have access to beaufix, you'd better use snowtools_git/tools/job_pert_forcings.sh
(parallelized + direct access to the archive via transfer jobs)
'''
from optparse import OptionParser
import os
from shutil import copyfile
import sys
from tools.job_gener_pert_forcings import multiprocess
from utilcrocO import split_list


def parse_options(arguments):
    parser = OptionParser(description = "Generate ensemble weather forcing.")
    parser.add_option("--vconf", dest = "vconf", type=str,
                      help = 'full name of the geometry')
    parser.add_option("-y", dest = "years", type=str,
                      help = 'comma separated years', action = 'callback', callback = split_list)
    parser.add_option('--nens', dest = 'nens', type = int,
                      help = 'ensemble size')
    parser.add_option('--param', dest = 'param', type = str,
                      help = ' path to the param file', default = os.environ['SNOWTOOLS_CEN'] + '/tools/param.txt')
    parser.add_option('--no_brutal_imp', dest = 'no_brutal_imp',
                      action = 'store_true', default = False, help = ' deactivate brutal factors on impurities')
    options, _ = parser.parse_args(arguments)
    return options


if __name__ == '__main__':
    parser = parse_options(sys.argv)

    machine = os.uname()[1]
    print(machine)
    if machine == 'sxcen':
        refPath = '/home/cluzetb/Data/safran/' + parser.vconf + '/ref/mb0000/meteo/'
        rootPath = '/cnrm/cen/users/NO_SAVE/cluzetb/vortex/safran/' + parser.vconf + '/'
    else:
        refPath = os.environ['CROCOPATH'] + '/safran/' + parser.vconf + '/ref/mb0000/meteo/'
        rootPath = os.environ['CROCOPATH'] + '/safran/' + parser.vconf + '/'

    for year in parser.years:
        forcName = refPath + 'FORCING_{0}080106_{1}080106.nc'.format(year, year + 1)
        outDir = rootPath + 'forcing_{0}{1}B_31D_11_t1500_160/'.format(year, year + 1)
        if not os.path.exists(outDir):
            os.mkdir(outDir)
        copyfile(parser.param, outDir + 'param.txt')
        # run the perturbation in parallel
        brutalImp = not parser.no_brutal_imp
        multi = multiprocess(forcName, parser.param, parser.nens, outDir, brutalImp)
