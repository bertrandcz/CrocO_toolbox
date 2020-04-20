# -*- coding: utf-8 -*-
'''
Created on 5 feb. 2019

@author: cluzetb

perform local tests/dev of SODA based on :
    - set of prep files from an openloop
    - opt : set of REAL observations OR generate synthetical observations.

'''
from optparse import OptionParser
import os
import sys
import time

from CramponPf import CramponPf
from postes.utilpostes import set_conf_everydate
from utilcrampon import read_conf


usage = 'crocO --opts'


def parse_options(arguments):

    parser = OptionParser(usage)

    parser.add_option("--xpid",
                      type="string", action="store", dest="xpid", default=None,
                      help="xpid of the crampon run")
    parser.add_option("--xpidol",
                      type="string", action="store", dest="xpidol", default=None,
                      help="xpid of the associated openloop crampon run (if any)")
    parser.add_option("--xpidobs",
                      type="string", action="store", dest="xpidobs", default=None,
                      help="xpid of the associated observations")
    parser.add_option("-d",
                      type = 'string', action = 'callback', callback = callvars, dest='dates', default = None,
                      help = "specify assimilation dates to try (yyyymmddhh) or all for all conf.assimdates.")
    # parser.add_option("--evid",
    #                   action="store", type="string", dest="evid", default="",
    #                   help="xpev (useful for pickle save")
    parser.add_option("--kind",
                      action="store", type="string", dest="kind", default='beaufixpp',
                      help='kind of the run run to perform (beaufixpp for CramponPp')
    parser.add_option("--vapp",
                      action="store", type="string", dest="vapp", default='s2m',
                      help="vapp of the soda-vortex task")
    parser.add_option("--vconf",
                      action="store", type="string", dest="vconf", default='12',
                      help="vconf of the soda-vortex task")
    # parser.add_option("--user",
    #                   action='store', type='string', dest='user', default='cluzetb',
    #                   help= 'user where to find the files')
    parser.add_option("--todo",
                      action = 'store', type = 'string', dest='todo', default = None,
                      help= 'type of evaluation')
    parser.add_option("--nmembers",
                      action = 'store', type = int, dest='nmembers', default = None,
                      help= 'nb members for post proc.')
    parser.add_option("--distr",
                      action = 'store_true', dest = 'distr', default = False,
                      help = 'specify if your simulation is distributed or not')
    parser.add_option("--vars", type = 'string', action="callback", callback=callvars, default = 'all',
                      help="specify assimilation variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snowdepth")
    parser.add_option("--ppvars", type = 'string', action="callback", callback=callvars, default = [],
                      help = 'specify study variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snowdepth")')
    parser.add_option("--fact", type = 'string', action = 'callback', callback = callvars, default = 1.,
                      help = ' set multiplicative factors for the observation errors')
    parser.add_option("--classesS", type = 'string', action="callback", callback=callvars, default = 'all',
                      help="specify analyzed slope classes ex : 0,20,40")
    parser.add_option("--classesA", type = 'string', action="callback", callback=callvars, default = 'all',
                      help="specify analyzed aspect classes ex : N,NE,E,SE,S,SO")
    parser.add_option("--classesE", type = 'string', action="callback", callback=callvars, default = 'all',
                      help="specify analyzed elev classes ex : 1800,2100,2400,2700")
    parser.add_option("--classesId", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed ids of classes separated by commas ex : 145,146,147")
    parser.add_option("-o",
                      action = 'store', type = 'string', dest='saverep', default = "",
                      help= 'name of the postproc/type_of_anal/ without /')
#     parser.add_option("--mpi",
#                       action = 'store_true', dest = 'mpi', default = False,
#                       help='activate MPI for soda')
    parser.add_option("--pf",
                      action = 'store', dest = 'pf', default = 'global', type = 'string',
                      help='choose pf agorithm : global, rlocal, klocal')
    parser.add_option("--nloc_pf", action = 'store', dest = 'nloc_pf', default = 0, type = int,
                      help='set the localization radius (1->+oo) or the k-localization counter (1->+oo)')
    parser.add_option("--neff", action = 'store', dest = 'neff', default = 0, type = int,
                      help='set the Neff target value for the inflation factor')
    parser.add_option('--synth',
                      action = 'store', dest = 'synth', type = int, default = None,
                      help ='choose member as a synthetic observation. TODO : 0 for random.')
    parser.add_option('--noise',
                      action = 'store', dest = 'noise', type = float, default = 0.,
                      help ='multiplicative gaussian noise to the synthetical obs.')
    parser.add_option('--sensor',
                      action = 'store', dest = 'sensor', type = 'string', default = None,
                      help ='provide sensor name for OBS reading (MODIS, pleiades...).')
    parser.add_option("--archive_synth",
                      action = 'store_true', dest = 'archive_synth', default = False,
                      help = 'specify if you want to archive the synthetical obs with its mask or not.')
    parser.add_option("--need_masking",
                      action = 'store_false', dest = 'need_masking', default = True,
                      help = 'specify if you need to mask the obs or use it as it is;')
    parser.add_option("--readprep",
                      action = 'store_true', dest = 'readprep', default = False,
                      help = 'read prep and obs')
    parser.add_option("--readaux",
                      action = 'store_true', dest = 'readaux', default = False,
                      help = 'read auxiliary files')
    parser.add_option("--readobs",
                      action = 'store_true', dest = 'readobs', default = False,
                      help = 'read observation files')
    parser.add_option("--readoper",
                      action = 'store_true', dest = 'readoper', default = False,
                      help = 'read ooper pickle file')
    parser.add_option("--notreadpro",
                      action = 'store_true', dest = 'notreadpro', default = False,
                      help = 'read observation files')
    parser.add_option("--clim",
                      action = 'store_true', dest = 'clim', default = False,
                      help = 'read the clim')
    (options, args) = parser.parse_args(arguments)

    del args

    return options


def callvars(option, opt, value, parser):
    if type(value) is not float:  # normal case
        setattr(parser.values, option.dest, value.split(','))
    else:  # bug with float args like --fact
        setattr(parser.values, option.dest, value)


def set_options(args, pathConf = None):
    options = parse_options(args)

    if 'VORTEXPATH' not in list(os.environ.keys()):
        raise Exception('you must export VORTEXPATH to the root of your local vortex xps.')
    else:
        options.vortexpath = os.environ['VORTEXPATH']
    if options.xpid is not None:
        options.xpiddir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpid + '/'
    if options.xpidobs is not None:
        options.xpidobsdir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/obs/' + options.xpidobs + '/'
    if options.xpidol is not None:
        options.xpidoldir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpidol + '/'
    if pathConf is None:
        try:
            confPath = options.xpiddir + '/conf/' + options.vapp + '_' + options.vconf + '.ini'
            conf = read_conf(confPath)
        except ValueError as e:
            # if no conf file can be found, fake it
            print('faking a new conf file.')
            if not os.path.exists(options.xpiddir + '/conf'):
                os.makedirs(options.xpiddir + '/conf', exist_ok = True)
            ok = set_conf_everydate(2013, 1, confPath, nens = 40, endY = 2018)
            conf = read_conf(confPath)
    else:
        conf = read_conf(pathConf)

    if 'all' in options.dates:
        # in the case we are pping, last date corresponds to the end of the simulation and is not interesting.
        if hasattr(conf, 'stopdates'):
            options.dates = conf.stopdates[0:-1]
        # otherwise:
        else:
            options.dates = conf.assimdates[0:-1]
    elif type(options.dates) is str:
        options.dates = [options.dates]
    return options, conf


def execute(args):
    options, conf = set_options(args)
    run = CramponPf(options, conf)
    if options.todo in ['corr', 'load', 'generobs']:
        pass
    elif options.todo in ['pf', 'pfpython']:
        run.run()

        run.post_proc(options)
    elif options.todo is not None:
        pb = run.post_proc(options)
        return pb
    else:
        run.run()


if __name__ == '__main__':
    start_time = time.time()

    execute(sys.argv)
    elapsed_time = time.time() - start_time
    print(elapsed_time)
