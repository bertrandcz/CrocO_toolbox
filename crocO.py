# -*- coding: utf-8 -*-
'''
Created on 5 feb. 2019

@author: cluzetb

perform local tests/dev of SODA based on :
    - set of prep files from an openloop
    - opt : set of REAL observations OR generate synthetical observations.

'''
from CrocOrun import CrocOrun
from optparse import OptionParser
import os
import random
import sys
import time
from utilcrocO import read_conf, ImmutableOpt, Pgd, area, parse_classes


usage = 'crocO --opts'


def parse_options(arguments):

    parser = OptionParser(usage)

    parser.add_option("--xpid",
                      type="string", action="store", dest="xpid", default=None,
                      help="xpid of the crampon run : must exist (if you are not post-processing, you must create it before launching the experiment")
    parser.add_option("--xpidol",
                      type="string", action="store", dest="xpidol", default=None,
                      help="xpid of the associated openloop crampon run (if any)")
    parser.add_option("--xpidobs",
                      type="string", action="store", dest="xpidobs", default=None,
                      help="xpid of the associated observations")
    parser.add_option("-d",
                      type = 'string', action = 'callback', callback = callvars, dest='dates', default = None,
                      help = "pre/post-proc : Specify a date to work on (PREP files on that date must exist). post-proc :must be a subset of conf.assimdates.")
    parser.add_option("--kind",
                      action="store", type="string", dest="kind", default='beaufixpp',
                      help='kind of the run run to perform (beaufixpp for CrocOpp')
    parser.add_option("--vapp",
                      action="store", type="string", dest="vapp", default='s2m',
                      help="vapp of the soda-vortex task")
    parser.add_option("--vconf",
                      action="store", type="string", dest="vconf", default='12',
                      help="vconf of the soda-vortex task")
    parser.add_option("--todo",
                      action = 'store', type = 'string', dest='todo', default = None,
                      help= 'type of evaluation')
    parser.add_option("--nmembers",
                      action = 'store', type = int, dest='nmembers', default = None,
                      help= 'specify the number of members.')
    parser.add_option("--vars", type = 'string', action="callback", callback=callvars, default = 'all',
                      help="specify assimilation variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snowdepth")
    parser.add_option("--ppvars", type = 'string', action="callback", callback=callvars, default = [],
                      help = 'specify study variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snowdepth")')
    parser.add_option("--fact", type = 'string', action = 'callback', callback = callvars, default = 1.,
                      help = ' set multiplicative factors for the observation errors')
    parser.add_option("--classesS", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed slope classes ex : 0,20,40")
    parser.add_option("--classesA", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed aspect classes ex : N,NE,E,SE,S,SO")
    parser.add_option("--classesE", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed elev classes ex : 1800,2100,2400,2700")
    parser.add_option("--classesId", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed ids of classes separated by commas ex : 145,146,147")
    parser.add_option("-o",
                      action = 'store', type = 'string', dest='saverep', default = "",
                      help= 'name of the postproc/type_of_anal/ without /')
    parser.add_option("--pf",
                      action = 'store', dest = 'pf', default = 'global', type = 'string',
                      help='choose pf algorithm : global, rlocal, klocal')
    parser.add_option("--nloc_pf", action = 'store', dest = 'nloc_pf', default = 0, type = int,
                      help='set the localization radius (1->+oo) or the k-localization counter (1->+oo)')
    parser.add_option("--neff", action = 'store', dest = 'neff', default = 0, type = int,
                      help='set the Neff target value for the inflation factor')
    parser.add_option('--synth',
                      action = 'store', dest = 'synth', type = int, default = None,
                      help ='choose member as a synthetic observation (for local tests only, or generation of observations.')
    parser.add_option('--noise',
                      action = 'store', dest = 'noise', type = float, default = 0.,
                      help ='multiplicativesas gaussian noise to the synthetical obs.')
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


def set_options(args, readConf = True, pathConf = None, pathPgd = None ):
    """
    This methods allows to parse, format and check the options of the crocO command and merge them with those written in the configuration file (pathConf, if any).
    The user-provided options with supersede the values written in the configuration file.
    The returned object has immutable attributes.
    - readConf : load a conf file, naively looking for a conf file in the xpiddir.
    - path conf : help finding the conf file.

    """

    options = parse_options(args)
    if 'VORTEXPATH' not in list(os.environ.keys()):
        raise Exception('you must export VORTEXPATH to the root of your local vortex xps.')
    else:
        options.vortexpath = os.environ['VORTEXPATH']
    try:
        options.xpiddir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpid + '/'
    except AttributeError:
        print('you must prescribe an xpid (root of your experiment), located at $VORTEXPATH/s2m/<geometry>/xpid/ and create it before the simulation')
    if options.xpidobs is not None:
        options.xpidobsdir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/obs/' + options.xpidobs + '/'
    if options.xpidol is not None:
        options.xpidoldir = options.vortexpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpidol + '/'

    # load the PGD (mandatory), it is key to describe the working geometry.
    if pathPgd is None:
        try:
            options.pathPgd = '/'.join([os.environ['VORTEXPATH'], options.vapp, options.vconf, 'spinup/pgd/PGD_'] + area(options.vconf) + '.nc')
            options.pgd = Pgd(options.pathPgd)
        except FileNotFoundError:
            raise Exception('I could not find the PGD in the spinup dir... \nhelp me with pathPgd=<path to conf file>')
    else:
        options.pathPgd = pathPgd
        options.pgd = Pgd(options.pathPgd)
    # convert the classes*args into more explicit variables
    options = parse_classes(options)
    # load the conf file (not mandatory)
    if readConf is True:
        if pathConf is None:
            try:
                confPath = options.xpiddir + '/conf/' + options.vapp + '_' + options.vconf + '.ini'
                print('opening conf file : ', confPath)
                conf = read_conf(confPath)
            except FileNotFoundError:
                raise Exception('I could not find the conf file by myself... \nhelp me with pathConf=<path to conf file>')
        else:
            conf = read_conf(pathConf)
        # convert the classes*args into more explicit variables
        conf = parse_classes(conf)
        # ensure assimdates is a list of string
        # assimdates is always set.
        if type(conf.assimdates) is str:
            conf.assimdates = [str(conf.assimdates)]
        else:
            conf.assimdates = list(map(str, conf.assimdates))
        if hasattr(conf, 'stopdates'):
            if type(conf.stopdates) is str:
                conf.stopdates = [str(conf.stopdates)]
            else:
                conf.stopdates = list(map(str, conf.stopdates))

        # in the case we are pping, and options.dates =='all', we must set options.dates to the assimdates of the conf file
        if options.dates == 'all':
            if hasattr(conf, 'stopdates'):
                options.dates = conf.stopdates[0:-1]
            else:
                options.dates = conf.assimdates[0:-1]
        # set some important parameters:
        if not hasattr(conf, 'openloop'):
            conf.openloop = 'off'
        # BC April 2020 not sure that's safe
        if not hasattr(options, 'openloop'):
            options.openloop = 'off'

        # set the observations dir
        if options.sensor is None and options.openloop == 'off':
            try:
                mb = 'mb{0:04d}'.format(options.synth)
                options.sensor = mb + '_v' + ''.join(options.vars)  \
                    + '_E' + ''.join(options.classesE) + '_A' + ''.join(options.classesA) \
                    + '_S' + ''.join(options.classesS) + '_N' + str(options.noise)
            except TypeError:
                raise Exception('if you dont specify obs, please specify a synth member to draw')
        elif options.sensor is not None:
            if hasattr(conf, 'sensor'):
                print('caution : you\'re overwriting conf.sensor with options.sensor.')

        # for local runs : if options.synth is not None, need to decrease the number of members by 1
        # if options.synth is set to 0, draw a member:
        if options.synth is not None:
            # draw a synthetical member
            if options.synth == 0:
                options.synth = random.choice(list(range(1, options.nmembers + 1)))
                print('drawn synthetic member :', options.synth)
            # reduce nmembers by 1
            options.nmembers -= 1

        # merge options and conf into an immutable object.
        # conf file values are overwritten by options
        # do not overwrite with None:
        opt_dict = options.__dict__
        for key, val in opt_dict.items():
            if val is None:
                del opt_dict[key]
        immutable_option = ImmutableOpt(**{**conf.__dict__, **opt_dict})
    else:
        print('not reading any conf file')
        if type(options.dates) is str:
            options.dates = [options.dates]
        # return an immutable object
        immutable_option = ImmutableOpt(**options.__dict__)

    return immutable_option


def execute(args):
    options = set_options(args)
    run = CrocOrun(options)
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
