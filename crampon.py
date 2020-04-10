# -*- coding: utf-8 -*-
'''
Created on 5 feb. 2019

@author: cluzetb

* perform local tests/dev of CRAMPON PF based on :
    - set of prep files from an openloop
    - opt : set of REAL observations OR generate synthetical observations.
* post-process outputs of CRAMPON assimilation sequence
* lauch local parallelized CRAMPON assimilation sequence
'''
from CramponPf import CramponPf, CramponObs
from CramponPp import CramponPp
from optparse import OptionParser, Values
import os
import random
import sys
import time
from utilcrampon import read_conf, Opt, ImmutableOpt, Pgd, area, parse_classes,\
    read_opts_in_namelist


usage = 'crampon --opts'


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
    parser.add_option("--vapp",
                      action="store", type="string", dest="vapp", default='s2m',
                      help="vapp (model) of the soda-vortex task")
    parser.add_option("--vconf",
                      action="store", type="string", dest="vconf", default='12',
                      help="vconf (geometry) of the soda-vortex task")
    parser.add_option("--todo",
                      action = 'store', type = 'choice', dest='todo', default = None,
                      choices=('generobs', 'localrun', 'localpp', 'beaufixpp', 'pf', 'pfpython'),
                      help= 'type of evaluation. localpp : post-process local run of the pf (previously ran --todo {pf,pfpython}). beaufixpp: load post-processing done on beaufix.')
    parser.add_option("--nmembers",
                      action = 'store', type = int, dest='nmembers', default = None,
                      help= 'specify the number of members.')
    parser.add_option("--vars", type = 'string', action="callback", callback=callvars, default = 'all',
                      help="specify assimilation variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snow depth, SWE for snow water equivalent")
    parser.add_option("--ppvars", type = 'string', action="callback", callback=callvars, default = [],
                      help = 'specify study variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snow depth, SWE for snow water equivalent')
    parser.add_option("--fact", type = 'string', action = 'callback', callback = callvars, default = 1.,
                      help = ' set multiplicative factors for the observation errors')
    parser.add_option("--classes_s", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed slope classes ex : 0,20,40")
    parser.add_option("--classes_a", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed aspect classes ex : N,NE,E,SE,S,SO")
    parser.add_option("--classes_e", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed elev classes ex : 1800,2100,2400,2700")
    parser.add_option("--classes_id", type = 'string', action="callback", callback=callvars, default = None,
                      help="specify analyzed ids of classes separated by commas ex : 145,146,147")
    parser.add_option("-o",
                      action = 'store', type = 'string', dest='saverep', default = "",
                      help= 'name of the postproc/type_of_anal/ without /')
    parser.add_option("--pf",
                      action = 'store', dest = 'pf', default = None, type = 'choice',
                      choices=('global', 'klocal', 'rlocal'),
                      help='choose pf algorithm : global, rlocal, klocal.')
    parser.add_option("--nloc_pf", action = 'store', dest = 'nloc_pf', default = None, type = int,
                      help='set the localization radius (1->+oo) or the k-localization counter (1->+oo).')
    parser.add_option("--neff", action = 'store', dest = 'neff', default = 0, type = int,
                      help='set the Neff target value for the inflation factor. 1 to deactivate localization.')
    parser.add_option('--synth',
                      action = 'store', dest = 'synth', type = int, default = None,
                      help ='choose member as a synthetic observation (for local tests only, or generation of observations.')
    parser.add_option('--noise',
                      action = 'store', dest = 'noise', type = float, default = 0.,
                      help ='perturb (synthetic) observations with a multiplicative gaussian noise')
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
                      help = 'read prep files or the pickles {bg,an,ol}_date.pkl')
    parser.add_option("--readaux",
                      action = 'store_true', dest = 'readaux', default = False,
                      help = 'read auxiliary files')
    parser.add_option("--readobs",
                      action = 'store_true', dest = 'readobs', default = False,
                      help = 'read observation files')
    parser.add_option("--readtruth",
                      action = 'store_true', dest = 'readtruth', default = False,
                      help = 'read truth from openloop run (synthetic experiments only. synth must be defined in either opts or conf file.')
    parser.add_option("--readoper",
                      action = 'store_true', dest = 'readoper', default = False,
                      help = 'read oper pickle file @TODO remove from master branch.')
    parser.add_option("--notreadpro",
                      action = 'store_true', dest = 'notreadpro', default = False,
                      help = 'do NOT read pro files or the pickles EnsPro*')
    parser.add_option("--clim",
                      action = 'store_true', dest = 'clim', default = False,
                      help = 'read the clim')
    opts_no_defaults = Values()
    __, args = parser.parse_args(arguments, values=opts_no_defaults)
    options = Values(parser.get_default_values().__dict__)
    options._update_careful(opts_no_defaults.__dict__)

    del args

    return options, opts_no_defaults


def callvars(option, opt, value, parser):
    if type(value) is not float:  # normal case
        setattr(parser.values, option.dest, value.split(','))
    else:  # bug with float args like --fact
        setattr(parser.values, option.dest, value)


def set_options(args, readConf = True, pathConf = None, pathPgd = None, mutable=False):
    """
    This methods allows to parse, format and check the options of the crampon command and merge them with those written in the configuration file (pathConf, if any).
    The user-provided options with supersede the values written in the configuration file which supersed over the default not in args options.
    The returned object has immutable attributes.
    - readConf : load a conf file, naively looking for a conf file in the xpiddir.
    - pathConf : help finding the conf file.
    - pathPgd : help finding the conf file.
    - mutable : NOT RECOMMENDED: if True, return a mutable option object (lazy useful solution for pp on beaufix)
    """
    options, no_default_opts = parse_options(args)
    if 'CRAMPONPATH' not in list(os.environ.keys()):
        raise Exception('you must export CRAMPONPATH to the root of your local experiments.')
    else:
        options.cramponpath = os.environ['CRAMPONPATH']
    try:
        options.xpiddir = options.cramponpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpid + '/'
    except AttributeError:
        print('you must prescribe an xpid (root of your experiment), located at $CRAMPONPATH/s2m/<geometry>/xpid/ and create it before the simulation')
    if options.xpidobs is not None:
        options.xpidobsdir = options.cramponpath + '/' + options.vapp + '/' + options.vconf + '/obs/' + options.xpidobs + '/'
    if options.xpidol is not None:
        options.xpidoldir = options.cramponpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpidol + '/'

    if options.todo != 'generobs' and options.synth and options.sensor:
        raise Exception(" either you are running a synth experiment from scratch (--synth) or you are using pre-generated observations (--sensor)")
    # load the PGD (mandatory), it is key to describe the working geometry.
    if pathPgd is None:
        try:
            options.pathPgd = '/'.join([os.environ['CRAMPONPATH'], options.vapp, options.vconf, 'spinup/pgd/PGD_']) + area(options.vconf) + '.nc'
            options.pgd = Pgd(options.pathPgd)
        except FileNotFoundError:
            raise Exception('I could not find the PGD in the spinup dir... \nhelp me with pathPgd=<path to conf file>')
    else:
        options.pathPgd = pathPgd
        options.pgd = Pgd(options.pathPgd)
    # convert the classes*args into more explicit variables
    options = parse_classes(options)
    if pathConf:
        # readConf is useful only when the path to the conf is implicit.
        readConf = True

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
        conf.pgd = options.pgd
        conf = parse_classes(conf)

        # populate the conf file with the default options THAT WHERE NOT provided by the user.
        # the conf file supersedes over theses not provided default options
        updict = {}
        for opt in options.__dict__:
            if opt not in no_default_opts.__dict__:
                updict[opt] = options.__dict__[opt]
        # unpack the conf dict in case the conf file was parsed using vortex
        if '_internal' in conf.__dict__:
            conf_dict = dict()
            for key in conf.__dict__['_internal']:
                conf_dict[key] = conf.__getattr__(key)
        else:
            conf_dict = conf.__dict__
        conf = Opt(**{**updict, **conf_dict})

        # set some important parameters:

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

        if options.dates and 'all' in options.dates:
            if hasattr(conf, 'stopdates'):
                options.dates = conf.stopdates[0:-1]
            else:
                options.dates = conf.assimdates[0:-1]
        # set openloop or not.
        if not hasattr(conf, 'openloop'):
            conf.openloop = 'off'
        # in case the pf arg has not been provided, or in case of local test of the pf,
        # do not overwrite conf.openloop
        elif options.pf is None or options.todo == 'pf':
            options.openloop = conf.openloop
        elif options.pf == 'ol':
            options.openloop = 'on'
        else:
            options.openloop = 'off'
        # set the observations dir
        if options.sensor is None:
            if hasattr(conf, 'sensor'):
                if conf.sensor is None:
                    try:
                        mb = 'mb{0:04d}'.format(options.synth)
                        print(' clclc', options.classes_e, options.classes_a, options.classes_s)
                        options.sensor = mb + '_v' + ''.join(options.vars)  \
                            + '_E' + ''.join(options.classes_e) + '_A' + ''.join(options.classes_a) \
                            + '_S' + ''.join(options.classes_s) + '_N' + str(options.noise)
                    except AttributeError:
                        raise Exception('if you dont specify obs, please specify a synth member to draw')
                else:
                    print("reading the obs from conf.sensor", conf.sensor)
            else:
                if options.readobs or options.readtruth:
                    raise Exception('No sensor found in the conf file',
                                    'please specify a sensor either in options or the conf file')
        else:
            if hasattr(conf, 'sensor'):
                print('caution : you\'re overwriting conf.sensor with options.sensor.')

        # for local pf runs : if options.synth is not None, need to decrease the number of members by 1
        # if options.synth is set to 0, draw a member:
        options.mblist = list(range(1, options.nmembers + 1))
        if options.synth is not None:
            # draw a synthetical member
            if options.synth == 0:
                options.synth = random.choice(options.mblist)
                print('drawn synthetic member :', options.synth)
            # reduce nmembers by 1
            if options.synth in options.mblist:
                options.mblist.remove(options.synth)
        options.nmembers = len(options.mblist)
        # merge options and conf into an immutable object.
        # conf file values are overwritten by options
        # do not overwrite with None:
        opt_dict = options.__dict__
        for key in list(opt_dict):
            if opt_dict[key] is None:
                del opt_dict[key]

        immutable_option = ImmutableOpt(**{**conf.__dict__, **opt_dict})
    else:
        print('not reading any conf file')

        # read the PF parameters in the namelist
        # only if they are not defined by the user and not in the conf file.
        # set openloop parameters
        if options.pf:
            options.openloop = 'off'
            options = read_opts_in_namelist(options)
        else:
            options.openloop = 'on'
        if type(options.dates) is str:
            options.dates = [options.dates]
        # return an immutable object
        immutable_option = ImmutableOpt(**options.__dict__)
    return immutable_option


def execute(args):
    """
    Use of crampon for pre/post-processing through the command line is possible but not handy.
    Basic sequences of execution are presented below.
    The recommended use of crampon is to launch it in a separate script/notebook (see some examples in pproc_scripts/ :
     - set the args
     - parse them with set_options.
     If you're post-processing some experiment, the conf file of the experiment will be read and overwritten by your options.
     - instantiate CramponPf or CramponPp
     - have fun with those handy objects.
    """

    options = set_options(args)
    if options.todo in ['localpp', 'beaufixpp']:
        run = CramponPp(options)
    elif options.todo in ['pf', 'pfpython']:
        run = CramponPf(options)
        pp = CramponPp(options)
    elif options.todo == 'generobs':
        run = CramponObs(options)


if __name__ == '__main__':
    start_time = time.time()

    execute(sys.argv)
    elapsed_time = time.time() - start_time
    print(elapsed_time)
