# !/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on 5 feb. 2019

@author: cluzetb

* perform local tests/dev of CROCO PF based on :
    - set of prep files from an openloop
    - opt : set of REAL observations OR generate synthetical observations.
* post-process outputs of CROCO assimilation sequence
* lauch local parallelized CROCO assimilation sequence
'''
from argparse import ArgumentParser, Namespace
import datetime
import os
import random
import sys
import time

from CrocoPf import CrocoPf, CrocoObs
from CrocoPp import CrocoPp
import numpy as np
from utilcrocO import Opt, ImmutableOpt, Pgd, area, parse_classes,\
    read_opts_in_namelist, set_sensor, set_provars, merge_two_dicts
from tools.read_conf import read_conf

# from optparse import OptionParser, Values
usage = 'crocO --opts'


def parse_args(arguments):

    parser = ArgumentParser(usage)

    parser.add_argument("--xpid",
                        type=str, action="store", dest="xpid", default=None,
                        help="xpid of the crocO run : must exist (if you are not post-processing, you must create it before launching the experiment")
    parser.add_argument("--xpidol",
                        type=str, action="store", dest="xpidol", default=None,
                        help="xpid of the associated openloop crocO run (if any)")
    parser.add_argument("-d",
                        type=callvars, dest='dates', default = None,
                        help = "pre/post-proc : Specify a date to work on (PREP files on that date must exist). post-proc :must be a subset of conf.assimdates.")
    parser.add_argument("--vapp",
                        action="store", type=str, dest="vapp", default='s2m',
                        help="vapp (model) of the soda-vortex task")
    parser.add_argument("--vconf",
                        action="store", type=str, dest="vconf", default='12',
                        help="vconf (geometry) of the soda-vortex task")
    parser.add_argument("--todo",
                        action = 'store', type = str, dest='todo', default = None,
                        choices=('generobs', 'localrun', 'pfpp', 'parallelpp', 'pf', 'pfpython', 'parallel'),
                        help= 'type of evaluation. pfpp : post-process local run of the pf (previously ran --todo {pf,pfpython}). parallelpp: load post-processing done on beaufix.')
    parser.add_argument("--nmembers",
                        action = 'store', type = int, dest='nmembers', default = None,
                        help= 'specify the number of members.')
    parser.add_argument("--vars", type=callvars, default = 'all',
                        help="specify assimilation variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snow depth, SWE for snow water equivalent")
    parser.add_argument("--ppvars", type=callvars, default = [],
                        help = 'specify study variables separated by commas : B1,...,B7 for MODIS bands, SCF, DEP for snow depth, SWE for snow water equivalent')
    parser.add_argument("--fact", type=callvars, default = 1.,
                        help = ' set multiplicative factors for the observation errors')
    parser.add_argument("--classes_s", type=callvars, default = None,
                        help="specify analyzed slope classes ex : 0,20,40")
    parser.add_argument("--classes_a", type=callvars, default = None,
                        help="specify analyzed aspect classes ex : N,NE,E,SE,S,SO")
    parser.add_argument("--classes_e", type=callvars, default = None,
                        help="specify analyzed elev classes ex : 1800,2100,2400,2700")
    parser.add_argument("--classes_id", type=callvars, default = None,
                        help="specify analyzed PYTHON (starting at 0) ids of classes separated by commas ex : 145,146,147. Exclusive with classes_e,a,s")
    parser.add_argument("-o",
                        action = 'store', type = str, dest='saverep', default = None,
                        help= 'name of the postproc/type_of_anal/ without /')
    parser.add_argument("--pf",
                        action = 'store', dest = 'pf', default = None, type = str,
                        choices=('global', 'klocal', 'rlocal', 'ol'),
                        help='choose pf algorithm : global, rlocal, klocal.')
    parser.add_argument("--xdloc_pf", type=callunits, default = None,
                        help = "specify the localisation angle (default unit in degrees) for the rlocal pf." +
                        " if not default specify the units after a semicolon e.g. : 1.9:m, 1.9:km, 1.9:deg, 1.9:rad")
    parser.add_argument("--lloo_pf",
                        action = "store_true", dest = 'lloo_pf', default = False,
                        help = "activate leave-one-out experiments for the localised PFs.")
    parser.add_argument("--neff", action = 'store', dest = 'neff', default = 0, type = int,
                        help='set the Neff target value for the inflation factor. 1 to deactivate localization.')
    parser.add_argument('--synth',
                        action = 'store', dest = 'synth', type = int, default = None,
                        help ='choose member as a synthetic observation (for local tests, generation of observations or synthetic data assimilation')
    parser.add_argument('--noise',
                        action = 'store', dest = 'noise', type = float, default = 0.,
                        help ='perturb (synthetic) observations with a multiplicative gaussian noise')
    parser.add_argument('--sensor',
                        action = 'store', dest = 'sensor', type = str, default = None,
                        help ='provide sensor name for OBS reading (MODIS, pleiades...).')
    parser.add_argument('--sensor_in',
                        action = 'store', dest = 'sensor_in', type = str, default = None,
                        help ='When masking obs from an existing sensor, provide it as an argument')
    parser.add_argument("--archive_synth",
                        action = 'store_true', dest = 'archive_synth',
                        default = True,  # BC 11/06/20 @TODO : fix local pf pp bugs with archive_synth=False
                        help = 'specify if you want to archive the synthetical obs with its mask or not.')
    parser.add_argument("--no_need_masking",
                        action = 'store_true', dest = 'no_need_masking', default = False,
                        help = 'specify if you do not need to mask the obs')
    parser.add_argument("--readprep",
                        action = 'store_true', dest = 'readprep', default = False,
                        help = 'read prep files or the pickles {bg,an,ol}_date.pkl')
    parser.add_argument("--readaux",
                        action = 'store_true', dest = 'readaux', default = False,
                        help = 'read auxiliary files')
    parser.add_argument("--readobs",
                        action = 'store_true', dest = 'readobs', default = False,
                        help = 'read observation files')
    parser.add_argument("--readtruth",
                        action = 'store_true', dest = 'readtruth', default = False,
                        help = 'read truth from openloop run (synthetic experiments only. synth must be defined in either opts or conf file.')
    parser.add_argument("--readoper",
                        action = 'store_true', dest = 'readoper', default = False,
                        help = 'read oper pickle file @TODO remove from master branch.')
    parser.add_argument("--notreadpro",
                        action = 'store_true', dest = 'notreadpro', default = False,
                        help = 'do NOT read pro files or the pickles EnsPro*')
    parser.add_argument("--clim",
                        action = 'store_true', dest = 'clim', default = False,
                        help = 'read the clim')
    # new opts for parallel run
    parser.add_argument("--arch",
                        type=str, action="store", dest="arch", default=None,
                        help=" absolute path to the archive (default=xpid)")
    parser.add_argument("-n",
                        type = str, dest = 'namelist',
                        default = None)
    parser.add_argument("-b", "--begin",
                        action="store", type=str, dest="datedeb", default=None,
                        help="Date to start the simulation (YYYYMMDD): MANDATORY OPTION")
    parser.add_argument("-e", "--end",
                        action="store", type=str, dest="datefin", default=None,
                        help="Date to finish the simulation (YYYYMMDD): MANDATORY OPTION (unless --oper)")
    parser.add_argument("-f", "--forcing",
                        action="store", type=str, dest="forcing", default=None,
                        help="path of the forcing file or of the directory with the forcing files - default: None")
    parser.add_argument("--escroc",
                        action="store", type=str, dest="escroc", default=None,
                        help="ESCROC subensemble")
    parser.add_argument("--spinup",
                        action='store', type=str, dest = 'spinup', default=None,
                        help='path to the pgd and prep files')
    parser.add_argument("--provars",
                        type=callvars, dest = 'provars', default = ['all_notartes'],
                        help = 'specify the list of variables to write down into the PRO files. (CSELECT in the namelist\n\
                      Default : TALB_ISBA (albedo), TS_ISBA (Snow Surface temperature), DSN_T_ISBA (Snow Depth) and WSN_T_ISBA (Snow Water Equivalent).')
    parser.add_argument("--pathConf",
                        action="store", type=str, dest="pathConf", default=None,
                        help="specify an unusual pathConf")

    opts_no_defaults = Namespace()
    args = parser.parse_args( args=arguments,
                              namespace=opts_no_defaults)
    # BC june 2020 mimicking parser.get_default_values
    # options = Values(parser.get_default_values().__dict__)
    options = Namespace(**{key: parser.get_default(key) for key in args.__dict__})
    # BC june 2020 mimicking Values._update_careful
    # options._update_careful(opts_no_defaults.__dict__)
    for key in options.__dict__.keys():
        if key in opts_no_defaults.__dict__.keys():
            dval = opts_no_defaults.__dict__[key]
            if dval is not None:
                options.__setattr__(key, dval)
    del args

    return options, opts_no_defaults


def callvars(value):
    if type(value) is not float:  # normal case
        return value.split(',')
    else:  # bug with float args like --fact
        return value


def callunits(value):
    """
    BC 09/06/20 : units converter for xdloc_pf.
    """
    print(value)
    if type(value) is not float:
        ll = value.split(':')
        ll[0] = float(ll[0])
        if ll[1] == 'deg':
            ret = ll[0]
        elif ll[1] == 'm':
            ret = ll[0] / 6371000.
        elif ll[1] == 'km':
            ret = ll[0] / 6371. * 180. / np.pi
        elif ll[1] == 'rad':
            ret = ll[0] * 180 / np.pi
        else:
            raise Exception('supported units :deg, m, km, rad.')
        return ret
    else:  # normal case
        return value


def set_options(args, readConf = True, useVortex = True, pathConf = None, pathPgd = None, mutable=False):
    """
    This methods allows to parse, format and check the options of the crocO command and merge them with those written in the configuration file (pathConf, if any).
    The user-provided options with supersede the values written in the configuration file which supersed over the default not in args options.
    The returned object has immutable attributes.
    - readConf : load a conf file, naively looking for a conf file in the xpiddir.
    - pathConf : help finding the conf file.
    - pathPgd : help finding the conf file.
    - mutable : NOT RECOMMENDED: if True, return a mutable option object (lazy useful solution for pp on beaufix)
    """
    options, no_default_opts = parse_args(args[1:])
    if 'CROCOPATH' not in list(os.environ.keys()):
        raise Exception('you must export CROCOPATH to the root of your local experiments.')
    else:
        options.crocOpath = os.environ['CROCOPATH']
    try:
        if not os.path.isabs(options.xpid):
            options.xpiddir = options.crocOpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpid + '/'
        else:
            options.xpiddir = options.xpid
    except AttributeError:
        print('you must prescribe an xpid (root of your experiment), located at $CROCOPATH/s2m/<geometry>/xpid/ and create it before the simulation')

    if (options.classes_e is not None or options.classes_a is not None or options.classes_s is not None) and options.classes_id is not None:
        raise Exception('classes_id is exclusive with the other classes args.')
    if options.sensor is not None:
        options.sensordir = options.crocOpath + '/' + options.vapp + '/' + options.vconf + '/obs/' + options.sensor + '/'
    if options.sensor_in is not None:
        options.sensor_in_dir = options.crocOpath + '/' + options.vapp + '/' + options.vconf + '/obs/' + options.sensor_in + '/'
    if options.xpidol is not None:
        options.xpidoldir = options.crocOpath + '/' + options.vapp + '/' + options.vconf + '/' + options.xpidol + '/'
    if options.namelist is None:
        options.namelist = os.environ['SNOWTOOLS_CEN'] + '/DATA/OPTIONS_V81_NEW_OUTPUTS_NC.nam'
        print("forcing the namelist path to its default value :", options.namelist)
        options.namelist_is_default = True
    else:
        options.namelist_is_default = False
    if options.todo != 'generobs' and options.synth and options.sensor:
        raise Exception(" ither you are running a synth experiment from scratch (--synth) or you are using pre-generated observations (--sensor), but not both !")
    # load the PGD (mandatory), it is key to describe the working geometry.
    if pathPgd is None:
        try:
            options.pathPgd = '/'.join([os.environ['CROCOPATH'], options.vapp, options.vconf, 'spinup/pgd/PGD_']) + area(options.vconf) + '.nc'
            options.pgd = Pgd(options.pathPgd)
        except IOError:
            raise Exception('I could not find the PGD in the spinup dir.',
                            'help me with pathPgd=<path to conf file>')
    else:
        options.pathPgd = pathPgd
        options.pgd = Pgd(options.pathPgd)
    # convert the classes args into more explicit variables
    options = parse_classes(options)
    # convert the provars arg
    options = set_provars(options)
    if pathConf:
        # readConf is useful only when the path to the conf is implicit.
        readConf = True
    elif options.pathConf:
        pathConf = options.pathConf
        readConf = True
    # load the conf file (not mandatory)
    if readConf is True:
        if pathConf is None:
            try:
                confPath = options.xpiddir + '/conf/' + options.vapp + '_' + options.vconf + '.ini'
                print('opening conf file : ', confPath)
                conf = read_conf(confPath, useVortex = useVortex)
            except IOError:
                raise Exception('I could not find the conf file by myself', 'help me with pathConf=<path to conf file>')
        else:
            conf = read_conf(pathConf, useVortex=useVortex)
        # convert the classes*args into more explicit variables
        conf.pgd = options.pgd
        conf = parse_classes(conf)

        # populate the conf file with the default options THAT WHERE NOT provided by the user.
        # the conf file supersedes over these not provided default options
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
        conf = Opt(**merge_two_dicts(updict, conf_dict))

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
            # conf.stopdates is useful in parallel mode
            if not hasattr(conf, 'stopdates'):
                if options.datefin:
                    conf.stopdates = conf.assimdates + [options.datefin]
                else:
                    try:
                        test = datetime.datetime.strptime(conf.datefin, '%Y%m%d%H')
                    except ValueError:
                        test = datetime.datetime.strptime(conf.datefin, "%Y-%m-%d %H:%M:%S")
                        print('there should be a bug when reading the conf file with vortex. Not handled a the moment')
                    except (TypeError, AttributeError) as _:
                        print('no datefin in conf file. probably generating obs.')
                        test = None
                    conf.stopdates = conf.assimdates + [test] if test else conf.assimdates

            options.dates = conf.assimdates
        # set openloop or not.
        if not hasattr(conf, 'openloop'):
            conf.openloop = 'off'
        # in case the pf arg has not been provided, or in case of local test of the pf,
        # do not overwrite conf.openloop
        if options.todo == 'parallelpp':
            options.openloop = conf.openloop
            options = read_opts_in_namelist(options)
        elif options.pf is None or options.todo == 'pf':
            options.openloop = conf.openloop
        elif options.pf == 'ol':
            options.openloop = 'on'
        else:
            options.openloop = 'off'
        # set the observations dir
        if options.sensor is None:
            if options.openloop == 'off':
                if hasattr(conf, 'sensor'):
                    if conf.sensor is None:
                        try:
                            options = set_sensor(options)
                        except TypeError:
                            raise Exception
                    else:
                        print("reading the sensor from conf.sensor", conf.sensor)
                else:
                    if options.readobs or options.readtruth:
                        raise Exception('No sensor found in the conf file',
                                        'please specify a sensor either in options or the conf file')
            else:
                try:
                    options = set_sensor(options)
                except (TypeError, ValueError) as _:
                    pass
        else:
            if hasattr(conf, 'sensor'):
                print('caution : you\'re overwriting conf.sensor with options.sensor.')

        # for local pf runs : if options.synth is not None, need to decrease the number of members by 1
        # if options.synth is set to 0, draw a member:
        if not options.nmembers:
            options.nmembers = conf.nmembers
        options.mblist = list(range(1, int(options.nmembers) + 1))
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
            if key in conf.__dict__ and opt_dict[key] is None:
                del opt_dict[key]
        if mutable is False:
            return_option = ImmutableOpt(**merge_two_dicts(conf.__dict__, opt_dict))
        else:
            return_option = Opt(**merge_two_dicts(conf.__dict__, opt_dict))

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
        if mutable is False:
            # return an immutable object
            return_option = ImmutableOpt(**options.__dict__)
        else:
            return_option = Opt(**options.__dict__)
    return return_option


def execute(args):
    """
    Use of crocO for pre/post-processing through the command line is possible but not handy.
    Basic sequences of execution are presented below.
    The recommended use of crocO is to launch it in a separate script/notebook (see some examples in pproc_scripts/ :
     - set the args
     - parse them with set_options.
     If you're post-processing some experiment, the conf file of the experiment will be read and overwritten by your options.
     - instantiate CrocoPf or CrocoPp
     - have fun with those handy objects.
    """

    options = set_options(args)
    if options.todo in ['pfpp', 'parallelpp']:
        run = CrocoPp(options)
    elif options.todo in ['pf', 'pfpython']:
        run = CrocoPf(options)
        run.run()
        pp = CrocoPp(options)
    elif options.todo == 'generobs':
        run = CrocoObs(options)

    return run


if __name__ == '__main__':
    start_time = time.time()

    run = execute(sys.argv)
    elapsed_time = time.time() - start_time
    print(elapsed_time)
