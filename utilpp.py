#! /usr/bin/python
# -*- coding: utf-8 -*-

'''
Created on 3 juin 2019

@author: cluzetb
'''

import calendar
import copy
import datetime
import pickle

import numpy as np


def read_part(options):
    gg3 = dict()
    for dd in options.dates:
        if options.todo == 'pfpp' or options.todo == 'pf':
            filename = dd + '/PART'
            f = open(filename, 'rb')
        else:  # =='parallelpp'
            try:
                filename = options.xpiddir + '/workSODA/PART_' + dd + '.txt.foo'
                f = open(filename, 'rb')
            except FileNotFoundError:
                filename = options.xpiddir + '/workSODA/PART_' + dd + '.txt'
                f = open(filename, 'rb')
        gg3[dd] = np.genfromtxt(f,
                                delimiter = ',',
                                usecols = list(range(0, options.nmembers)))
    return gg3


def read_mask(options):
    npts = 187
    imask = dict()
    for dd in options.dates:
        if options.todo == 'pfpp' or options.todo == 'pf':
            filename = dd + '/IMASK'
            f = open(filename, 'r')

        else:  # =='parallelpp'
            try:
                filename = options.xpiddir + '/workSODA/IMASK_' + dd + '.txt.foo'
                f = open(filename, 'rb')
            except FileNotFoundError:
                filename = options.xpiddir + '/workSODA/IMASK_' + dd + '.txt'
                f = open(filename, 'rb')

        imask[dd] = dict()
        for var in options.vars:
            imask[dd][var] = np.empty((npts, options.nloc_pf), dtype = int)
        il = 0
        il4 = 0
        il5 = 0
        for line in f:
            data = np.array(line.split(',')[0:-1])
            if il % 2 == 0:
                imask[dd]['B4'][il4, :] = data  # convert into PYTHON class indices
                imask[dd]['B4'][il4, :] = imask[dd]['B4'][il4, :] - 1.
                il4 += 1
            else:
                imask[dd]['B5'][il5, :] = data  # convert into PYTHON class indices
                imask[dd]['B5'][il5, :] = imask[dd]['B5'][il5, :] - 1.
                il5 += 1
            il += 1
        imask[dd]['B4'] = np.ma.masked_where(imask[dd]['B4'] == -1., imask[dd]['B4'] )
        imask[dd]['B5'] = np.ma.masked_where(imask[dd]['B5'] == -1., imask[dd]['B5'] )
    return imask


def read_BG(options):
    """
    be careful, vars in options.vars MUST be put in the same order they are written in the file, i.e. the order in the assim namelist.
    """

    bg = dict()
    for dd in options.dates:
        if options.todo == 'pfpp' or options.todo == 'pf':
            filename = dd + '/BG_CORR'
            f = open(filename, 'r')
        else:  # =='parallelpp'
            try:
                filename = options.xpiddir + '/workSODA/BG_CORR_' + dd + '.txt.foo'
                f = open(filename, 'rb')
            except FileNotFoundError:
                filename = options.xpiddir + '/workSODA/BG_CORR_' + dd + '.txt'
                f = open(filename, 'rb')

        bg[dd] = dict()
        data = np.array([[float(v) for v in line.split(',')[0:-1]] for line in f])
        for iv, var in enumerate(options.vars):
            bg[dd][var] = data[iv::len(options.vars)]
            bg[dd][var] = np.ma.masked_where(bg[dd][var] == 1.e+20, bg[dd][var])
    return bg


def read_alpha(options):
    alpha = dict()
    for dd in options.dates:
        if options.todo == 'pfpp' or options.todo == 'pf' or options.todo == 'pf':
            filename = dd + '/ALPHA'
            f = open(filename, 'r')
        else:  # =='parallelpp'
            try:
                filename = options.xpiddir + '/workSODA/ALPHA_' + dd + '.txt.foo'
                f = open(filename, 'r')
            except FileNotFoundError:
                filename = options.xpiddir + '/workSODA/ALPHA_' + dd + '.txt'
                f = open(filename, 'r')
        for line in f:
            if ',' in line:
                alpha[dd] = np.array([float(a) for a in line.split(',')[0:-1]])
            else:
                alpha[dd] = float(line)

    return alpha


def set_itimes(run, clim =False, fromOl = False):
    '''
    return a mask of valid dates (after sept of starting year + only midday)
    '''
    # be careful sometimes times differ btw openloop and an
    if hasattr(run, 'ensProOl') and not hasattr(run, 'ensProAn') or fromOl:
        ttttimes = run.ensProOl['time']
    else:
        ttttimes = run.ensProAn['time']
    times = [not((t.year == run.datedeb.year and t.month <= 9) or t.hour != 0 or (t.year == run.datefin.year and t.month > 6))
             for t in ttttimes]
    itimes = [i for i, t in enumerate(times) if t]

    if clim is False:
        return itimes
    else:
        if calendar.isleap(run.datefin.year):
            start_year = 2015
        else:
            start_year = 2001

        ttC = [datetime.datetime(start_year, 8, 1) + datetime.timedelta(d) for d in range(0, 365 if start_year == 2001 else 366)]
        timesC = [not((t.year == start_year and t.month <= 9) or (t.year == start_year + 1 and t.month > 6)) for t in ttC]
        itimesC = [i for i, t in enumerate(timesC) if t]
        return itimes, itimesC


def load_pickle2(pathPkl):
    '''
    BC Feb 2020
    pickle files can be generated by python2 on beaufix.
    When reading it with python 3.7.2+, use encoding = 'latin1' for compatibility
    https://stackoverflow.com/questions/28218466/unpickling-a-python-2-object-with-python-3
    '''
    with open(pathPkl, 'rb') as f:
        print(('loading from pickle ! ', pathPkl))
        gg = pickle.load(f, encoding = 'latin1')
    return gg


def old_read_truth(run, var, baseline = False):
    """
    BC 24/10/19 : function to read the truth in the OL from the corresponding run.
    - Sometimes, if the OL is a subset of a bigger OL, and the truth is in this bigger OL, need to read in it.
    k is the variable.
    - find the baseline experiment if necessary.
    """
    if hasattr(run, 'truth'):
        return run.truth
    if not baseline:
        try:
            truth = run.ensProOl[var][:, :, run.mbsynth]
        except Exception:
            print('\n\nWARNING : there is no corresponding mbsynth in this openloop not enough members\n looking in the bigger OL xp.\n\n')
            print('loading ' + run.xpidoldir[0:-4] + '/crocO/' + run.options.saverep + '/ensProOl.pkl')
            with open(run.xpidoldir[0:-4] + '/crocO/' + run.options.saverep + '/ensProOl.pkl', 'rb') as f:
                gg = pickle.load(f)
            truth = gg[var][:, :, run.mbsynth]
        return truth
    else:
        import CrocoPp

        # bc do not use prev. block because need to set itimes
        opts = copy.copy(run.options)
        print('doesnot work, read in file')
        opts.xpid = 'baseline_{0}'.format(run.assimdates[0].strftime('%Y')) + '/'
        opts.xpiddir = opts.crocOpath + '/' + opts.vapp + '/' + opts.vconf + '/' + opts.xpid
        opts.nmembers = 1
        base = CrocoPp.CrocoPp(opts, pathConf='{0}baseline_{1}/conf/s2m_12.ini'.format(run.rootdir, run.assimdates[0].strftime('%Y')))
        truth = base.ensProOl[var]
        itimes = set_itimes(base, fromOl = True)
        if truth.shape[-1] > 1:
            raise Exception('youre not readin a truth.')
        return truth, itimes


def RMSE(ens, truth, aggrTime = False, aggrDomain = True):
    """
    (aggrtime = False, aggrDomain = False): time-variant RMSE of an ensemble median over a domain.
    (aggrtime = True, aggrDomain = False): domain variant, time averaged RMSE

    IN: ens (ndate, npts, nmembers)
        truth    (ndate, npts)
    OUT: RMSE(ndate)
    """
    if aggrDomain:
        if not aggrTime:
            return np.sqrt(np.mean(np.square(np.mean(ens, axis = 2) - truth), axis = 1))
        else:
            return np.sqrt(np.mean(np.mean(np.square(np.mean(ens, axis = 2) - truth), axis = 1)))
    else:
        if not aggrTime:
            return np.sqrt(np.square(np.mean(ens, axis = 2) - truth))
        else:
            return np.sqrt(np.mean(np.square(np.mean(ens, axis = 2) - truth), axis = 0))


def spread(ens, aggrTime = False, aggrDomain = True):
    """
    time-variant spread of an ensemble over a domain.
    IN: ENSEMBLE (ndate, npts, nmembers)
    OUT: spread(ndate)
    """
    if aggrDomain:
        if not aggrTime:
            return np.sqrt(np.mean(np.array(
                [np.mean((m - np.expand_dims(meanPt, axis = 1))**2, axis = 1)
                 for (m, meanPt) in zip(np.rollaxis(ens, 1), np.mean(ens, axis = 2).T)]
            ), axis = 0))
        else:
            return np.sqrt(np.mean(np.mean(np.array(
                [np.mean((m - np.expand_dims(meanPt, axis = 1))**2, axis = 1)
                 for (m, meanPt) in zip(np.rollaxis(ens, 1), np.mean(ens, axis = 2).T)]
            ), axis = 0)))
    else:
        if not aggrTime:
            return np.sqrt(np.array(
                [np.mean((m - np.expand_dims(meanPt, axis = 1))**2, axis = 1)
                 for (m, meanPt) in zip(np.rollaxis(ens, 1), np.mean(ens, axis = 2).T)]
            )).T
        else:
            return np.sqrt(np.array(np.mean(
                [np.mean((m - np.expand_dims(meanPt, axis = 1))**2, axis = 1)
                 for (m, meanPt) in zip(np.rollaxis(ens, 1), np.mean(ens, axis = 2).T)], axis=1))).T
