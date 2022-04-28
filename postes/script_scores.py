'''
Created on 13 mai 2020

@author: cluzetb
'''

from snowtools.scores.ensemble import EnsembleScores

import numpy as np
from utilpp import RMSE, spread, bias
from snowtools.scores.generic import rankDiagram


def compute_CRPSS_ol_oper(runs, xps, selMassif, alwaysXclude=[]):
    '''
    compute the CRPSS of each element of runs versus the openloop and the oper run
    return matrices of CRPSS with rows and cols sorted by not excluded/excluded and then by elevation
    xps is a dict of runs simplified ids.
    selMassif is a submassif of the geom (if any)
    alwaysXclude is a list of posts of this massif you excluded from the assimilation
    '''

    # iterate over classesId in order to sort by alX /not alX and then by elevation
    r0 = runs[list(runs.keys())[0]]
    pgd = r0.pgd
    year = r0.ensProAn['time'][0].year

    if not hasattr(r0, 'ensProOl') or not hasattr(r0, 'obsTs') or not hasattr(r0, 'oper'):
        raise Exception('please use --readoper, --xpidol and --readobs')

    classesId_noAlX = np.array([i for i, m in enumerate(pgd.massif) if m in selMassif and str(pgd.station[i]) not in alwaysXclude])
    classesId_AlX = np.array([i for i, m in enumerate(pgd.massif) if m in selMassif and str(pgd.station[i]) in alwaysXclude])
    classesId = np.concatenate((classesId_noAlX[np.argsort(pgd.elev[classesId_noAlX])], classesId_AlX[np.argsort(pgd.elev[classesId_AlX])]), axis = 0)
    CRPSS_ol = np.empty((len(runs.keys()), len(classesId) ))
    ReliS_ol = np.empty((len(runs.keys()), len(classesId)))
    CRPSS_oper = np.empty((len(runs.keys()), len(classesId) ))
    ReliS_oper = np.empty((len(runs.keys()), len(classesId)))
    list_enum = [str(pgd.station[c]) for c in classesId if str(pgd.station[c]) in xps.keys() and str(pgd.station[c]) not in alwaysXclude]

    if len(alwaysXclude) > 0:
        list_enum.append('alX')
    list_enum.append('all')

    for ii, station in enumerate(list_enum):

        run = runs[station]
        CRPSS_oper[ii, :], ReliS_oper[ii, :], CRPSS_ol[ii, :], ReliS_ol[ii, :] = CRPSS_ol_oper(run, year, classesId)

    return CRPSS_ol, ReliS_ol, CRPSS_oper, ReliS_oper, classesId, classesId_AlX, classesId_noAlX


def frame_year_classes(run, classesId=None ):
    """
    reduce time to the winter season and extract only dates in common
    """
    if classesId is None:
        classesId = range(run.pgd.npts)
    year = int(run.ensProAn['time'][0].strftime('%Y'))
    mobs = [i for i, t in enumerate(run.obsTs['time'])
            if (t.hour == 6 and ((t.year == year and t.month > 9) or (t.month < 7 and t.year == year + 1)))]
    mensAn = [i for i, t in enumerate(run.ensProAn['time']) if (t.hour == 6 and (t.month > 9 or t.month < 7))]
    mensOl = [i for i, t in enumerate(run.ensProOl['time']) if (t.hour == 6 and (t.month > 9 or t.month < 7))]
    moper = [i for i, t in enumerate(run.oper['time']) if (t.hour == 6 and (t.month > 9 or t.month < 7))]
    # mask out nans in obs
    obs = np.ma.masked_invalid(run.obsTs['DEP'][mobs, :][:, classesId])
    ensAn = run.ensProAn['DEP'][mensAn, :, :][:, classesId, :]
    ensOl = run.ensProOl['DEP'][mensOl, :, :][:, classesId, :]
    oper = run.oper['DEP'][moper, :][:, classesId]
    return obs, ensAn, ensOl, oper


def CRPSS_ol_oper(run, classesId=None, aggrSpace = False, aggrTime = True, relativeScore = False,
                  no_oper = False, no_decomp = False):
    if not hasattr(run, 'fake'):
        obs, ensAn, ensOl, oper = frame_year_classes(run, classesId=classesId)

    else:
        obs = run.obs
        ensAn = run.ensAn
        ensOl = run.ensOl
        oper = run.oper
    if classesId is None:
        classesId = range(run.pgd.npts)
    if aggrTime is True and aggrSpace is False:
        if not no_decomp:  # slower ?
            CRPSAn = np.array([EnsembleScores(list(range(ensAn.shape[0])),
                                              list(range(ensAn.shape[0])),
                                              obs[:, cl],
                                              ensAn[:, cl, :].T,
                                              ).CRPS_decomp() for cl in range(len(classesId))])
            CRPS_ol = np.array([EnsembleScores(list(range(ensOl.shape[0])),
                                               list(range(ensOl.shape[0])),
                                               obs[:, cl],
                                               ensOl[:, cl, :].T,
                                               ).CRPS_decomp() for cl in range(len(classesId))])
            if not no_oper:
                CRPS_oper = np.array([EnsembleScores(list(range(oper.shape[0])),
                                                     list(range(oper.shape[0])),
                                                     obs[:, cl],
                                                     oper[:, cl, :].T,
                                                     ).CRPS_decomp() for cl in range(len(classesId))])
        else:
            CRPSAn = np.array([EnsembleScores(list(range(ensAn.shape[0])),
                                              list(range(ensAn.shape[0])),
                                              obs[:, cl],
                                              ensAn[:, cl, :].T,
                                              ).CRPS() for cl in range(len(classesId))])
            print()
            CRPS_ol = np.array([EnsembleScores(list(range(ensOl.shape[0])),
                                               list(range(ensOl.shape[0])),
                                               obs[:, cl],
                                               ensOl[:, cl, :].T,
                                               ).CRPS() for cl in range(len(classesId))])
            if not no_oper:
                CRPS_oper = np.array([EnsembleScores(list(range(oper.shape[0])),
                                                     list(range(oper.shape[0])),
                                                     obs[:, cl],
                                                     oper[:, cl, :].T,
                                                     ).CRPS() for cl in range(len(classesId))])

        if relativeScore is False:
            if not no_decomp:
                if not no_oper:
                    CRPSS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_oper[i, 0] if CRPS_oper[i, 0] > 0. else np.nan
                                                                for i in range(CRPSAn.shape[0])]))
                    ReliS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 1] / CRPS_oper[i, 1] if CRPS_oper[i, 1] > 0. else np.nan
                                                                for i in range(CRPSAn.shape[0])]))
                CRPSS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_ol[i, 0] if CRPS_ol[i, 0] > 0. else np.nan
                                                          for i in range(CRPSAn.shape[0])]))
                ReliS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 1] / CRPS_ol[i, 1] if CRPS_ol[i, 1] > 0. else np.nan
                                                          for i in range(CRPSAn.shape[0])]))
            else:
                if not no_oper:
                    CRPSS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i] / CRPS_oper[i] if CRPS_oper[i] > 0. else np.nan
                                                                for i in range(CRPSAn.shape[0])]))
                CRPSS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i] / CRPS_ol[i] if CRPS_ol[i] > 0. else np.nan
                                                          for i in range(CRPSAn.shape[0])]))

        else:
            if not no_decomp:
                if not no_oper:
                    CRPSS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_oper[i, 0] if (CRPS_oper[i, 0] > 0. and CRPS_oper[i, 0] >= CRPSAn[i, 0])
                                                                else CRPS_oper[i, 0] / CRPSAn[i, 0] - 1 if (CRPS_oper[i, 0] > 0. and CRPS_oper[i, 0] < CRPSAn[i, 0])
                                                                else np.nan
                                                                for i in range(CRPSAn.shape[0])]))
                    ReliS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 1] / CRPS_oper[i, 1] if (CRPS_oper[i, 1] > 0. and CRPS_oper[i, 1] >= CRPSAn[i, 1])
                                                                else CRPS_oper[i, 1] / CRPSAn[i, 1] - 1 if (CRPS_oper[i, 1] > 0. and CRPS_oper[i, 1] < CRPSAn[i, 1])
                                                                else np.nan
                                                                for i in range(CRPSAn.shape[0])]))
                CRPSS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_ol[i, 0] if (CRPS_ol[i, 0] > 0. and CRPS_ol[i, 0] >= CRPSAn[i, 0])
                                                          else CRPS_ol[i, 0] / CRPSAn[i, 0] - 1 if (CRPS_ol[i, 0] > 0. and CRPS_ol[i, 0] < CRPSAn[i, 0])
                                                          else np.nan
                                                          for i in range(CRPSAn.shape[0])]))
                ReliS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 1] / CRPS_ol[i, 1] if (CRPS_ol[i, 1] > 0. and CRPS_ol[i, 1] >= CRPSAn[i, 1])
                                                          else CRPS_ol[i, 1] / CRPSAn[i, 1] - 1 if (CRPS_ol[i, 1] > 0. and CRPS_ol[i, 1] < CRPSAn[i, 1])
                                                          else np.nan
                                                          for i in range(CRPSAn.shape[0])]))
            else:
                if not no_oper:
                    CRPSS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_oper[i, 0] if (CRPS_oper[i, 0] > 0. and CRPS_oper[i, 0] >= CRPSAn[i, 0])
                                                                else CRPS_oper[i, 0] / CRPSAn[i, 0] - 1 if (CRPS_oper[i, 0] > 0. and CRPS_oper[i, 0] < CRPSAn[i, 0])
                                                                else np.nan
                                                                for i in range(CRPSAn.shape[0])]))
                CRPSS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i] / CRPS_ol[i] if (CRPS_ol[i] > 0. and CRPS_ol[i] >= CRPSAn[i])
                                                          else CRPS_ol[i] / CRPSAn[i] - 1 if (CRPS_ol[i] > 0. and CRPS_ol[i] < CRPSAn[i])
                                                          else np.nan
                                                          for i in range(CRPSAn.shape[0])]))

    elif aggrTime and aggrSpace:
        ensAn_flat = ensAn.reshape(-1, ensAn.shape[-1])
        ensOl_flat = ensOl.reshape(-1, ensOl.shape[-1])
        oper_flat = oper.reshape(-1, oper.shape[-1])
        obs_flat = obs.flatten()
        if not no_decomp:
            if not no_oper:
                CRPS_oper = EnsembleScores(list(range(oper_flat.shape[0])),
                                           list(range(oper_flat.shape[0])),
                                           obs_flat,
                                           oper_flat.T,
                                           ).CRPS_decomp()
            CRPS_an = EnsembleScores(list(range(ensAn_flat.shape[0])),
                                     list(range(ensAn_flat.shape[0])),
                                     obs_flat,
                                     ensAn_flat.T,
                                     ).CRPS_decomp()
            CRPS_ol = EnsembleScores(list(range(ensOl_flat.shape[0])),
                                     list(range(ensOl_flat.shape[0])),
                                     obs_flat,
                                     ensOl_flat.T,
                                     ).CRPS_decomp()
            if not no_oper:
                CRPSS_oper = 1 - CRPS_an[0] / CRPS_oper[0]
                ReliS_oper = 1 - CRPS_an[1] / CRPS_oper[1]
            CRPSS_ol = 1 - CRPS_an[0] / CRPS_ol[0]
            ReliS_ol = 1 - CRPS_an[1] / CRPS_ol[1]

        else:
            if not no_oper:
                CRPS_oper = EnsembleScores(list(range(oper_flat.shape[0])),
                                           list(range(oper_flat.shape[0])),
                                           obs_flat,
                                           oper_flat.T,
                                           ).CRPS()
            CRPS_an = EnsembleScores(list(range(ensAn_flat.shape[0])),
                                     list(range(ensAn_flat.shape[0])),
                                     obs_flat,
                                     ensAn_flat.T,
                                     ).CRPS()
            CRPS_ol = EnsembleScores(list(range(ensOl_flat.shape[0])),
                                     list(range(ensOl_flat.shape[0])),
                                     obs_flat,
                                     ensOl_flat.T,
                                     ).CRPS()
            if not no_oper:
                CRPSS_oper = 1 - CRPS_an / CRPS_oper
            CRPSS_ol = 1 - CRPS_an / CRPS_ol
    if not no_decomp:
        if not no_oper:
            return CRPSS_ol, ReliS_ol, CRPSS_oper, ReliS_oper
        else:
            return CRPSS_ol, ReliS_ol
    else:
        if not no_oper:
            return CRPSS_ol, CRPSS_oper
        else:
            return CRPSS_ol


def RMSE_spread_bias_ol_oper(run, classesId=None, aggrSpace = False, aggrTime = True):
    if classesId is None:
        classesId = range(run.pgd.npts)
    obs, ensAn, ensOl, oper = frame_year_classes(run, classesId)

    RMSE_an = RMSE(ensAn, obs, aggrTime = aggrTime, aggrDomain = aggrSpace)
    RMSE_ol = RMSE(ensOl, obs, aggrTime = aggrTime, aggrDomain = aggrSpace)
    RMSE_oper = RMSE(oper, obs, aggrTime = aggrTime, aggrDomain = aggrSpace)
    bias_an = bias(ensAn, obs, aggrTime = aggrTime, aggrDomain = aggrSpace)
    bias_ol = bias(ensOl, obs, aggrTime = aggrTime, aggrDomain = aggrSpace)
    bias_oper = bias(oper, obs, aggrTime = aggrTime, aggrDomain = aggrSpace)
    spread_an = spread(ensAn, aggrTime = aggrTime, aggrDomain = aggrSpace)
    spread_ol = spread(ensOl, aggrTime = aggrTime, aggrDomain = aggrSpace)

    return RMSE_an, RMSE_ol, RMSE_oper, bias_an, bias_ol, bias_oper, spread_an, spread_ol


def compute_rankdiag(run, nbins=None):
    _, mobs, mens, _ = set_itimes_posts(run)
    if nbins is None:
        nbins = run.options.nmembers + 1
    if run.options.openloop == 'on':
        ens = run.ensProOl['DEP'][mens, :, :]
    else:
        ens = run.ensProAn['DEP'][mens, :, :]

    obs = np.ma.masked_invalid(run.obsTs['DEP'][mobs, :]).flatten()

    rd = rankDiagram(ens.reshape(ens.shape[0] * ens.shape[1],
                                 ens.shape[2]).T,
                     obs,
                     nbins = nbins)
    return rd


def compute_RMSE_spread_bias(run, no_oper = False, aggrTime = True, aggrDomain = True):
    """
    copute domain and time aggregated metrics.
    be careful: the spread is computed on the same dates as bias and RMSE for a proper assessment of the spread skill.
    """
    _, mobs, mens, moper = set_itimes_posts(run)
    if run.options.openloop == 'on':
        ens = run.ensProOl['DEP'][mens, :, :]
    else:
        ens = run.ensProAn['DEP'][mens, :, :]
    rmse = RMSE(ens, np.ma.masked_invalid(
        run.obsTs['DEP'][mobs, :]), aggrTime=aggrTime, aggrDomain=aggrDomain)
    sspread = spread(ens, truth = run.obsTs['DEP'][mobs, :], aggrTime=aggrTime, aggrDomain=aggrDomain)
    bbias = bias(ens, np.ma.masked_invalid(
        run.obsTs['DEP'][mobs, :]), aggrTime=aggrTime, aggrDomain=aggrDomain)
    if no_oper is False:
        rmse_oper = RMSE(run.oper['DEP'][moper, :], np.ma.masked_invalid(
            run.obsTs['DEP'][mobs, :]), aggrTime=aggrTime, aggrDomain=aggrDomain)
        bias_oper = bias(run.oper['DEP'][moper, :], np.ma.masked_invalid(
            run.obsTs['DEP'][mobs, :]), aggrTime=aggrTime, aggrDomain=aggrDomain)

        return rmse_oper, bias_oper, rmse, sspread, bbias
    else:
        return rmse, sspread, bbias


def set_itimes_posts(run):
    """
    set temporal masks to put in common an ensemble (an or ol), its obs timeseries and the associated oper run.
    """
    if run.options.openloop == 'on':
        timesEns = run.ensProOl['time']
    else:
        timesEns = run.ensProAn['time']
    year = timesEns[0].year
    mobs = [i for i, t in enumerate(run.obsTs['time'])
            if (t.hour == 6 and ((t.year == year and t.month > 9) or (t.month < 7 and t.year == year + 1)))]
    mens = [i for i, t in enumerate(timesEns) if (
        t.hour == 6 and (t.month > 9 or t.month < 7))]
    moper = [i for i, t in enumerate(run.oper['time']) if (
        t.hour == 6 and (t.month > 9 or t.month < 7))]

    return year, mobs, mens, moper
