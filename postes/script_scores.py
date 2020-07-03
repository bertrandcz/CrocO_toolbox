'''
Created on 13 mai 2020

@author: cluzetb
'''

from scores.ensemble import EnsembleScores

import numpy as np

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


def CRPSS_ol_oper(run, year, classesId):
    mobs = [i for i, t in enumerate(run.obsTs['time'])
            if (t.hour == 6 and ((t.year == year and t.month > 9) or (t.month < 7 and t.year == year + 1)))]
    mensAn = [i for i, t in enumerate(run.ensProAn['time']) if (t.hour == 6 and (t.month > 9 or t.month < 7))]
    mensOl = [i for i, t in enumerate(run.ensProOl['time']) if (t.hour == 6 and (t.month > 9 or t.month < 7))]
    moper = [i for i, t in enumerate(run.oper['time']) if (t.hour == 6 and (t.month > 9 or t.month < 7))]
    obs = np.ma.masked_invalid(run.obsTs['DEP'][mobs, :][:, classesId])
    ensAn = run.ensProAn['DEP'][mensAn, :, :][:, classesId, :]
    ensOl = run.ensProOl['DEP'][mensOl, :, :][:, classesId, :]
    oper = run.oper['DEP'][moper, :][:, classesId]
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
    CRPS_oper = np.array([EnsembleScores(list(range(oper.shape[0])),
                                         list(range(oper.shape[0])),
                                         obs[:, cl],
                                         oper[:, cl, :].T,
                                         ).CRPS_decomp() for cl in range(len(classesId))])
    CRPSS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_oper[i, 0] if CRPS_oper[i, 0] > 0. else np.nan
                                                for i in range(CRPSAn.shape[0])]))
    ReliS_oper = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 1] / CRPS_oper[i, 1] if CRPS_oper[i, 0] > 0. else np.nan
                                                for i in range(CRPSAn.shape[0])]))
    CRPSS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 0] / CRPS_ol[i, 0] if CRPS_ol[i, 0] > 0. else np.nan
                                              for i in range(CRPSAn.shape[0])]))
    ReliS_ol = np.ma.masked_invalid(np.array([1 - CRPSAn[i, 1] / CRPS_ol[i, 1] if CRPS_ol[i, 0] > 0. else np.nan
                                              for i in range(CRPSAn.shape[0])]))

    return CRPSS_ol, ReliS_ol, CRPSS_oper, ReliS_oper
