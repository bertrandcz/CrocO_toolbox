'''
Created on 15 oct. 2019

@author: cluzetb

full eval on SWE, aggrclasses + time
'''
from CrocOpp import CrocOpp
from crocO import set_options
from scores.ensemble import EnsembleScores
import sys
import time
from utilcrocO import setSubsetclasses
from utilpp import set_itimes, RMSE, spread

import numpy as np
import pandas as pd


np.set_printoptions(threshold=sys.maxsize)
start_time = time.time()
# params
# ##############################################

# read the args
aggrClasses = True
aggrTime = True

years = [2013, 2014, 2015, 2016]
dictmembers = {2013: [66, 12, 122, 149],
               2014: [69, 28, 2, 122],
               2015: [92, 97, 14, 141],
               2016: [50, 153, 90, 117]}

nens = 40
# assimvars = 'B4,B5,SCF'
assimvars = 'DEP'
# assimvars = 'B4,B5,DEP,SWE'
ppvars = 'DEP,SWE'
k = 'SWE'
# Neff = ['1', '10', '7']
Neff = ['7']
kind = 'fromol'
scores = ['CRPS']
# ###############################################
suffix = '' if nens == 160 else '_{0}'.format(nens) if 'DEP' not in assimvars else '_{0}_{1}'.format(nens, assimvars)
if len(Neff) > 0 and kind != 'postes':
    suffixs = [suffix + '_' + n  if n is not '7' else suffix for n in Neff]
    suffix = suffixs
    SUFFIX_NAME = '_40_DEP_neffstudy'
elif kind == 'postes':
    suffixs = [suffix + '_' + n for n in Neff]
    suffix = suffixs
    SUFFIX_NAME = ''
# SUFFIX_NAME = 'k2013_2014'
runs = ['global', 'klocal5', 'rlocal']
# runs = ['klocal5', 'klocal1', 'klocal10']

runs = [r + s for r in runs for s in suffix]
pdruns = runs + [ 'ol', ]  # 'cl']

if kind == 'fromol' or kind == 'postes':
    dictSynth = dictmembers

elif kind == 'baseline':
    dictSynth = {y: [-1] for y in years}
# initialize the score
index = pd.MultiIndex.from_tuples([(year, mbsynth, f)
                                   for year in years for mbsynth in dictSynth[year] for f in ['f', 'nf']])
if 'CRPS' in scores:
    score_type = ['CRPS', 'Reli', 'Resol']
    dfScores_CRPS = pd.DataFrame(np.empty((len(index), len(pdruns) * len(score_type) )), columns = ['_'.join([r, s]) for r in pdruns for s in score_type], index = index)
if 'RMSE' in scores:
    score_type = ['RMSE', 'spread', 'ss']
    dfScores_RMSE = pd.DataFrame(np.empty((len(index), len(pdruns) * len(score_type)  )), columns = ['_'.join([r, s]) for r in pdruns for s in score_type], index = index)

# load OLs
OL = dict()
fTruth = dict()
nfTruth = dict()
itimes = dict()
itimesC = dict()

for year in years:
    xp = 'art2_OL_{0}_t1500'.format(year)
    args = [
        '/home/cluzetb/assim/crocO.py',
        '--xpid', xp,
        '-d', 'all',
        '--vars', assimvars,
        '--ppvars', ppvars,
        '-o', 'commitGMDtest',
        '--classesE', '1800,2100,2400,2700,3000,3300,3600',
        '--classesS', '0,20',
        '--classesA', 'SW,S,SE,E',
        # '--clim'
    ]
    options, conf = set_options(args)

    OL[year] = CrocOpp(options, conf)  # full loading, expensive but necessary in exploration mode.
    pgd = OL[year].pgd
    # set time and focus selection
    focusCl = setSubsetclasses(pgd, options.classesE,
                               options.classesA, options.classesS)[0]
    nfocusCl = [p for p in range(pgd.npts) if p not in focusCl ]

    itimes[year], itimesC[year] = set_itimes(OL[year], clim = True)

    fEnsProOl = OL[year].ensProOl[k][:, focusCl, :][itimes[year], :, :]
    # fEnsProCl = OL[year].ensProClim[k][:, focusCl, :][itimesC[year], :, :]
    nfEnsProOl = OL[year].ensProOl[k][:, nfocusCl, :][itimes[year], :, :]
    # nfEnsProCl = OL[year].ensProClim[k][:, nfocusCl, :][itimesC[year], :, :]
    if nens < fEnsProOl.shape[-1]:
        fEnsProOl = fEnsProOl[:, :, 0: nens]
        nfEnsProOl = nfEnsProOl[:, :, 0: nens]

    fTruth[year] = dict()
    nfTruth[year] = dict()
    for mbsynth in dictSynth[year]:
        # tell read_truth where to read the truth
        OL[year].mbsynth = mbsynth - 1
        OL[year].readTruth()
        fTruth[year][mbsynth] = OL[year].truth[k][:, focusCl][itimes[year], :]
        nfTruth[year][mbsynth] = OL[year].truth[k][:, nfocusCl][itimes[year], :]
        print('shptruth', np.shape(fTruth[year][mbsynth]), np.shape(nfTruth[year][mbsynth]))
    for mbsynth in dictSynth[year]:
        if 'CRPS' in scores:
            fEnsProOl_flat = fEnsProOl.reshape(-1, fEnsProOl.shape[-1])
            # fEnsProCl_flat = fEnsProCl.reshape(-1, fEnsProCl.shape[-1])
            nfEnsProOl_flat = nfEnsProOl.reshape(-1, nfEnsProOl.shape[-1])
            # nfEnsProCl_flat = nfEnsProCl.reshape(-1, nfEnsProCl.shape[-1])
            dfScores_CRPS.loc[(year, mbsynth, 'f')]['ol' + '_CRPS', 'ol' + '_Reli', 'ol' + '_Resol'] = EnsembleScores(list(range(fEnsProOl_flat.shape[0])),
                                                                                                                      list(range(fEnsProOl_flat.shape[0])),
                                                                                                                      fTruth[year][mbsynth].flatten(),
                                                                                                                      fEnsProOl_flat.T,
                                                                                                                      ).CRPS_decomp()
            dfScores_CRPS.loc[(year, mbsynth, 'nf')]['ol' + '_CRPS', 'ol' + '_Reli', 'ol' + '_Resol'] = EnsembleScores(list(range(nfEnsProOl_flat.shape[0])),
                                                                                                                       list(range(nfEnsProOl_flat.shape[0])),
                                                                                                                       nfTruth[year][mbsynth].flatten(),
                                                                                                                       nfEnsProOl_flat.T,
                                                                                                                       ).CRPS_decomp()
            # dfScores_CRPS.loc[(year, mbsynth, 'f')]['cl' + '_CRPS', 'cl' + '_Reli', 'cl' + '_Resol'] = EnsembleScores(list(range(fEnsProCl_flat.shape[0])),
            #                                                                                                            list(range(fEnsProCl_flat.shape[0])),
            #                                                                                                           fTruth[year][mbsynth].flatten(),
            #                                                                                                           fEnsProCl_flat.T,
            #                                                                                                           ).CRPS_decomp()
            # dfScores_CRPS.loc[(year, mbsynth, 'nf')]['cl' + '_CRPS', 'cl' + '_Reli', 'cl' + '_Resol'] = EnsembleScores(list(range(nfEnsProCl_flat.shape[0])),
            #                                                                                                            list(range(nfEnsProCl_flat.shape[0])),
            #                                                                                                            nfTruth[year][mbsynth].flatten(),
            #                                                                                                            nfEnsProCl_flat.T,
            #                                                                                                            ).CRPS_decomp()
        if 'RMSE' in scores:
            dfScores_RMSE.loc[(year, mbsynth, 'f')]['ol' + '_RMSE', 'ol' + '_spread']  = RMSE(fEnsProOl, fTruth[year][mbsynth], aggrTime = aggrTime), spread(fEnsProOl, aggrTime = aggrTime)
            dfScores_RMSE.loc[(year, mbsynth, 'nf')]['ol' + '_RMSE', 'ol' + '_spread'] = RMSE(nfEnsProOl, nfTruth[year][mbsynth], aggrTime = aggrTime), spread(nfEnsProOl, aggrTime = aggrTime)
            # dfScores_RMSE.loc[(year, mbsynth, 'f')]['cl' + '_RMSE', 'cl' + '_spread']  = RMSE(fEnsProCl, fTruth[year][mbsynth], aggrTime = aggrTime), spread(fEnsProCl, aggrTime = aggrTime)
            # dfScores_RMSE.loc[(year, mbsynth, 'nf')]['cl' + '_RMSE', 'cl' + '_spread'] = RMSE(nfEnsProCl, nfTruth[year][mbsynth], aggrTime = aggrTime), spread(nfEnsProCl, aggrTime = aggrTime)

itimesAn = dict()
for year in years:
    for mbsynth in dictSynth[year]:
        for rr in runs:
            xp = '{0}_{1}_{2}'.format(year,
                                      mbsynth if (kind == 'fromol' and mbsynth != -1)
                                      else '{0}_{1}'.format(mbsynth, kind) if kind == 'postes'
                                      else kind if kind == 'baseline' else None,
                                      rr)
            print('xp', xp)
            args = [
                '/home/cluzetb/snowtools_git/assim/crocO.py',
                '--xpid', xp,
                '-d', 'all',
                '--vars', 'B4,B5,SCF' if kind != 'postes' else 'SCF',
                '--ppvars', 'B4,B5,DEP,SWE'if kind != 'postes' else 'DEP,SWE',
                '-o', 'commitGMDtest',
                '--classesE', '1800,2100,2400,2700,3000,3300,3600'if kind != 'postes' else '1200,1500,1800,2100,2400',
                '--classesS', '0,20' if kind != 'postes' else '0',
                '--classesA', 'SW,S,SE,E',
            ]
            options, conf = set_options(args)
            # try:
            run = CrocOpp(options, conf)
            # except Exception:
            # print('bug with xp {0}, maybe it doesnot exist'.format(xp))
            # continue
            itimesAn[year] = set_itimes(run)
            # check
            # if baseline is False and mbsynth - 1 != run.mbsynth:
            #    raise Exception('mbsynth confusion !!!')
            # else:
            fEnsProAn = run.ensProAn[k][:, focusCl, :][itimesAn[year], :, :]
            nfEnsProAn = run.ensProAn[k][:, nfocusCl, :][itimesAn[year], :, :]
            print('shpan', fEnsProAn.shape, nfEnsProAn.shape)
            if 'CRPS' in scores:
                fEnsProAn_flat = fEnsProAn.reshape(-1, run.ensProAn[k].shape[-1])
                nfEnsProAn_flat = nfEnsProAn.reshape(-1, run.ensProAn[k].shape[-1])

                dfScores_CRPS.loc[(year, mbsynth, 'f')][rr + '_CRPS', rr + '_Reli', rr + '_Resol'] = EnsembleScores(list(range(fEnsProAn_flat.shape[0])),
                                                                                                                    list(range(fEnsProAn_flat.shape[0])),
                                                                                                                    fTruth[year][mbsynth].flatten(),
                                                                                                                    fEnsProAn_flat.T,
                                                                                                                    ).CRPS_decomp()
                dfScores_CRPS.loc[(year, mbsynth, 'nf')][rr + '_CRPS', rr + '_Reli', rr + '_Resol'] = EnsembleScores(list(range(nfEnsProAn_flat.shape[0])),
                                                                                                                     list(range(nfEnsProAn_flat.shape[0])),
                                                                                                                     nfTruth[year][mbsynth].flatten(),
                                                                                                                     nfEnsProAn_flat.T,
                                                                                                                     ).CRPS_decomp()

            if 'RMSE' in scores:
                dfScores_RMSE.loc[(year, mbsynth, 'f')][rr + '_RMSE', rr + '_spread'] = RMSE(fEnsProAn, fTruth[year][mbsynth], aggrTime = aggrTime), spread(fEnsProAn, aggrTime = aggrTime)
                dfScores_RMSE.loc[(year, mbsynth, 'nf')][rr + '_RMSE', rr + '_spread'] = RMSE(nfEnsProAn, nfTruth[year][mbsynth], aggrTime = aggrTime), spread(nfEnsProAn, aggrTime = aggrTime)
for rr in pdruns:
    if 'RMSE' in scores:
        dfScores_RMSE[rr + '_ss'] = dfScores_RMSE[rr + '_spread'].div(dfScores_RMSE[rr + '_RMSE'])
# save the pickles
if 'CRPS' in scores:
    dfScores_CRPS.to_pickle('/home/cluzetb/notebooks/dfScores_{0}_{1}{2}.pkl'.format('CRPS', 'quant' if kind == 'fromol'
                                                                                     else 'quant_{0}'.format(kind) if kind == 'postes'
                                                                                     else kind if kind == 'baseline'else None, SUFFIX_NAME))
if 'RMSE' in scores:
    dfScores_RMSE.to_pickle('/home/cluzetb/notebooks/dfScores_{0}_{1}{2}.pkl'.format('RMSE', 'quant' if kind == 'fromol'
                                                                                     else 'quant_{0}'.format(kind) if kind == 'postes'
                                                                                     else kind if kind == 'baseline'else None, SUFFIX_NAME))
elapsed_time = time.time() - start_time
print('elapsed_time ', elapsed_time)
