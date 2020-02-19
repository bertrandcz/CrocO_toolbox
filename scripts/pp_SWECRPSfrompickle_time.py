# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:26:45 2019

@author: cluzetb
"""

import datetime
from snowtools_git.scores.generic import rankDiagram
import time

from matplotlib.dates import date2num

from CrocOpp import CrocOpp
from Ensemble import *
from Operators import *
from SemiDistributed import *
from crocO import set_options
import matplotlib.pyplot as plt
from plotcrocO import *
from utilcrocO import setSubsetclasses
from utilpp import RMSE, spread


start_time = time.time()
years = [2013, 2014, 2015, 2016]
dictmembers = {2013: [66, 12, 122, 149],
               2014: [69, 28, 2, 122],
               2015: [92, 97, 14, 141],
               2016: [50,
                      153,
                      90,
                      117
                      ]}

nens = 40
assimvars = 'DEP'
Neff = ['7']  # ['1', '10', '7']
runs = ['global', 'klocal1', 'rlocal']
# kind = 'fromol'
kind = 'postes'
score = 'CRPS'
aggrClasses = False
aggrTime = True
savefig = False
suffix = '' if nens == 160 else '_{0}'.format(nens) if 'DEP' not in assimvars else '_{0}_{1}'.format(nens, assimvars)
if len(Neff) > 0 and kind != 'postes':
    suffixs = [suffix + '_' + n  if n is not '7' else suffix for n in Neff]
    suffix = suffixs
    SUFFIX_NAME = '_40_DEP_neffstudy'
elif kind == 'postes':
    suffixs = [suffix + '_' + n for n in Neff]
    suffix = suffixs
    SUFFIX_NAME = ''

#SUFFIX_NAME = 'k2013_2014'

runssuf = [r + s for r in runs for s in suffix]
# aggregated scores.

if score == 'rankdiag':
    aggrClasses = True
    aggrTime = True
nbins = 41
dcolors = {'klocal1': 'y', 'klocal': 'c', 'klocal5': 'b', 'rlocal': 'g', 'rlocal_40': 'g', 'global': 'r', 'global_40': 'r', 'klocal5_40': 'b'}
RUN = dict()
for year in years:
    for mbsynth in dictmembers[year]:
        print('#########', year, mbsynth, '###########')
        for ik, key in enumerate(sorted(runssuf)):
            xp = '{0}_{1}_{2}'.format(year,
                                      mbsynth if (kind == 'fromol' and mbsynth != -1)
                                      else '{0}_{1}'.format(mbsynth, kind) if kind == 'postes'
                                      else kind if kind == 'baseline' else None,
                                      key)
            args = [
                '/home/cluzetb/snowtools_git/assim/crocO.py',
                '--xpid', xp,
                '--xpidol', 'art2_OL_{0}_t1500'.format(year),
                '-d', 'all',
                '--vars', assimvars,
                '--ppvars', 'SWE',
                '-o', 'compar2to3',
                '--classesE', '1800,2100,2400,2700,3000,3300,3600' if kind != 'postes' else '1200,1500,1800,2100,2400',
                '--classesS', '0,20,40',  # cheat with the 3 slopes for the pie plots
                '--classesA', 'SW,S,SE,E',
            ]
            options, conf = set_options(args)

            RUN[xp] = CrocOpp(options, conf)
            pgd = RUN[xp].pgd
            # set time and focus selection
            if kind != 'postes':
                focusCl = setSubsetclasses(pgd, options.classesE,
                                           options.classesA, ['0', '20'])[0]
                nfocusCl = [p for p in range(pgd.npts) if p not in focusCl ]
            else:
                focusCl = setSubsetclasses(pgd, options.classesE,
                                           options.classesA, ['0'])[0]
                nfocusCl = [p for p in range(pgd.npts) if p not in focusCl ]
            if aggrClasses:
                selFocus = focusCl
            else:
                selFocus = list(range(pgd.npts))

            itimes, itimesC = set_itimes(RUN[xp], clim = True, fromOl = True)
            fEnsProAn = dict()
            nfEnsProAn = dict()
            fAnCRPS = dict()
            nfAnCRPS = dict()

            if not os.path.exists(options.xpiddir + '/crocO/' + options.saverep):
                os.mkdir(options.xpiddir + '/crocO/' + options.saverep)
            os.chdir(options.xpiddir + '/crocO/' + options.saverep)
            # BC 07/01/20 dirty dirty
            os.chdir('/home/cluzetb/notebooks/articleGMD//pie/')
            if score == 'CRPS' or score == 'rankdiag':
                if ik == 0:
                    # remove beginning of the season and hours different from 12 and 29 feb if any

                    fTruth = dict()
                    fEnsProOl = dict()
                    fEnsProClim = dict()
                    nfTruth = dict()
                    nfEnsProOl = dict()
                    nfEnsProClim = dict()
                    fClimCRPS = dict()
                    fClimReli = dict()
                    fClimResol = dict()
                    fOlCRPS2 = dict()
                    fOlReli = dict()
                    fOlResol = dict()
                    nfClimCRPS = dict()
                    nfClimReli = dict()
                    nfClimResol = dict()
                    nfOlCRPS2 = dict()
                    nfOlReli = dict()
                    nfOlResol = dict()
                    fAnCRPS2 = dict()
                    fAnReli = dict()
                    fAnResol = dict()
                    nfAnCRPS2 = dict()
                    nfAnReli = dict()
                    nfAnResol = dict()
                    nfOlCRPS = dict()
                    frankOl = dict()
                    nfrankOl = dict()
                    frankAn = dict()
                    nfrankAn = dict()
                    OlCRPS = PrepAbs('all', options, ptinom='olcrps')
                    OlCRPS.data = dict()

                if type(options.ppvars) is str:
                    options.ppvars = [options.ppvars]

                for k in ['SWE']:
                    print('treated variable', k)

                    # OPENLOOP
                    if ik == 0:

                        fTruth[k] = RUN[xp].truth[k][:, selFocus][itimes, :]
                        fEnsProOl[k] = RUN[xp].ensProOl[k][:, selFocus, :][itimes, :, :]
                        if nens < fEnsProOl[k].shape[-1]:
                            fEnsProOl[k] = fEnsProOl[k][:, :, 0: nens]

                        if hasattr(RUN[xp], 'ensProClim'):
                            fEnsProClim[k] = RUN[xp].ensProClim[k][:, selFocus, :][itimesC, :, :]

                        # timeseries or scalars
                        if aggrClasses:

                            nfTruth[k] = RUN[xp].truth[k][:, nfocusCl][itimes, :]
                            nfEnsProOl[k] = RUN[xp].ensProOl[k][:, nfocusCl, :][itimes, :, :]
                            if nens < fEnsProOl[k].shape[-1]:
                                nfEnsProOl[k] = nfEnsProOl[k][:, :, 0: nens]
                            if hasattr(RUN[xp], 'ensProClim'):
                                nfEnsProClim[k] = RUN[xp].ensProClim[k][:, nfocusCl, :]
                                nfEnsProClim[k] = nfEnsProClim[k][itimesC, :, :]
                            # scalars
                            if aggrTime:

                                fTruth[k] = fTruth[k].flatten()
                                nfTruth[k] = nfTruth[k].flatten()

                                fEnsProOl[k] = fEnsProOl[k].reshape(-1, fEnsProOl[k].shape[-1])
                                nfEnsProOl[k] = nfEnsProOl[k].reshape(-1, nfEnsProOl[k].shape[-1])
                                if score == 'CRPS':

                                    fOlCRPS2[k], fOlReli[k], fOlResol[k] = EnsembleScores(list(range(fEnsProOl[k].shape[0])),
                                                                                          list(range(fEnsProOl[k].shape[0])),
                                                                                          fTruth[k],
                                                                                          fEnsProOl[k].T,
                                                                                          ).CRPS_decomp()
                                    nfOlCRPS[k] = EnsembleScores(list(range(nfEnsProOl[k].shape[0])),
                                                                 list(range(nfEnsProOl[k].shape[0])),
                                                                 nfTruth[k],
                                                                 nfEnsProOl[k].T,
                                                                 ).CRPS()
                                elif score == 'rankdiag':
                                    frankOl[k] = rankDiagram(fEnsProOl[k].T, fTruth[k], isSorted = False, nbins = nbins)
                                    nfrankOl[k] = rankDiagram(nfEnsProOl[k].T, nfTruth[k], isSorted = False, nbins = nbins)
                                if hasattr(RUN[xp], 'ensProClim'):
                                    fEnsProClim[k] = fEnsProClim[k].reshape(-1, fEnsProClim[k].shape[-1])
                                    nfEnsProClim[k] = nfEnsProClim[k].reshape(-1, nfEnsProClim[k].shape[-1])
                                    fClimCRPS[k], fClimReli[k], fClimResol[k] = EnsembleScores(list(range(fEnsProClim[k].shape[0])),
                                                                                               list(range(fEnsProClim[k].shape[0])),
                                                                                               fTruth[k],
                                                                                               fEnsProClim[k].T,
                                                                                               ).CRPS_decomp()
                                    nfClimCRPS[k], nfClimReli[k], nfClimResol[k] = EnsembleScores(list(range(nfEnsProClim[k].shape[0])),
                                                                                                  list(range(nfEnsProClim[k].shape[0])),
                                                                                                  fTruth[k],
                                                                                                  nfEnsProClim[k].T,
                                                                                                  ).CRPS_decomp()

                            # timeseries
                            else:
                                # OPENLOOP
                                lu = np.array([EnsembleScores(list(range(fEnsProOl[k].shape[1])),
                                                              list(range(fEnsProOl[k].shape[1])),
                                                              fTruth[k][it, :],
                                                              fEnsProOl[k][it, :, :].T,
                                                              ).CRPS_decomp() for it in range(len(itimes))])
                                fOlCRPS2[k] = np.array([l[0] for l in lu])
                                fOlReli[k] = np.array([l[1] for l in lu])
                                fOlResol[k] = np.array([l[2] for l in lu])
                                lu = np.array([EnsembleScores(list(range(nfEnsProOl[k].shape[1])),
                                                              list(range(nfEnsProOl[k].shape[1])),
                                                              nfTruth[k][it, :],
                                                              nfEnsProOl[k][it, :, :].T,
                                                              ).CRPS_decomp() for it in range(len(itimes))])
                                nfOlCRPS2[k] = np.array([l[0] for l in lu])
                                nfOlReli[k] = np.array([l[1] for l in lu])
                                nfOlResol[k] = np.array([l[2] for l in lu])

                                # CLIM
                                if hasattr(RUN[xp], 'ensProClim'):

                                    lu = np.array([EnsembleScores(list(range(fEnsProClim[k].shape[1])),
                                                                  list(range(fEnsProClim[k].shape[1])),
                                                                  fTruth[k][it, :],
                                                                  fEnsProClim[k][it, :, :].T,
                                                                  ).CRPS_decomp() for it in range(len(itimesC))])
                                    fClimCRPS[k] = np.array([l[0] for l in lu])
                                    fClimReli[k] = np.array([l[1] for l in lu])
                                    fClimResol[k] = np.array([l[2] for l in lu])
                                    lu = np.array([EnsembleScores(list(range(nfEnsProClim[k].shape[1])),
                                                                  list(range(nfEnsProClim[k].shape[1])),
                                                                  nfTruth[k][it, :],
                                                                  nfEnsProClim[k][it, :, :].T,
                                                                  ).CRPS_decomp() for it in range(len(itimesC))])
                                    nfClimCRPS[k] = np.array([l[0] for l in lu])
                                    nfClimReli[k] = np.array([l[1] for l in lu])
                                    nfClimResol[k] = np.array([l[2] for l in lu])
                        else:
                            OlCRPS.data[k] = np.array([EnsembleScores(list(range(fEnsProOl[k].shape[0])),
                                                                      list(range(fEnsProOl[k].shape[0])),
                                                                      fTruth[k][:, cl],
                                                                      fEnsProOl[k][:, cl, :].T,
                                                                      ).CRPS() for cl in range(pgd.npts)])

                    itimesAn = set_itimes(RUN[xp], fromOl = False)
                    fEnsProAn[k] = RUN[xp].ensProAn[k][:, selFocus, :]
                    fEnsProAn[k] = fEnsProAn[k][itimesAn, :, :]

                    if aggrClasses:
                        nfEnsProAn[k] = RUN[xp].ensProAn[k][:, nfocusCl, :]
                        nfEnsProAn[k] = nfEnsProAn[k][itimesAn, :, :]

                        if aggrTime:

                            fEnsProAn[k] = fEnsProAn[k].reshape(-1, fEnsProAn[k].shape[-1])
                            nfEnsProAn[k] = nfEnsProAn[k].reshape(-1, nfEnsProAn[k].shape[-1])
                            if score == 'CRPS':
                                fAnCRPS[k] = EnsembleScores(list(range(fEnsProAn[k].shape[0])),
                                                            list(range(fEnsProAn[k].shape[0])),
                                                            fTruth[k],
                                                            fEnsProAn[k].T,
                                                            ).CRPS_decomp()

                                nfAnCRPS[k] = EnsembleScores(list(range(nfEnsProAn[k].shape[0])),
                                                             list(range(nfEnsProAn[k].shape[0])),
                                                             nfTruth[k],
                                                             nfEnsProAn[k].T,
                                                             ).CRPS_decomp()
                            elif score == 'rankdiag':
                                frankAn[k] = rankDiagram(fEnsProAn[k].T, fTruth[k], isSorted = False, nbins = nbins)
                                nfrankAn[k] = rankDiagram(nfEnsProAn[k].T, nfTruth[k], isSorted = False, nbins = nbins)
                        else:
                            lu = np.array([EnsembleScores(list(range(fEnsProAn[k].shape[1])),
                                                          list(range(fEnsProAn[k].shape[1])),
                                                          fTruth[k][it, :],
                                                          fEnsProAn[k][it, :, :].T,
                                                          ).CRPS_decomp() for it in range(len(itimesAn))])
                            fAnCRPS2[k] = np.array([l[0] for l in lu])
                            fAnReli[k] = np.array([l[1] for l in lu])
                            fAnResol[k] = np.array([l[2] for l in lu])
                            lu = np.array([EnsembleScores(list(range(nfEnsProAn[k].shape[1])),
                                                          list(range(nfEnsProAn[k].shape[1])),
                                                          nfTruth[k][it, :],
                                                          nfEnsProAn[k][it, :, :].T,
                                                          ).CRPS_decomp() for it in range(len(itimesAn))])
                            nfAnCRPS2[k] = np.array([l[0] for l in lu])
                            nfAnReli[k] = np.array([l[1] for l in lu])
                            nfAnResol[k] = np.array([l[2] for l in lu])
                    else:
                        AnCRPS = PrepAbs('all', options, ptinom=xp)

                        AnCRPS.data = dict()
                        AnCRPS_skill = PrepAbs('all', options, ptinom=xp)
                        AnCRPS_skill.data = dict()
                        AnCRPS.data[k] = np.array([EnsembleScores(list(range(fEnsProOl[k].shape[0])),
                                                                  list(range(fEnsProOl[k].shape[0])),
                                                                  fTruth[k][:, cl],
                                                                  fEnsProAn[k][:, cl, :].T,
                                                                  ).CRPS() for cl in range(pgd.npts)])
                        AnCRPS_skill.data[k] = 1 - AnCRPS.data[k] / OlCRPS.data[k]
                if aggrClasses and not aggrTime:

                    if ik == 0:
                        fig, axes = plt.subplots(1, 2, figsize = (9, 3.5), sharey = True)
                        fig.subplots_adjust(top=0.93, left = 0.12, right = 0.99)

                        for ax in axes:
                            for tt in RUN[xp].conf.assimdates:
                                ax.axvline(x=tt, ls= ':', color = '0.7')
                        # observed classes OL
                        axes[0].plot(RUN[xp].ensProOl['time'][itimes], fOlCRPS2[k], color = 'k', label = 'openloop', linewidth = 1)
                        # axes[0].plot(RUN[xp].ensProOl['time'][itimes], fOlReli[k], color = 'k', label = 'reli', linewidth = 1)
                        # axes[0].plot(RUN[xp].ensProOl['time'][itimes], fOlResol[k], color = 'k', label = 'resol', linewidth = 1, ls = '--')

                        # CLIM
                        if hasattr(RUN[xp], 'ensProClim'):

                            axes[0].plot(RUN[xp].ensProOl['time'][itimes], fClimCRPS[k], color = 'brown', label = 'clim', linewidth = 1)
                            # axes[0].plot(RUN[xp].ensProOl['time'][itimes], fClimReli[k], color = 'brown', linewidth = 1)
                            # axes[0].plot(RUN[xp].ensProOl['time'][itimes], fClimResol[k], color = 'brown', linewidth = 1, ls = '--')

                        axes[0].set_title('Observed classes')
                        axes[0].set_ylabel(k + ' CRPS [$\mathrm{\mathsf{kgm^{-2}}}$]')
                        # non-observed classes
                        axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfOlCRPS2[k], color = 'k', linewidth = 1)
                        # axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfOlReli[k], color = 'k', linewidth = 1)
                        # axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfOlResol[k], color = 'k', linewidth = 1, ls = '--')

                        # CLIM
                        if hasattr(RUN[xp], 'ensProClim'):

                            axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfClimCRPS[k], color = 'brown', label = 'clim', linewidth = 1)
                            # axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfClimReli[k], color = 'brown', linewidth = 1)
                            # axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfClimResol[k], color = 'brown', linewidth = 1, ls = '--')
                        axes[1].set_title('Not observed classes')
                        fig.autofmt_xdate()

                    axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnCRPS2[k], color = dcolors[key], label = key, lw = 1)
                    # axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnReli[k], color = dcolors[key], lw=1)
                    # axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnResol[k], color = dcolors[key], ls = '--')
                    axes[0].legend()

                    axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnCRPS2[k], color = dcolors[key], label = key, lw = 1)
                    # axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnReli[k], color = dcolors[key], lw=1)
                    # axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnResol[k], color = dcolors[key], ls = '--')
                    axes[0].legend()

                elif not aggrClasses and aggrTime:
                    """
                    if ik == 0:
                        OlCRPS.isloaded = True
                        piefOl = Pie(OlCRPS,
                                     focusCl = focusCl,
                                     )
                        piefOl.plot(cmap = 'Reds', clims = [0, 100],
                                    title = 'openloop CRPS',
                                    savefig = True)

                    AnCRPS.isloaded = True
                    piefAn = Pie(AnCRPS,
                                 focusCl = focusCl,
                                 )
                    piefAn.plot(cmap = 'Reds', clims = [0, 100],
                                title = key + ' CRPS',
                                savefig = True)
                    """
                    AnCRPS_skill.isloaded = True
                    piefAn_sk = Pie(AnCRPS_skill,
                                    focusCl = focusCl, maketitle = True
                                    )
                    piefAn_sk.plot(cmap = 'RdBu', clims = [-1, 1],
                                   title = sorted(runs)[ik],
                                   colortitle = k + ' CRPSS', savefig = savefig)

                elif aggrClasses and aggrTime:
                    print('-------', xp, '-------')
                    if score == 'CRPS':
                        print(fAnCRPS)
                        print(nfAnCRPS)
                        print(fOlCRPS2)
                        print(nfOlCRPS)
                    elif score == 'rankdiag':
                        if not os.path.exists(RUN[xp].options.xpidoldir + '/crocO/pie/'):
                            os.mkdir(RUN[xp].options.xpidoldir + '/crocO/pie/')
                        fig, axes = plt.subplots(1, 2, figsize = (7, 3), sharey = True)
                        fig.subplots_adjust(wspace = 0.05)
                        axes[0].set_ylim([0, 1])
                        axes[0].set_xlim([-0.5, nbins - 0.5])
                        axes[1].set_xlim([-0.5, nbins - 0.5])
                        axes[0].bar(np.arange(len(frankOl[k][0])), frankOl[k][0], color = 'k', alpha = 0.5)
                        axes[0].bar(np.arange(len(frankAn[k][0])), frankAn[k][0], color = 'b', alpha = 0.5)
                        axes[1].bar(np.arange(len(nfrankOl[k][0])), nfrankOl[k][0], color = 'k', alpha = 0.5)
                        axes[1].bar(np.arange(len(nfrankAn[k][0])), nfrankAn[k][0], color = 'b', alpha = 0.5)
                        axes[0].set_ylabel('frequency', fontsize = 14)
                        axes[0].text(0.5, -0.15, 'observed', va  = 'center', ha = 'center', fontsize = 14, transform = axes[0].transAxes)
                        axes[1].text(0.5, -0.15, 'not observed', va  = 'center', ha = 'center', fontsize = 14, transform = axes[1].transAxes)
                        plt.show()

                        # fig.savefig('../pie/{0}_{1}_{2}.png'.format(score, year, mbsynth))
                '''BC 10/12/19 smells like it's deprecated
                elif score == 'RMSE':
                    """
                    timseseries of the ensemble mean RMSE aggregated over the domain points.
                    """
                    if ik == 0:
                        fTruth = dict()
                        nfTruth = dict()
                        fEnsProOl = dict()
                        nfEnsProOl = dict()
                        fEnsProAn = dict()
                        nfEnsProAn = dict()
                        fOldisp = dict()
                        nfOldisp = dict()
                        fOlRMSE = dict()
                        nfOlRMSE = dict()
                        fAnRMSE = dict()
                        fAndisp = dict()
                        nfAnRMSE = dict()
                        nfAndisp = dict()
                    for k in['SWE']:
                        if ik == 0:
                            fTruth[k] = RUN[xp].ensProOl[k][:, selFocus, RUN[xp].mbsynth][itimes, :]
                            nfTruth[k] = RUN[xp].ensProOl[k][:, nfocusCl, RUN[xp].mbsynth][itimes, :]
                            fEnsProOl[k] = RUN[xp].ensProOl[k][:, selFocus, :][itimes, :, :]
                            nfEnsProOl[k] = RUN[xp].ensProOl[k][:, nfocusCl, :][itimes, :, :]
                            # compute RMSE and spread.
                            # roll axis brings nth dimension to the front to iterate on it.
    
                            fOlRMSE[k] = RMSE(fEnsProOl[k], fTruth[k])
                            fOldisp[k] = spread(fEnsProOl[k])
                            nfOlRMSE[k] = RMSE(nfEnsProOl[k], nfTruth[k])
                            nfOldisp[k] = spread(nfEnsProOl[k])
                        if ik == 0:
                            fig, axes = plt.subplots(1, 2, figsize = (9, 3.5))
                            fig.subplots_adjust(top=0.93, left = 0.09, right = 0.99)
                            for ax in axes:
                                for tt in RUN[xp].conf.assimdates:
                                    ax.axvline(x=tt, ls= ':', color = '0.7')
                            axes[0].plot(RUN[xp].ensProOl['time'][itimes], fOlRMSE[k], color = 'k', label = '$\mathrm{\mathsf{RMSE_{OL}}}$', linewidth = 1, zorder = 2)
                            axes[0].plot(RUN[xp].ensProOl['time'][itimes], fOldisp[k], color = 'k', label = '$\sigma_{OL}$', linewidth = 1, ls = '--', zorder = 1)
                            axes[1].plot(RUN[xp].ensProOl['time'][itimes], fOlRMSE[k], color = 'k', linewidth = 1, zorder = 2)
                            axes[1].plot(RUN[xp].ensProOl['time'][itimes], nfOldisp[k], color = 'k', linewidth = 1, ls = '--', zorder = 1)
                            fig.autofmt_xdate()
                            axes[0].set_ylabel('SWE [$\mathrm{\mathsf{kgm^{-2}}}$]')
                            axes[0].set_title('Observed classes')
                            axes[1].set_title(' Not Observed classes')
                        itimesAn = set_itimes(RUN[xp], fromOl = False)
                        fEnsProAn[k] = RUN[xp].ensProAn[k][:, selFocus, :][itimesAn, :, :]
                        nfEnsProAn[k] = RUN[xp].ensProAn[k][:, nfocusCl, :][itimesAn, :, :]
                        fAnRMSE[k] = RMSE(fEnsProAn[k], fTruth[k])
                        fAndisp[k] = spread(fEnsProAn[k])
                        nfAnRMSE[k] = RMSE(nfEnsProAn[k], nfTruth[k])
                        nfAndisp[k] = spread(nfEnsProAn[k])
    
                        axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnRMSE[k], color = dcolors[key], label = key, linewidth = 1, zorder = 2)
                        axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAndisp[k], color = dcolors[key], linewidth = 1, ls = '--', zorder = 1)
                        axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnRMSE[k], color = dcolors[key], linewidth = 1, zorder = 2)
                        axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAndisp[k], color = dcolors[key], linewidth = 1, ls = '--', zorder = 1)
                    axes[0].legend()
                '''

        if aggrClasses and not aggrTime:
            if not os.path.exists(RUN[xp].options.xpidoldir + '/crocO/pie/'):
                os.mkdir(RUN[xp].options.xpidoldir + '/crocO/pie/')
            axes[0].set_ylim([0, max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])])
            axes[1].set_ylim([0, max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])])
            plt.show()
            fig.savefig(RUN[xp].options.xpidoldir + '/crocO/pie/{0}_timeseries_{1}_{2}.png'.format(score, year, mbsynth))

el_time = time.time() - start_time
print(el_time)
