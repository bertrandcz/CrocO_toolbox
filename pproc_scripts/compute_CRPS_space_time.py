# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:26:45 2019

@author: cluzetb
"""

from CramponPp import CramponPp
from SemiDistributed import PrepAbs
from consts import CRAMPON
from crampon import set_options
import os
from plotcrampon import Pie
from scores.ensemble import EnsembleScores
import time
from utilcrampon import setSubsetclasses
from utilplot import set_colors, set_title
from utilpp import set_itimes

import matplotlib.pyplot as plt
import numpy as np


def compute_CRPS_space_time(nens, dictmembers, runs, assimvars, saverep = 'gmd', kind='fromol', aggrClasses = False, aggrTime = True, savefig=False, savepath = None):
    '''
    Compute SWE CRPS over a set of experiments.
    CRPS can be aggregated over time or space (discriminating between observed and non observed classes.)
    for aggregation over both time and space, use multipp.py
    saverep : the crampon subdir where to look for pickle files
    '''
    start_time = time.time()

    if aggrClasses and aggrTime:
        raise Exception("for aggregation over time and space, use multipp.py")
    if not aggrClasses and not aggrTime:
        raise Exception("please specify one dimension to aggregate along")

    RUN = dict()
    k = 'SWE'
    for year in dictmembers:
        for mbsynth in dictmembers[year]:
            print('#########', year, mbsynth, '###########')
            for ik, key in enumerate(runs):
                xp = '{0}_{1}_{2}'.format(year,
                                          mbsynth if (kind == 'fromol' and mbsynth != -1)
                                          else '{0}_{1}'.format(mbsynth, kind) if kind == 'postes'
                                          else kind if kind == 'baseline' else None,
                                          key)
                args = [
                    CRAMPON + '/crampon.py',
                    '--xpid', xp,
                    '--xpidol', 'art2_OL_{0}_t1500'.format(year),
                    '-d', 'all',
                    '--vars', assimvars,  # useful if the pickles do not exist
                    '--ppvars', 'SWE',
                    '-o', saverep,
                    '--nmembers', str(nens),
                    '--classes_e', '1800,2100,2400,2700,3000,3300,3600' if kind != 'postes' else '1200,1500,1800,2100,2400',
                    '--classes_s', '0,20,40',  # cheat with the 3 slopes for the pie plots
                    '--classes_a', 'SW,S,SE,E',
                    '--readtruth'
                ]
                options = set_options(args)

                RUN[xp] = CramponPp(options)
                pgd = RUN[xp].pgd
                # set time and focus selection
                if kind != 'postes':
                    focusCl = setSubsetclasses(pgd, options.classes_e,
                                               options.classes_a, ['0', '20'])[0]
                    nfocusCl = [p for p in range(pgd.npts) if p not in focusCl ]
                else:
                    focusCl = setSubsetclasses(pgd, options.classes_e,
                                               options.classes_a, ['0'])[0]
                    nfocusCl = [p for p in range(pgd.npts) if p not in focusCl ]
                if aggrClasses:
                    selFocus = focusCl
                elif aggrTime:
                    selFocus = list(range(pgd.npts))

                itimes, itimesC = set_itimes(RUN[xp], clim = True, fromOl = True)
                fEnsProAn = dict()
                nfEnsProAn = dict()

                if not os.path.exists(options.xpiddir + '/crampon/' + options.saverep):
                    os.mkdir(options.xpiddir + '/crampon/' + options.saverep)
                os.chdir(options.xpiddir + '/crampon/' + options.saverep)
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
                    OlCRPS = PrepAbs('all', options, ptinom='olcrps')
                    OlCRPS.data = dict()

                    fTruth[k] = RUN[xp].truth[k][:, selFocus][itimes, :]
                    fEnsProOl[k] = RUN[xp].ensProOl[k][:, selFocus, :][itimes, :, :]
                    if nens < fEnsProOl[k].shape[-1]:
                        fEnsProOl[k] = fEnsProOl[k][:, :, 0: nens]

                    if hasattr(RUN[xp], 'ensProClim'):
                        fEnsProClim[k] = RUN[xp].ensProClim[k][:, selFocus, :][itimesC, :, :]

                    # timeseries
                    if aggrClasses:

                        nfTruth[k] = RUN[xp].truth[k][:, nfocusCl][itimes, :]
                        nfEnsProOl[k] = RUN[xp].ensProOl[k][:, nfocusCl, :][itimes, :, :]
                        if nens < fEnsProOl[k].shape[-1]:
                            nfEnsProOl[k] = nfEnsProOl[k][:, :, 0: nens]
                        if hasattr(RUN[xp], 'ensProClim'):
                            nfEnsProClim[k] = RUN[xp].ensProClim[k][:, nfocusCl, :]
                            nfEnsProClim[k] = nfEnsProClim[k][itimesC, :, :]
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

                    elif aggrTime:
                        OlCRPS.data[k] = np.array([EnsembleScores(list(range(fEnsProOl[k].shape[0])),
                                                                  list(range(fEnsProOl[k].shape[0])),
                                                                  fTruth[k][:, cl],
                                                                  fEnsProOl[k][:, cl, :].T,
                                                                  ).CRPS() for cl in range(pgd.npts)])

                # every ik
                itimesAn = set_itimes(RUN[xp], fromOl = False)
                fEnsProAn[k] = RUN[xp].ensProAn[k][:, selFocus, :]
                fEnsProAn[k] = fEnsProAn[k][itimesAn, :, :]

                if aggrClasses:
                    nfEnsProAn[k] = RUN[xp].ensProAn[k][:, nfocusCl, :]
                    nfEnsProAn[k] = nfEnsProAn[k][itimesAn, :, :]

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
                elif aggrTime:
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

            # plotting
                if aggrClasses:

                    if ik == 0:
                        fig, axes = plt.subplots(1, 2, figsize = (9, 3.5), sharey = True)
                        fig.subplots_adjust(top=0.93, left = 0.12, right = 0.99)

                        for ax in axes:
                            for tt in RUN[xp].assimdates:
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

                    axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnCRPS2[k], color = set_colors(key), label = key, lw = 1)
                    # axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnReli[k], color = set_colors(key), lw=1)
                    # axes[0].plot(RUN[xp].ensProAn['time'][itimesAn], fAnResol[k], color = set_colors(key), ls = '--')
                    axes[0].legend()

                    axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnCRPS2[k], color = set_colors(key), label = key, lw = 1)
                    # axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnReli[k], color = set_colors(key), lw=1)
                    # axes[1].plot(RUN[xp].ensProAn['time'][itimesAn], nfAnResol[k], color = set_colors(key), ls = '--')
                    axes[0].legend()

                elif aggrTime:
                    print('...plotting...', xp)
                    AnCRPS_skill.isloaded = True
                    piefAn_sk = Pie(AnCRPS_skill,
                                    focusCl = focusCl, maketitle = True
                                    )
                    piefAn_sk.plot(cmap = 'RdBu', clims = [-1, 1],
                                   title = set_title(runs[ik]),
                                   colortitle = k + ' CRPSS', savefig = savefig,
                                   savepath = (savepath[0:-4] + '_' + set_title(runs[ik]) + savepath[-4:]) if savepath else savepath)

    if aggrClasses:
        if not os.path.exists(RUN[xp].options.xpidoldir + '/crampon/pie/'):
            os.mkdir(RUN[xp].options.xpidoldir + '/crampon/pie/')
        axes[0].set_ylim([0, max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])])
        axes[1].set_ylim([0, max(axes[0].get_ylim()[1], axes[1].get_ylim()[1])])
        plt.show()
        fig.savefig(RUN[xp].options.xpidoldir + '/crampon/pie/CRPS_timeseries_{0}_{1}.png'.format(year, mbsynth))
    el_time = time.time() - start_time
    print(el_time)
