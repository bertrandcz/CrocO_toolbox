#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 11 févr. 2019

@author: cluzetb
'''
import os

from Ensemble import PrepEnsBg, PrepEnsAn
from Operators import PrepEnsOperator
from SemiDistributed import PrepBg, PrepAbs, Synthetic, Real
import matplotlib.pyplot as plt
import numpy as np
from plotcrampon import Pie
from utilcrampon import setSubsetclasses


class PostCroco(object):
    '''
    post processing class
    '''

    def __init__(self, xpiddir, xpidobsdir, options, **kwargs):
        '''

        '''
        self.options = options
        self.xpiddir = xpiddir
        self.xpidobsdir = xpidobsdir

        for key, value in list(kwargs.items()):  # user can add any kwargument to this object
            setattr(self, key, value)

    def run(self, directFromXp = True):
        if self.options.todo == 'pievar':
            for dd in self.options.dates:
                os.chdir(dd)
                # for mb in range(1, self.options.nmembers + 1)
                for mb in [5]:
                    p = PrepBg(dd, mb, self.options)
                    pie = Pie(p)
                    pie.plot()
                os.chdir('..')
        # ((CORRELATION PLOT (only works with 1 class)))
        if self.options.todo == 'corr':
            for dd in self.options.dates:
                os.chdir(dd)
                self.pb = PrepEnsBg(self.options, dd)
                self.pa = PrepEnsAn(self.options, dd, directFromXp = False)
                # pa = PrepEnsAn(self.options, dd)
                self.obs = Synthetic(self.xpiddir, dd, self.options)  # load pre-existing Synthetic, don't generate it.
                op = PrepEnsOperator(self.pb, sdObj=self.obs)
                # corr = op.m21diff(reverse=True)
                subset, _ = setSubsetclasses(self.pb.ens[1].pgd, self.options.classesE, self.options.classesA, self.options.classesS)
                for cl in subset:
                    print(('printing in', cl))
                    corSd = op.bgcovariance(cl)
                    piecor = Pie(corSd, focusCl=cl,)
                    piecor.plot( cmap = 'RdBu', savefig=True, clims = [-1., 1.])

                    # plt.savefig('../pie/corr_' + str(int(pb.ens[1].pgd.elevClass[cl])) + '_' + str(int(pb.ens[1].pgd.aspectClass[cl])) +
                    #            '_' + str(int(np.arctan(pb.ens[1].pgd.slopeClass[cl]) * 180. / np.pi)) + '_' + dd + '.png')
                os.chdir('..')
        if self.options.todo == 'pfpython':
            if not hasattr(self, 'pf'):
                raise Exception
            for dd in self.options.dates:
                os.chdir(dd)
                self.pf.pb.calcmedian(ptinombase = 'bg')
                self.pf.pag.calcmedian(ptinombase = 'ang')
                self.pf.pal.calcmedian(ptinombase = 'anl')
                clims = [0., 1.]
                piemedb = Pie(self.pf.pb.median)
                piemedb.plot(savefig = True, clims = clims)

                piemedag = Pie(self.pf.pag.median)
                piemedag.plot(savefig = True, clims = clims)
                piemedal = Pie(self.pf.pal.median)
                piemedal.plot(savefig = True, clims = clims)
                # pieplot obs
                pieobs = Pie(self.pf.obs)
                pieobs.plot(savefig = True, clims = clims)
                # for i, a in enumerate(ax[:, 0]):
                #    a.set_ylabel([self.pf.pb.median.data.keys()[i]])

                # plt.savefig('pie/pfanalysis.png')
                innov = PrepAbs(dd, self.options, ptinom='innov')
                incrg = PrepAbs(dd, self.options, ptinom='incrg')
                incrl = PrepAbs(dd, self.options, ptinom='incrl')
                resg = PrepAbs(dd, self.options, ptinom='resg')
                resl = PrepAbs(dd, self.options, ptinom='resl')
                innov.data = dict()
                incrg.data = dict()
                resg.data = dict()
                incrl.data = dict()
                resl.data = dict()

                for k in list(self.pf.pb.median.data.keys()):
                    innov.data[k] = self.pf.obs.data[k] - self.pf.pb.median.data[k]

                    incrg.data[k] = self.pf.pag.median.data[k] - self.pf.pb.median.data[k]
                    resg.data[k] = self.pf.obs.data[k] - self.pf.pag.median.data[k]
                    incrl.data[k] = self.pf.pal.median.data[k] - self.pf.pb.median.data[k]
                    resl.data[k] = self.pf.obs.data[k] - self.pf.pal.median.data[k]

                innov.isloaded = True
                incrg.isloaded = True
                resg.isloaded = True
                incrl.isloaded = True
                resl.isloaded = True
                piemedinnov = Pie(innov)

                piemedinnov.plot(savefig = True, cmap = 'RdBu', clims = [-0.1, 0.1])

                piemedincrg = Pie(incrg)
                piemedincrg.plot(savefig = True, cmap = 'RdBu', clims = [-0.1, 0.1])

                # pieplot obs
                piemedresg = Pie(resg)
                piemedresg.plot(savefig = True, cmap = 'RdBu', clims = [-0.1, 0.1] )

                piemedincrl = Pie(incrl)
                piemedincrl.plot(savefig = True, cmap = 'RdBu', clims = [-0.1, 0.1])

                # pieplot obs
                piemedresl = Pie(resl)
                piemedresl.plot(savefig = True, cmap = 'RdBu', clims = [-0.1, 0.1] )

                os.chdir('..')

        if self.options.todo == 'analysis':
            if not self.options.distr:  # case semi-distributed simulation
                for dd in self.options.dates:
                    os.chdir(dd)
                    # pieplot bg
                    pb = PrepEnsBg(self.options, dd)
                    pb.calcmedian(ptinombase = 'bg')

                    # pieplot an
                    pa = PrepEnsAn(self.options, dd)
                    pa.calcmedian(ptinombase = 'an')

                    obs = Synthetic(self.xpiddir, dd, self.options)

                    _, _ = plt.subplots(len(list(pb.median.data.keys())), 3, subplot_kw=dict(polar=True))
                    plt.close()

                    '''
                    piemedb = Pie(pb.median)

                    # piemedb.plot(ax = ax[:, 0])
                    piemedb.plot()

                    piemeda = Pie(pa.median)
                    #piemeda.plot(ax = ax[:, 2])
                    piemeda.plot()
                    # pieplot obs
                    pieobs = Pie(obs)
                    #pieobs.plot(ax = ax[:, 1])
                    pieobs.plot()
                    '''
                    # for i, a in enumerate(ax[:, 0]):
                    #    a.set_ylabel([pb.median.data.keys()[i]])

                    # plt.savefig('pie/analysis.png')
                    # ######################""
                    _, _ = plt.subplots(len(list(pb.median.data.keys())), 3, subplot_kw=dict(polar=True))
                    plt.close()
                    innov = PrepAbs(dd, self.options, ptinom='innov')
                    incr = PrepAbs(dd, self.options, ptinom='incr')
                    res = PrepAbs(dd, self.options, ptinom='res')
                    innov.data = dict()
                    incr.data = dict()
                    res.data = dict()
                    obs.load()
                    for k in list(pb.median.data.keys()):
                        innov.data[k] = obs.data[k] - pb.median.data[k]

                        incr.data[k] = pa.median.data[k] - pb.median.data[k]
                        res.data[k] = obs.data[k] - pa.median.data[k]

                    innov.isloaded = True
                    incr.isloaded = True
                    res.isloaded = True
                    piemedinnov = Pie(innov)
#                     piemedinnov.plot(ax = ax2[:,0], cmap = 'RdBu')
                    piemedinnov.plot(cmap = 'RdBu', clims = [-0.1, 0.1], title = 'Innovation')

                    piemedincr = Pie(incr)
#                     piemedincr.plot(ax = ax2[:,1], cmap = 'RdBu')
                    piemedincr.plot(cmap = 'RdBu', clims = [-0.1, 0.1], title = 'Increments')
                    # pieplot obs
                    piemedres = Pie(res)

#                     piemedres.plot(ax = ax2[:,2], cmap = 'RdBu')
                    piemedres.plot(cmap = 'RdBu', clims = [-0.1, 0.1], title = 'Residuals')


#                     for i, a in enumerate(ax[:,0]):
#                         a.set_ylabel([pb.median.data.keys()[i]])
#
#                     plt.savefig('../pie/deriv.png')
#                     plt.close()

            if self.options.distr:  # case semi-distributed simulation

                for dd in self.options.dates:
                    os.chdir(dd)

                    # Read the bg

                    BG = PrepEnsBg(self.options, dd)
                    BG.stackit()

                    # Read the analysis
                    AN = PrepEnsAn(self.options, dd)
                    AN.stackit()

                    # Read the obs
                    if (self.options.synth >= 0):
                        OBS = Synthetic(self.xpidobsdir, dd, self.options)
                        OBS.load()
                    else:
                        OBS = Real(self.xpidobsdir, dd, self.options)
                        OBS.load()

                    stepZ = 100
                    minZ = np.floor(np.min(OBS.pgd.elev) / stepZ) * stepZ
                    maxZ = np.ceil(np.max(OBS.pgd.elev) / stepZ) * stepZ
                    binZ = np.arange( minZ, maxZ, stepZ)

                    BG_mean = np.zeros((len(binZ), BG.nmembers))
                    AN_mean = np.zeros((len(binZ), BG.nmembers))
                    OBS_mean = np.zeros((len(binZ)))

                    # extracts the mean OBS, BG and AN per elevation band
                    for idx, bin_Z_lower in enumerate(binZ):
                        # missing a filter on no data ! which is set to min(elev) in PGD.elev I think
                        f = np.where( (OBS.pgd.elev > bin_Z_lower) & (OBS.pgd.elev < bin_Z_lower + stepZ) )
                        tmp = OBS.data['DEP'][0, f[0]]
                        OBS_mean[idx] = np.mean(tmp)
                        for mb in range(1, BG.nmembers + 1):
                            tmp = BG.stack['DEP'][mb - 1, f]
                            BG_mean[idx, mb - 1] = np.mean(tmp)

                            tmp = AN.stack['DEP'][mb - 1, f]
                            AN_mean[idx, mb - 1] = np.mean(tmp)

                    plt.figure()
                    plt.plot(BG_mean, color=(0.4, 0.4, 0.4), alpha=0.5)
                    plt.plot(AN_mean, 'b')
                    plt.plot(OBS_mean, 'g', linewidth=4.0)
                    plt.savefig('../var_elev.png')
                    plt.close()
