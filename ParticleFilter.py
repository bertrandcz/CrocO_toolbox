# -*- coding: utf-8 -*-
'''
Created on 14 mars 2019

@author: cluzetb
Emulator to develop a localized particle filter in a friendly environment.
'''

from Ensemble import PrepEnsBg, PrepEnsAn
from Operators import AnalysisOperator
from Operators import PrepEnsOperator
from SemiDistributed import Synthetic
from utilcrocO import set_errors

import matplotlib.pyplot as plt
import numpy as np


class ParticleFilter(object):
    '''
    classdocs
    '''

    def __init__(self, xpiddir, options, dd):
        '''
        Constructor
        '''
        self.options = options
        self.xpiddir = xpiddir
        self.dd = dd  # specific date
        self.load()

    def load(self):
        '''
        strong duplicate from PostCroco.run()
        -> some day, mutalize
        '''

        self.pb = PrepEnsBg(self.options, self.dd)
        self.obs = Synthetic(self.xpiddir, self.dd, self.options)  # load pre-existing Synthetic (OR REAL), don't generate it.
        self.pa = PrepEnsAn(self.options, self.dd)
        self.obs.load()

    def localize(self, k=None, errobs_factor = 1., plot = False):
        """
        if k is specified, only run test for Ndes =k
        """

        # 1- Compute covariance matrix
        # - compute it
        # - mask it where the number of matching members is too low (<10)

        op = PrepEnsOperator(self.pb, sdObj=self.obs)
        corr, counts = op.bgcovariance(countIt=True)
        for var in list(corr.keys()):
            # masking corr where it is masked and when counts < 10.
            corr[var] = np.ma.masked_array(corr[var], mask = np.ma.getmask(corr[var]) | np.ma.getmask(np.ma.masked_where(counts < 10., counts)))

        if 'sd' in list(self.pb.stack.keys()):
            list(self.pb.stack.keys()).remove('sd')
            list(corr.keys()).remove('sd')

        # set global mask
        if not hasattr(self, "globMask"):
            self.setGlobMask()

        # 3 -loop on points :             | parallelizable (shared memory)
        # -> set the local mask for each var
        # -> perform the particle filter.

        # localized weights
        weights = np.empty((self.pb.ens[1].pgd.npts, self.pb.nmembers))
        resample = np.empty((self.pb.ens[1].pgd.npts, self.pb.nmembers), dtype = int)
        effweights = np.empty((self.pb.ens[1].pgd.npts,))
        if k is None:
            testlen = np.sum(np.invert(self.globMask))  # go until globmask length
            tests = list(range(1, testlen + 1))
        else:
            testlen = 1
            tests = [k]
        print(('ndesmax:', testlen))
        for ic, ndesire in enumerate(tests):
            print(('ndes:', ndesire))
            for pt in range(self.pb.ens[1].pgd.npts):

                sel = dict()   # selection depends on var.
                for var in list(self.pb.stack.keys()):  # sd has already been removed
                    corrlist = np.ma.masked_array(corr[var][pt, :], mask = self.globMask | np.ma.getmask(corr[var])[pt, :])

                    # mask the covariance when significance level of +/-0.6 is not reached and increase it iteratively
                    ceil = 0.6
                    nsel = 0.

                    while(nsel < ndesire and ceil > 0.2):
                        sel[var] = np.empty((0,), dtype = int)
                        corrlistIt = np.ma.masked_where(np.abs(corrlist) < ceil, corrlist)
                        assimMask = np.invert(np.ma.getmask(corrlistIt))  # true when correlation is Ok
                        # print('corrlistIt', corrlistIt)
                        nsel = np.sum(assimMask)
                        # print('nsel', nsel)

                        if nsel > ndesire:
                            while len(sel[var]) < ndesire:
                                imax = np.ma.masked_array.argmax(corrlistIt)
                                sel[var] = np.append(sel[var], imax)
                                corrlistIt[imax] = np.ma.masked
                            nsel = ndesire
                        else:
                            sel[var] = np.where(assimMask)[0]  # where returns a tuple
                            # print('ceil', ceil)
                        ceil = ceil / 1.2

                # perform the local filter.
                weights[pt, :], resample[pt, :] = self.pf(sel, errobs_factor=errobs_factor)
                effweights[pt] = self.efficientWeights(weights[pt, :])
            if plot:
                cmap4 = plt.cm.get_cmap(name = "jet")
                # colors = cmap4(np.linspace(0, 1, testlen))

                # old way:
                # plt.gca().set_color_cycle(colors)
                print((self.pb.nmembers))

                _ = plt.hist(effweights, bins=0.5 + np.arange(self.pb.nmembers + 1), label = 'ndes = ' + str(ndesire),
                             color = cmap4(np.float(ic) / np.float(testlen)), alpha = 0.6, zorder = 10)
            return weights, resample

    def setGlobMask(self):
        """
        set the global mask based on first reflectance variable
        """
        # 2- set global coincidence masks globMask.
        # HYP : -glob doesn't depend on variable (e.g. sd is rejected from the study)
        #       -> glob is then defined on the first key
        # coincidence mask is False when :
        # - obs is valid
        # - AND ens is valid (e.g.more than 10 members or 10% of members are defined)
        # - AND  and obs-ens are comparable (no condition on that so far)

        maskVar = list(self.pb.stack.keys())[0]  # use mask of first reflectance variable
        print(('mask var', maskVar))
        print(( 'shape pb', np.shape(self.pb.stack[maskVar])))
        self.globMask = np.ma.getmask(np.ma.masked_invalid(self.obs.data[maskVar])) | (np.ma.count(np.ma.masked_invalid(self.pb.stack[maskVar]), axis = 1) < 10.)

    def globalize(self, errobs_factor = 1., plot = False):

        if not hasattr(self, 'globMask'):
            self.pb.stackit()
            self.setGlobMask()

        # global weights
        gsel = dict()
        # perform the global filter
        for var in list(self.pb.stack.keys()):
            gsel[var] = np.where(np.invert(self.globMask))[0]

        # print('gsel', gsel, np.sum(np.invert(globMask)))
        gweights, resample = self.pf(gsel, errobs_factor=errobs_factor)
        gresample = np.tile(resample.T, (self.pb.ens[1].pgd.npts, 1))
        print(('greshp', gresample.shape))
        geffweights = self.efficientWeights(gweights)
        print(('global assim (npts, Neff):', np.sum(np.invert(self.globMask)), geffweights))

        if plot:
            ylims = plt.gca().get_ylim()
            plt.bar(geffweights - 0.25, ylims[1], width = 0.5, zorder = 0, color = 'k', label= 'global, Npts: ' + str(np.sum(np.invert(self.globMask))))
            # plot results

            plt.legend().set_zorder(20)

            plt.xlim([0.5, self.pb.nmembers + 0.5])
            plt.ylim(ylims)
            plt.ylabel('Number of occurences', fontsize = 25)
            plt.xlabel('Neff', fontsize = 25)

        return gweights, gresample

    def pf(self, sel, errobs_factor = 1.):
        """
        perform the actual particle filter on a selection of points
        sel is a dict containing the assimilation points for each var.

        """

        # prepare observation and model vectors:
        listObs = []

        # prepare for err

        for var in list(self.pb.stack.keys()):

            listObs = listObs + list(map(int, sel[var].tolist()))
        if len(listObs) > 0:
            vectObs = np.empty((0,))
            matEns = np.empty((0, self.pb.nmembers ))  # 35 columns forced, put 0. when ensemble is not def.
            r_1 = np.empty((0,))
            for var in list(self.pb.stack.keys()):
                vectObs = np.hstack((vectObs, np.take(self.obs.data[var], sel[var], axis = 0)))

                r_1 = np.append(r_1, 1. / (set_errors([var])[0]) * np.ones((len(sel[var]),)))   # set errors must receive a list
                matEns = np.ma.masked_invalid(np.vstack((matEns, np.take(self.pb.stack[var], sel[var], axis = 0))))    # take the sel[var] lines of stack

            R_1 = np.diag(r_1) / errobs_factor
            # replace NaNs by 0.
            matEns[np.where(np.isnan(matEns))] = 0.

            # print('vectObs', vectObs)
            likelihoods = np.array([np.exp(-0.5 * np.matmul(np.matmul((matEns[:, i] - vectObs).T, R_1), (matEns[:, i] - vectObs)))
                                    for i in range(self.pb.nmembers)])
            weights = likelihoods / np.sum(likelihoods)

            # print(np.shape(weights))
            # print(weights)
        else:
            weights = np.ones((self.pb.nmembers,)) / self.pb.nmembers
        resample = self.resampling(weights)
        return weights, resample

    def resampling(self, weights):
        """
        PF resampling, replicated from assim_nature_isba_pf.F90
        /!\ resample is in python indices, e.g. 0 for member 1
        """
        resample = np.empty((self.pb.nmembers,), dtype = int)
        zweightcumul = np.cumsum(weights)
        zrand = np.random.random() / self.pb.nmembers
        for ires in range(0, self.pb.nmembers):
            if zrand <= zweightcumul[0]:  # zweightcumul[0] = poids du mb 1
                resample[ires] = 1
            for rk in range(1, self.pb.nmembers):
                if (zweightcumul[rk - 1] < zrand) and (zrand <= zweightcumul[rk]):
                    resample[ires] = rk + 1

            zrand += 1. / self.pb.nmembers
        # print(weights)
        # print(resample)
        # plt.hist(resample, bins = 0.5 + np.arange(self.pb.nmembers + 1))
        # ax2 = plt.gca().twinx()
        # ax2.plot(range(1, self.pb.nmembers + 1), weights, color = 'k', marker = '+', linestyle = 'None', markersize = 10)
        # plt.show()
        # !!!!!!!!!!!!!!!!!!!!! -1 for python indexing
        print((resample - 1))
        return resample - 1    # -1 for pythonindexing indexing !!!!!!!!!!!!!!!

    def analyze(self, weightsMat):

        anOp = AnalysisOperator(self.pb, weightsMat)
        an = anOp.analyze()
        return an

    def efficientWeights(self, weights):
        # reinvented formula 20/03
        effweights = np.around(1. / (self.pb.nmembers * np.mean([w**2 for w in weights])))

        return effweights
