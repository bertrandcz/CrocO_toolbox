#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 8 févr. 2019

@author: cluzetb
'''

from Ensemble import PrepEnsAbs
from SemiDistributed import PrepAbs
import numpy as np
from scores.ensemble import EnsembleScores


class SdOperator(object):
    '''
    basic class for operations on SemiDistributed objects
    either from real or abstract ones
    building generally abstract ones
    '''

    _abstract = True

    def __init__(self):
        '''
        Constructor
        '''


class EnsOperator(object):
    '''
    basic class for operations on Ensemble objects

    '''
    _abstract = True

    def __init__(self):
        '''
        Constructor
        '''


class PrepEnsOperator(EnsOperator):
    '''
    basic class for operations on PrepEns objects
    '''

    def __init__(self, ens1, ens2=None, sdObj = None):
        '''
        Constructor
        '''
        EnsOperator.__init__(self)
        self.ens1 = ens1
        self.ens1.stackit()
        if not hasattr(self.ens1, 'listvar'):
            self.ens1.listvar = list(self.ens1.stack.keys())
        if ens2 is not None:
            self.ens2 = ens2
            self.ens2.stackit()
            if not hasattr(self.ens2, 'listvar'):
                self.ens2.listvar = list(self.ens2.stack.keys())
        if sdObj is not None:
            self.sdObj = sdObj
            self.sdObj.load()

    def m2mdiff(self):
        # first, stack ensembles (more natural for computations)
        diff = PrepEnsAbs(self.ens1.options, self.ens1.dates)
        diff.listvar = self.ens1.listvar
        for var in self.diff:
            diff.stack[var] = self.ens1.stack[var] - self.ens2.stack[var]

    def m21diff(self, reverse = False):
        """
        compute ens - sdObj (reverse for innovations) and put result into ens2
        """
        fact = 1
        if reverse:
            fact = -1
        self.ens2 = PrepEnsAbs(self.ens1.options, self.ens1.date)
        self.ens2.listvar = self.ens1.listvar
        for var in self.ens2.listvar:
            self.ens2.stack[var] = fact * (self.ens1.stack[var] - self.sdObj.data[var])
        self.ens2.isloaded = True
        self.ens2.isstacked = True
        return self.ens2

    def CRPS(self):
        '''
        compute CRPS between the obs. and the ensemble.
        '''
        crps = dict()
        for var in self.ens1.listvar:
            crps[var] = np.array([EnsembleScores(1, 1, np.array([self.sdObj.data[var][pt]]),
                                                 np.expand_dims(self.ens1.stack[var][pt, :], axis=1)).CRPS()
                                  for pt in range(np.shape(self.ens1.stack[var])[0])])
        return crps

    def bgcovariance(self, clA=None, countIt=False):
        '''
        compute covariance in bg ensemble ens1
        if clA, restrict it to covariance between ens1 in clA and ens1 in all classes (then it is plottable)
        /!\ BC 18/12/19 : deprecated (need a transpose everywhere!!)

        '''
        if clA is None:
            # TODO : write full covariance computation.
            corr = dict()
            for var in self.ens1.listvar:
                ensMat = np.ma.masked_invalid(self.ens1.stack[var])  # (npts, nens)
                """ useless parallel (20/03/19)
                start_1 = time.time()
                corr[var] = np.array(compute_corr_parallel(ensMat))
                """
                corr[var] = np.ma.corrcoef(ensMat)  # removed the transpose bcz now the stack ois (npts, nens) (more natural)

            if countIt is True:
                counts = dict()

                # /!\!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!only compute for first band (same for all) be careful here...
                for var in self.ens1.listvar[0]:   # only compute for first band (same for all) be careful here...

                    # print('test1', np.shape(np.sum(np.invert(np.ma.getmask(ensMat[140, :]) | np.ma.getmask(ensMat[141, :])))))
                    counts[var] = np.array([[np.sum(np.invert(np.ma.getmask(ensMat[i, :]) | np.ma.getmask(ensMat[j, :])))
                                             for i in range(np.shape(ensMat)[0])] for j in range(np.shape(ensMat)[0])])
                    # print('counts')
                    # print(counts)
                    # print np.shape(counts)

                    # plt.figure()
                    # plt.imshow(counts[var])
                    # plt.show()

                return corr, counts
            else:
                return corr
        else:
            corr_cl = PrepAbs(self.ens1.date, self.ens1.options, ptinom='corr_')
            corr_cl.data = dict()
            for var in self.ens1.listvar:
                corr_cl.data[var] = np.ma.empty((np.shape(self.ens1.stack[var])[0]))
                for cl in range(np.shape(self.ens1.stack[var])[0]):
                    gg1 = np.ma.masked_invalid(self.ens1.stack[var][clA, :])
                    gg2 = np.ma.masked_invalid(self.ens1.stack[var][cl, :])
                    corr_cl.data[var][cl] = np.ma.corrcoef(gg1, gg2)[1, 0]
            corr_cl.isloaded = True
            return corr_cl


class AnalysisOperator(PrepEnsOperator):
    """
    Operator to perform a soda-like analysis,
    based on a PrepEnsBg and list of matrix to replicate. (Npts, nmembers)
    """

    def __init__(self, ens1, resampleMat):
        PrepEnsOperator.__init__(self, ens1, ens2 = None, sdObj = None)
        self.weightsMat = resampleMat
        self.ens2 = PrepEnsAbs(self.ens1.options, self.ens1.date)
        self.ens2.listvar = self.ens1.listvar

    def analyze(self):

        for var in self.ens2.listvar:

            self.ens2.stack[var] = np.array([np.squeeze(self.ens1.stack[var][pt, self.weightsMat[pt, :]]) for pt in range(self.ens1.stack[var].shape[0])])

        self.ens2.isloaded = True
        self.ens2.isstacked = True
        return self.ens2
