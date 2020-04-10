#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on 6 févr. 2019
Objects describing ensembles of Prep files.
@author: cluzetb
'''
from SemiDistributed import PrepBg, PrepAn, PrepAbs
from SemiDistributed import PrepOl
import numpy as np


class Ensemble(object):
    '''
    class to describe an ensemble of objects.
    '''
    _abstract = True


class PrepEns(Ensemble):
    """
    ensemble of prep.
    Data is stored twice and instanciating takes time : to ameliorate in the future.
    """
    _abstract = True

    def __init__(self, options, date):
        self.options = options
        self.date = date
        self.nmembers = options.nmembers
        self.ens = dict()
        self.stack = dict()
        self.isloaded = False
        self.isstacked = False

    def load(self):
        for mb in range(1, self.nmembers + 1):
            self.ens[mb].load()
        self.listvar = []
        for var in list(self.ens[1].data.keys()):
            self.listvar.append(var)
        self.pgd = self.ens[1].pgd
        self.isloaded = True

    def stackit(self):
        # for median calc and manipulation, to stack the members.
        if not self.isstacked:
            # self.data = dict()
            if not self.isloaded:
                self.load()

            for var in self.listvar:
                # print ('shape before transposing stack', np.shape(self.ens[1].data[var]))
                self.stack[var] = self.ens[1].data[var]
                for mb in range(2, self.nmembers + 1):
                    self.stack[var] = np.vstack((self.stack[var], self.ens[mb].data[var]))
                self.stack[var] = self.stack[var].T  # change in the dimension of the stack (pts as lines, mbs as cols)
                # self.stack[var] = np.ma.masked_where(self.stack[var] == self.stack[var][:].fill_value, self.stack[var])
        self.isstacked = True

    def calcmedian(self, ptinombase = ''):
        '''
        truc un peu batard, je surcharge en qq sorte la méthode median
        '''
        self.stackit()
        self.median = PrepAbs(self.date, self.options, ptinombase)
        self.median.ptinom = ptinombase + 'med'  # overwrite
        self.median.data = dict()
        for var in self.listvar:
            self.median.data[var] = np.ma.median(self.stack[var], axis = 1)
        self.median.isloaded = True


class PrepEnsOl(PrepEns):

    def __init__(self, options, date, directFromXp = True, isOl = False):
        PrepEns.__init__(self, options, date)
        for mb in range(1, self.nmembers + 1):
            self.ens[mb] = PrepOl(date, mb, options, directFromXp = directFromXp, isOl = isOl)


class PrepEnsBg(PrepEns):

    def __init__(self, options, date, directFromXp = True):
        PrepEns.__init__(self, options, date)

        for mb in range(1, self.nmembers + 1):
            self.ens[mb] = PrepBg(date, mb, options, directFromXp = directFromXp)


class PrepEnsAn(PrepEns):

    def __init__(self, options, date, directFromXp = True):
        PrepEns.__init__(self, options, date)
        for mb in range(1, self.nmembers + 1):
            self.ens[mb] = PrepAn(date, mb, options, directFromXp = directFromXp)


class PrepEnsAbs(PrepEns):

    def __init__(self, options, date):
        PrepEns.__init__(self, options, date)
        for mb in range(1, self.nmembers + 1):
            self.ens[mb] = PrepAbs(date, options, ptinom = 'abs' + str(mb))
