# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:51:07 2019

@author: cluzetb
"""
import copy
import time

from scipy import stats

from CramponPp import CramponPp
from Operators import PrepEnsOperator
from SemiDistributed import PrepAbs
from crampon import set_options
import matplotlib.pyplot as plt
import numpy as np
from plotcrampon import Pie, spaghettispace
from utilcrampon import setSubsetclasses


plt.close('all')
start_time = time.time()
kinds = {  # 'klocal_f10': 'art2_klocal_f10_infl_t1500@cluzetb',
    # 'klocal_f1': 'art2_klocal_f1_infl_t1500@cluzetb',
    # 'klocal': 'art2_klocal_infl_t1500@cluzetb',
    'global': 'art2_2013_global_infl_t1500@cluzetb',
    # 'rlocal': 'art2_rlocal_infl_t1500@cluzetb'
}

dcolors = {'klocal_f10': 'c', 'klocal_f1': 'm', 'klocal': 'b', 'rlocal': 'g', 'global': 'r'}

for ik, key in enumerate(sorted(kinds.keys())):
    args = [
        '/home/cluzetb/snowtools_git/assim/crampon.py',
        '--xpid', 'OL_youness_160',
        '--vconf', 'lautaretreduc',
        '-d', 'all',
        '--vars', 'B2,B4,B5,R54,R52,SCF',
        '--ppvars', 'B2,B4,B5,DEP,SWE,R54,R52',
        '-o', 'compar2to3',
    ]
    options, conf = set_options(args)

    run = CramponPp(options, conf)
    # run.readEns()


def whatrr():
    crpsa = dict()
    for dd in options.dates:
        opb = PrepEnsOperator(run.ensBg[dd], sdObj=run.obsArch[dd])
        ccb = opb.CRPS()
        opa = PrepEnsOperator(run.ensAn[dd], sdObj=run.obsArch[dd])
        cca = opa.CRPS()
        optcp = copy.copy(options)
        optcp.classesS = ['0', '20']
        optcp.classesE = 'all'
        optcp.classesA = 'all'
        optcp.synth = 11

        crpsa[dd] = PrepAbs(dd, optcp, ptinom='crpsa')

        data = dict()
        for k in options.ppvars:

            data[k] = 1 - cca[k] / ccb[k]

        crpsa[dd].data = data
        crpsa[dd].isloaded = True

        focusCl = setSubsetclasses(run.obsArch[dd].pgd, options.classesE,
                                   options.classesA, options.classesS)[0]
        percentiles = np.array([stats.percentileofscore(run.ensBg[dd].stack['DEP'][cl, :], run.obsArch[dd].data['DEP'][cl]) for cl in focusCl])
        quartiles = np.floor(percentiles / 25.) / 3.  # /3 for colors

        obsArchDep = np.array([run.obsArch[dd].data['DEP'][cl] for cl in focusCl])
        percentilesCol = np.where(obsArchDep < 0.000001, -100., quartiles)
        piecrpsa = Pie(crpsa[dd],
                       focusCl = focusCl,
                       focusColors = percentilesCol,)
        piecrpsa.plot(cmap = 'RdBu', clims = [-1, 1],
                      title = 'an_CRPSS (' + options.pf + '), ' + dd,
                      colortitle = 'CRPSS', savefig = True)

        spaghettispace('DEP', run.ensBg[dd], run.ensAn[dd], obs = run.obsArch[dd], savePath='spag_' + 'DEP' + dd + '_' + '.png')


elapsed_time = time.time() - start_time
print(elapsed_time)
