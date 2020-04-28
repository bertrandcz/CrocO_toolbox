'''
Created on 6 avr. 2020

@author: cluzetb
Some useful functions for plotting
'''
from utilcrocO import dictsAspect
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


def colorbar(mappable, ax = None):
    """
    from http://joseph-long.com/writing/colorbars/
    """
    if ax is None:
        ax = mappable.axes
    else:
        ax = ax
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def cm2inch(w, h):
    return(0.393701 * w, 0.393701 * h)


def set_colors(run_name):
    """
    run names can be complicated with all the suffixes.
    set colors finds the appropriate color corresponding to the run
    """
    if 'global' in run_name:
        return 'r'
    elif 'rlocal' in run_name:
        return 'g'
    elif 'klocal' in run_name:
        if 'klocal5' in run_name:
            return 'b'
        elif 'klocal1' in run_name:
            return 'y'
        elif 'klocal' in run_name:
            return 'b'


def set_title(run_name):
    """
    run names can be complicated with all the suffixes.
    set title sets a nice title corresponding to the run
    """
    if 'global' in run_name:
        return 'global'
    elif 'rlocal' in run_name:
        return 'rlocal'
    elif 'klocal' in run_name:
        if 'klocal5' in run_name:
            return 'klocal'
        elif 'klocal1' in run_name:
            return 'klocal1'
        elif 'klocal' in run_name:
            return 'klocal'


def niceName(pgd, cl, tolist = False):
    _, revdictAsp = dictsAspect()
    return str(int(pgd.elev[cl])) + '_' + revdictAsp[pgd.aspect[cl]] + '_' + str(int(np.arctan(pgd.slope[cl]) * 180. / np.pi))


def niceLabel(var, score = None, printunits=True):

    if score is None:
        sc = ''
    else:
        sc = score
    units = {'SWE': '[$\mathrm{\mathsf{kgm^{-2}}}$]',
             'DEP': '[$\mathrm{\mathsf{m}}$]',
             'B5': '',
             'B4': '',
             }
    if printunits:
        u = units[var]
    else:
        u = ''

    ddict = {'SWE': 'SWE {0} {1}'.format(sc, u),
             'DEP': 'HS {0} {1}'.format(sc, u),
             'B5': 'Band 5 {0} {1}'.format(sc, u),
             'B4': 'Band 4 {0} {1}'.format(sc, u),
             }
    return ddict[var]
