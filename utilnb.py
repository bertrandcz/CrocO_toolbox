#! /usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on April 2020
Some useful fonction for factorization of notebooks using the pandas dataframes holding the scores (dfScores*)
@author: cluzetb
'''
import pandas as pd


def colors(dc, col, colkey):
    if 'global' not in col:
        ret = dc[colkey]
    elif '_1_' in col:
        ret = 'm'
    elif '_10_' in col:
        ret = 'b'
    else:
        ret = dc[colkey]
    return ret


def colors_legs(xp):
    if xp.endswith('_10'):
        col = 'm'
        lab = '$N_{eff}^* = 10$'
    elif xp.endswith('_1'):
        col = 'darkorange'
        lab = '$N_{eff}^* = 1$'
    else:
        col = 'royalblue'
        lab = '$N_{eff}^* = 7$'
    return col, lab


def set_col_names(df, dc, s, suffix):
    algs = ['ol'] + sorted(list(set([k for k in dc.keys()
                                     if (len([c.startswith(k) for c in df.columns if c.startswith(k)]) > 0 and k != 'ol')])))
    colnames = [a + suffix if a != 'ol' else a for a in algs]

    cols_scores = [c + '_' + s for c in colnames]
    cols_r = [c + '_r' for c in cols_scores if 'ol' not in c]
    return algs, colnames, cols_scores, cols_r


def boxplot(df, fc, cols, pos, dc, ax, hatched = False):
    """
    boxplot columns of a dfScores.
    """

    bp_dict_D = dict()
    for i, col in enumerate(cols):
        colkey = [k for k in dc.keys() if col.startswith(k)][0]
        bp_dict_D[col] = ax.boxplot(df.loc[pd.IndexSlice[:, :, fc], :][col], 0, 'k+',
                                    positions = [pos[i]],
                                    boxprops=dict(facecolor=colors(dc, col, colkey) if hatched is False else 'white',
                                                  color = colors(dc, col, colkey),
                                                  hatch = '////' if hatched is True else None,
                                                  # showcaps = False, showfliers = False, whiskerprops = dict(color = 'None'),
                                                  # whis = [5, 95], widths = 3,
                                                  ),
                                    patch_artist=True, )
    return bp_dict_D


# function to compute skill scores
def compute_rel(df, cols, s):
    '''
    this function computes a skill score wrt the openloop
    '''
    for col in cols:
        if 'ol' not in col:
            df[col + '_r'] = 1 - df[col] / df['ol' + '_' + s]
    return df
