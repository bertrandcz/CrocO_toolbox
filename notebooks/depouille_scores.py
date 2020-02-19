
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import difflib
from utilcrocO import cm2inch
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.patches as mpatches


# In[34]:


score = 'CRPS'
kind = 'quant'
nens = 40
dep = 'DEP'
if nens == 160:
    suffix = ''
else:
    suffix = '_{0}'.format(nens)
if dep =='DEP':
    dep = '_DEP'
suffix = suffix + dep
suffixName=suffix
#suffixName = '_40_DEP_neffstudy'
#suffixName = suffix + 'k2013_2014'
name = 'dfScores_{0}_{1}{2}.pkl'.format(score, kind, suffixName)
dfScores = pd.read_pickle(name)
#dfScores
if dep == '':
    dfScores = pd.read_pickle('dfScores_CRPS_quant_40.pkl')
    dfScores160 = pd.read_pickle('dfScores_CRPS_quant.pkl')
dfScores


# In[35]:


dc = {'rlocal': 'g', 'klocal5': 'b', 'global': 'r',
          'ol': '0.7', 'cl': 'm', 'klocal10': 'c', 'klocal1': 'y'}
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
years = list(dfScores.loc[pd.IndexSlice[:, :, 'f'], :].index.levels[0])


# In[36]:


sscores = ['CRPS', 'Reli']

def set_col_names(df, dc, s, suffix):
    algs = ['ol'] + sorted(list(set([k for k in dc.keys()
      if (len([c.startswith(k) for c in df.columns if c.startswith(k)]) > 0 and k is not 'ol')])))
    colnames = [a + suffix if a is not 'ol' else a for a in algs]

    cols_scores = [c + '_' + s for c in colnames]
    cols_r = [c + '_r' for c in cols_scores if 'ol' not in c]
    return algs, colnames, cols_scores, cols_r
def compute_rel(df, cols, s):
    for col in cols:
        if 'ol' not in col:
            df[col + '_r'] = 1 - df[col] / df['ol' + '_' + s]
    return df
def boxplot(df, fc, cols, pos, dc, ax, hatched = False):
    bp_dict_D = dict()
    for i, col in enumerate(cols):
        colkey = [k for k in dc.keys() if col.startswith(k)][0]
        bp_dict_D[col] = ax.boxplot(df.loc[pd.IndexSlice[:, :, fc], :][col],0,'k+',
                                    positions = [pos[i]],
                                    boxprops=dict(facecolor=colors(dc, col, colkey) if hatched is False else 'white',
                                                  color = colors(dc, col, colkey),
                                                  hatch = '////' if hatched is True else None,
                                      #showcaps = False, showfliers = False, whiskerprops = dict(color = 'None'),
                                      #whis = [5, 95], widths = 3,
                                      ),
                         patch_artist=True, )
    return bp_dict_D


# In[37]:




fig, axs = plt.subplots(nrows=4, ncols=2, sharey='row',
                        figsize=(cm2inch(8.3, 13)))
plt.subplots_adjust(wspace=0.05 if score == 'CRPS' else 0.33,
                    top=0.97, left=0.20, bottom=0.2, hspace=0.2, right=0.99)
for iS, sc in enumerate(sscores):
    algs, colnames, cols, cols_r = set_col_names(dfScores, dc, sc, suffix)
    dfScores =compute_rel(dfScores, cols, sc)
    if dep =='':
        algs160, colnames160, cols160, cols_r160 = set_col_names(dfScores160, dc, sc,'')

        dfScores160 = compute_rel(dfScores160, cols160, sc)
    for iF, fc in enumerate(['f', 'nf']):
        # normal score
        ax = axs[2*iS, iF]

        boxplot(dfScores, fc, cols, np.arange(0.2, 2.2,0.6), dc, ax)
        if dep == '':
            boxplot(dfScores160, fc, cols160, np.arange(0.4, 2.5,0.6), dc, ax, hatched = True)
        
        ax.set_xlim([0,2.4])
        ax.set_xticklabels([])
        ax.set_xticks([])

        ax.set_ylim([0, 170])
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=2)
        if iF == 0:
            ax.set_ylabel('SWE %s [$\mathrm{\mathsf{kgm^{-2}}}$]' %
                          sc, fontsize=9)

        ax.grid(which = 'minor', alpha = 0.3)
        ax.grid(which = 'major')        # relative score
        ax = axs[2*iS+1, iF]
        # set relative scores

        boxplot(dfScores, fc, cols_r, 0.6+ np.arange(0.2, 2.2,0.6),dc, ax)
        ax.axhline(y=0, color = 'k', ls = '--')
        if dep == '':
            boxplot(dfScores160, fc, cols_r160, 0.6+np.arange(0.4, 2.5,0.6), dc,ax, hatched = True)
        ax.set_xlim([0,2.4])
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_ylim([-2, 1.])
        ax.set_yticks([-2, -1, 0, 1])
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=2)
        ax.axhline(y=1, color='0.7', ls='--', zorder=-1)
        if iF == 0:
            ax.set_ylabel('SWE %sS' % sc, fontsize=9)
        if iS == 1:
            ax.text(0.5, -0.25, 'observed' if fc == 'f' else 'not observed',
                    va='center', ha='center',
                    fontsize=12,
                    transform=ax.transAxes)
        ax.grid(which = 'minor', alpha = 0.3)
        ax.grid(which = 'major')
patches = []
for col in cols:
    colkey = [k for k in dc.keys() if col.startswith(k)][0]
    patches.append(mpatches.Patch(
        color=colors(dc, col, colkey), label=colkey))
fig.legend(handles=patches, bbox_to_anchor=(0.5, 0.08), ncol=2,
           bbox_transform=fig.transFigure, loc='center')
fig.savefig('aggr_DEP.pdf')
fig.savefig('aggr_DEP.png')


# In[6]:


dfScores160


# In[7]:


pwd


# In[ ]:





# In[8]:


sscores = ['CRPS', 'Reli']

for sc in sscores:
    fig, axs = plt.subplots(nrows=2, ncols=2, sharey='row',
                            figsize=(cm2inch(8.3, 8.3)))
    plt.subplots_adjust(wspace=0.05 if score == 'CRPS' else 0.33,
                        top=0.97, left=0.20, bottom=0.3, hspace=0.2, right=0.99)
    for i, fc in enumerate(['f', 'nf']):
        # normal score
        ax = axs[0, i]
        algs = ['ol'] + sorted(list(set([k for k in dc.keys()
                                         if (len([c.startswith(k) for c in dfScores.columns if c.startswith(k)]) > 0 and k is not 'ol')])))
        colnames = [a + suffix if a is not 'ol' else a for a in algs]
        cols = [c + '_' + sc for c in colnames]

        def boxplot(dfScores, fc, cols, ax, labs):
            bp_dict = dfScores.loc[pd.IndexSlice[:, :, fc], :].boxplot(column=cols, ax=ax, return_type='both',
                                                                       patch_artist=True, )
            for i, box in enumerate(bp_dict[1]['boxes']):
                colkey = [k for k in dc.keys() if cols[i].startswith(k)][0]
                plt.setp(box, color=colors(dc, cols[i], colkey))
            for i, box in enumerate(bp_dict[1]['medians']):
                colkey = [k for k in dc.keys() if cols[i].startswith(k)][0]
                plt.setp(box, color='k')
        boxplot(dfScores, fc, cols, ax, algs)
        ax.set_xticklabels([])
        ax.set_ylim([0, 170])
        # ax.set_yticks(range(0,175,25))
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=2)
        if i == 0:
            ax.set_ylabel('SWE %s [$\mathrm{\mathsf{kgm^{-2}}}$]' %
                          sc, fontsize=9)


        # skill score       
        ax = axs[1, i]
        for col in cols:
            if 'ol' not in col:
                dfScores[col + '_r'] = 1 - dfScores[col] /                     dfScores['ol' + '_' + sc]
        cols_r = [c + '_r' for c in cols if 'ol' not in c]

        
        
        ax.axhline(y=1, color='0.7', ls='--', zorder=-1)

        labs = [c for c in algs if 'ol' not in c]
        boxplot(dfScores, fc, cols_r, ax, labs)
        ax.set_xticklabels([])
        ax.set_ylim([-2, 1.1])
        ax.set_yticks([-2, -1, 0, 1])
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=2)
        ax.axhline(y=1, color='0.7', ls='--', zorder=-1)
        if i == 0:
            ax.set_ylabel('SWE %sS' % sc, fontsize=9)
        ax.text(0.5, -0.25, 'observed' if fc == 'f' else 'not observed',
                va='center', ha='center',
                fontsize=12,
                transform=ax.transAxes)
        patches = []
        for col in cols:
            colkey = [k for k in dc.keys() if col.startswith(k)][0]
            patches.append(mpatches.Patch(
                color=colors(dc, col, colkey), label=colkey))
    axs[0, 0].yaxis.set_label_coords(-0.285, 0.5)
    axs[1, 0].yaxis.set_label_coords(-0.31, 0.5)
    fig.legend(handles=patches, bbox_to_anchor=(0.5, 0.1), ncol=2,
               bbox_transform=fig.transFigure, loc='center')


# In[9]:



if kind == 'baseline':
    dict_q = {i: '{0}'.format(i[0])
              for i in dfScores.loc[pd.IndexSlice[:, :, 'f'], :].index}
else:
    dict_q = {i: 'q{1}'.format(i[0], quant)
              for (i, quant) in zip(dfScores.loc[pd.IndexSlice[:, :, 'f'], :].index, 
                                    list(range(20, 100, 20))*len(years))}

xlabels = ['RMSE', 'spread', 'ss'] if score == 'RMSE' else [
    'CRPS', 'Reli', 'Resol']
fig, axs = plt.subplots(nrows=2, ncols=3, sharex=True,
                       sharey=True if score == 'CRPS' else False, figsize = (15,5))
plt.subplots_adjust(wspace=0.05 if score =='CRPS' else 0.33, top=0.97, left=0.15, bottom=0.25, right=0.8)
for i, fc in enumerate(['f', 'nf']):
    for j, sc in enumerate(xlabels):
        ax = axs[i, j]
        sel = [c for c in dfScores.columns if sc in c]
        for col in sel:
            colkey = [k for k in dc.keys() if col.startswith(k)][0]
            ax.plot(range(len(dict_q.keys())),
                    dfScores.loc[pd.IndexSlice[:, :, fc], :][col],
                    #ax = ax,
                    marker='o',
                    markersize =4,
                    alpha=1 if col.startswith('ol') else 0.8,
                    zorder=0 if col.startswith(
                        'ol') or col.startswith('cl') else 1,
                    lw=3 if col.startswith('ol') else 1,
                    color=colors(dc,col, colkey),
                    label=colkey
                    )
        if kind =='quant':
            for iy in range(len(years)-1):
                ax.axvline(x = 3.5 + (4*iy), ls = '--', color = 'k', zorder = 0)
        if i == 1:
            ax.text(0.5, -0.55, xlabels[j], va ='center', ha = 'center', fontsize=16, transform=ax.transAxes)
            ax.set_xticks(range(len(dict_q.keys())))
            ax.set_xticklabels(
                [dict_q[i] for i in dfScores.loc[pd.IndexSlice[:, :, 'f'], :].index], rotation=90)
            for iy, year in enumerate(years):
                ax.text(1./(2.*len(years))+iy*1./len(years), -0.38, str(year),
                        va = 'center', ha = 'center',
                        fontsize = 14,
                        transform=ax.transAxes)
        if j == 2:
            if score == 'RMSE':
                ax.set_ylim([0, 3])
                ax.axhline(y=1, ls='--', zorder=-2, color='0.8')
        else:
            if score == 'RMSE':
                ax.set_ylim([0, 300])
            if score == 'CRPS':
                ax.set_ylim([0, 220])
axs[0, 0].set_ylabel('Observed')
axs[1, 0].set_ylabel('Not Observed')
lu = ax.legend(loc='center left', bbox_to_anchor=(1, 1))
gg = fig.text(0.01, 0.65, 'SWE %s [$\mathrm{\mathsf{kgm^{-2}}}$]' %
              score, va='center', rotation='vertical', fontsize=16)

fig.savefig('depouille_{0}_{1}{2}.png'.format(score, kind, suffixName), dpi=200)


# In[10]:


pwd


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




