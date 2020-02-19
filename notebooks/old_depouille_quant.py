
# coding: utf-8

# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# In[6]:


score = 'RMSE'
kind = 'quant'
name = 'dfScores_{0}_{1}.pkl'.format(score, kind)

dfScores = pd.read_pickle(name)


# In[7]:


dfScores


# In[9]:


dc = {p : c for p,c in zip(dfScores.columns, ['g','b','r','k','m'])}
dict_q = {i:'{0}_q{1}'.format(i[0],quant) 
      for (i,quant) in zip(dfScores.loc[pd.IndexSlice[:,:,'f'],:].index,range(20,100,20))}

xlabels = ['RMSE', 'spread', 'ss'] if score == 'RMSE' else [
    'CRPS', 'Reli', 'Resol']

fig, axs = plt.subplots(nrows = 2, ncols = 3,figsize = (20,10), sharex = True, sharey = True)
plt.subplots_adjust(wspace = 0.09,top = 0.97, left = 0.1, bottom = 0.3, right = 0.88)
plt.ylim([0,150])
for i, fc in enumerate(['f', 'nf']):
    for j, sc in enumerate(xlabels):
        sel = [c for c in dfScores.columns if sc in c]
        
gg = fig.text(0.01,0.65,'SWE CRPS [$\mathrm{\mathsf{kgm^{-2}}}$]', va = 'center', rotation = 'vertical', fontsize = 16)

fig.savefig('depouille_{0}_{1}.png'.format(score, kind))


# In[ ]:




