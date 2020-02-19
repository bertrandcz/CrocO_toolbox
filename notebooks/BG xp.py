
# coding: utf-8

# In[5]:


import time
from utilpp import read_BG
get_ipython().run_line_magic('matplotlib', 'inline')


# In[6]:


get_ipython().run_line_magic('run', '/home/cluzetb/assim/scripts/pp_SWECRPSfrompickle_time.py')


# In[7]:


bg50 = read_BG(RUN['2016_50_klocal5'].options)
bg153 = read_BG(RUN['2016_153_klocal5'].options)


# In[9]:


for dd in bg.keys():
    plt.close('all')
    fig, axs = plt.subplots(nrows = 2, ncols = 2,figsize =(5,5))
    axs[0,0].imshow(bg50[dd]['B4'], cmap = 'RdBu', vmin = -1, vmax = 1)
    axs[0,1].imshow(bg50[dd]['B5'], cmap = 'RdBu', vmin = -1, vmax = 1)
    axs[1,0].imshow(bg153[dd]['B4'], cmap = 'RdBu', vmin = -1, vmax = 1)
    axs[1,1].imshow(bg153[dd]['B5'], cmap = 'RdBu', vmin = -1, vmax = 1)
    fig.suptitle(dd)
    plt.show()


# In[ ]:




