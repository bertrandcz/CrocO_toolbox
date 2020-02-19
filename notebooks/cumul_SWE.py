
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import datetime
from scipy import stats
import matplotlib.pyplot as plt


# In[2]:


get_ipython().run_line_magic('run', '/home/cluzetb/assim/scripts/load_multi_xp.py')


# In[3]:


sumSWE = np.sum(RUN[xp].ensProOl['SWE'], axis=(0,1))
gg = np.array([stats.percentileofscore(sumSWE, s) for s in sumSWE])
gg[8]


# In[4]:


perc =  range(20,100,20)
mins = [np.argmin(abs(gg-p)) for p in perc]


# In[5]:


print(RUN[xp].ensProOl['SWE'].shape)
sumSWE_time = np.cumsum(np.sum(RUN[xp].ensProOl['SWE'], axis =1), axis = 0)
print(sumSWE_time.shape)


# In[6]:


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
# plot the openloop
[plt.plot(RUN[xp].ensProOl['time'],sumSWE_time[:,i], color ='0.7', alpha = 0.2, label = 'openloop' if i==0 else '' )
 for i in range(sumSWE_time.shape[1])]
ax.plot(RUN[xp].ensProOl['time'], np.median(sumSWE_time, axis=1), color = 'k', lw = 2, label ='median')
[ax.plot(RUN[xp].ensProOl['time'],sumSWE_time[:,mb], lw = 2, label = 'mb{0:04d}, perc:{1}'.format(mb+1,p)) for (mb,p) in zip(mins,perc)]
ax.set_ylabel('integrated SWE [$kgm^{-2}$]')
l = ax.legend()
fig.savefig('integrated_SWE_' + str(RUN[xp].ensProOl['time'][0].year) + '.png', bbox_inches = 'tight')


# In[ ]:





# In[ ]:




