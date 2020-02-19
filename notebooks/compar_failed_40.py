
# coding: utf-8

# In[45]:


get_ipython().run_line_magic('matplotlib', 'inline')
from plotcrocO import prepare_spaghetti_Marie
get_ipython().run_line_magic('run', '/home/cluzetb/assim/scripts/load_multi_xp.py postesdebug')


# In[46]:


focusCl = setSubsetclasses( RUN[xp].pgd, options.classesE, options.classesA, options.classesS)[0]
nfocusCl = [p for p in range(RUN[xp].pgd.npts) if p not in focusCl ]


# In[47]:


#for cl in range(RUN[xp].pgd.npts):
itimes, itimesC = set_itimes(RUN[xp], clim=True, fromOl=True)
itimesAn = set_itimes(RUN[xp])
xtimes = RUN[xp].ensProOl['time'][itimes]
for cl in focusCl:
    prepare_spaghetti_Marie(RUN[xp], 'SWE', cl)
    plt.plot(xtimes,RUN[xp].ensProOl['SWE'][itimes,cl, 65], lw = 3, color = 'g')
plt.show()


# In[31]:


RUN[xp].mbsynth


# In[26]:


from plotcrocO import snsscatter_2vars


# In[30]:


date = '2013122910'
for cl in focusCl:
    snsscatter_2vars(RUN[xp].options,RUN[xp].ensBg[date], RUN[xp].ensAn[date], RUN[xp].obsArch[date], 'DEP', 'SWE', cl)


# In[36]:


RUN[xp].ensProOl['SWE'].shape


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[8]:


from scipy.stats import norm
cl = 124
dd = '2017041410'
prepare_spaghetti_Marie(RUN[xp], 'SWE', cl)

for var in ['SWE', 'B4','B5']:
    if dd not in list(RUN['2016_50_rlocal_40'].obsSynth.keys()):
        print('possible dates', RUN['2016_50_rlocal_40'].obsSynth.keys())
    print('obs ' + var, RUN['2016_50_rlocal_40'].obsSynth[dd].data[var][cl])
    plot_part_cesar(RUN['2016_50_rlocal_40'].ensBg[dd],
          RUN['2016_50_rlocal_40'].ensAn[dd],
          RUN['2016_50_rlocal_40'].obsSynth[dd], cl, var, alpha = RUN['2016_50_rlocal_40'].alpha[dd])
for var in ['B4', 'B5']:
    snsscatter_2vars(RUN['2016_50_rlocal_40'].ensBg[dd].options,
                 RUN['2016_50_rlocal_40'].ensBg[dd],
                 RUN['2016_50_rlocal_40'].ensAn[dd],
                 RUN['2016_50_rlocal_40'].obsSynth[dd], 'SWE', var,cl )


# In[37]:


RUN[xp].ensBg['2017041410'].stack['SWE'].shape


# In[ ]:


plt.plot(RUN[xp].ensBg['2017041410'].stack['SWE'][124,:])
plt.plot(RUN[xp].ensAn['2017041410'].stack['SWE'][124,:])


# In[ ]:





# In[ ]:




