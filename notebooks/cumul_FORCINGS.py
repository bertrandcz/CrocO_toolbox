
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


from snowtools.utils.prosimu import prosimu
import os
from netCDF4 import Dataset
import numpy as np
import datetime
from scipy import stats
import matplotlib.pyplot as plt
import pickle
import pandas as pd
years = range(2013, 2017)
#years = [2013]
nens = 160
cl = 102
def pickleforcings(year, nens, listvars = ['IMPWET1', 'IMPDRY1', 'IMPWET2', 'IMPDRY2', 'PRECIP']):
    rootForcing = '/home/cluzetb/vortex/safran/12/'
    os.chdir('{0}forcing_{1}{2}B_31D_11_t1500_160'.format(
        rootForcing, year, year+1))
    forcEns = dict()
    for imb, mb in enumerate(map('mb{:04d}'.format, range(1, nens+1))):
        print(imb)
        forc = Dataset(
            '{0}/meteo/FORCING_{1}080106_{2}080106.nc'.format(mb, year, year+1))
        # print(list(forc.variables.keys()))
        for var in listvars:
            if var == 'PRECIP':
                gg = np.expand_dims(
                    forc.variables['Rainf'][:, cl] + forc.variables['Snowf'][:, cl], axis=1)
            else:
                gg = np.expand_dims(forc.variables[var][:, cl], axis=1)
            if imb == 0:
                forcEns[var] = gg
            else:
                forcEns[var] = np.concatenate((
                    forcEns[var], gg), axis=1)
        forc.close()
    with open('forcEns.pkl', 'wb') as f:
        pickle.dump(forcEns, f)


# In[32]:


years = [2013, 2016]
dictmembers = {2013: [66, 12, 122, 149],
               2016: [50, 153, 90, 117]}
listvars = ['IMPWET1', 'IMPDRY1', 'IMPWET2', 'IMPDRY2', 'PRECIP']
index = pd.MultiIndex.from_tuples([(year, mbsynth)
                                   for year in years for mbsynth in dictmembers[year]])
print(index)
dfforc = pd.DataFrame(np.empty((len(index), len(listvars))), columns = listvars, index = index)
import datetime

for year in years:
    

    f = '/home/cluzetb/vortex/safran/12/forcing_{0}{1}B_31D_11_t1500_160/forcEns.pkl'.format(year, year+1)
    if not os.path.exists(f):
        pickleforcings(year, nens)
    with open(f, 'rb') as f:
        forcEns = pickle.load(f)
    vecTime = np.array([datetime.datetime(year,8,1,6,0,0) + datetime.timedelta(hours=h) for h in range(forcEns[var].shape[0])])
    mask = [v >= datetime.datetime(year, 10,1,0,0,0) and v < datetime.datetime(year + 1,7,1,0,0,0) for v in vecTime]
    imask = [i for i in range(len(vecTime)) if mask[i] is True]
    vecTime = vecTime[mask]
    for var in dfforc.columns:
        forcEns[var] = forcEns[var][imask, :]
    fig = plt.figure(figsize = (5,5))
    gg = plt.plot(vecTime, np.cumsum(forcEns['PRECIP'], axis=0), color='0.7', alpha = 0.5)
    for mbsynth in dictmembers[year]: 
        for var in dfforc.columns:
            dfforc.loc[(year,mbsynth)][var] = stats.percentileofscore(np.sum(forcEns[var], axis = 0),
                                                               np.sum(forcEns[var][:,mbsynth - 1], axis = 0)) 

        plt.plot(vecTime, np.cumsum(forcEns['PRECIP'][:, mbsynth-1], axis=0), label = 'mb{0:04d}'.format(mbsynth))
        plt.ylabel('Cumul precip class 102')
        plt.xlabel('time')
        fig.autofmt_xdate()
    plt.legend()    

    


# In[34]:


dfforc.plot(figsize = (8,4),y = ['IMPWET1', 'IMPDRY1', 'IMPWET2', 'IMPDRY2','PRECIP'])


# In[ ]:




