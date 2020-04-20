
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from utils.prosimu import prosimu
import glob
pgd = netCDF4.Dataset('/cnrm/cen/users/NO_SAVE/cluzetb/vortex/s2m/12/2016_50_global/crampon/compar2to3/PGD.nc')
slope = pgd.variables['SSO_SLOPE'][:]


# In[22]:


for prep in sorted(glob.glob('/cnrm/cen/users/NO_SAVE/cluzetb/vortex/s2m/12/2016_50_global/mb0001/bg/PREP_2017041410.nc')):
    dd= prosimu(prep)
    
    sum_tot =  np.nansum([dd.read('WSN_VEG' + str(i + 1)) / dd.read('RSN_VEG' + str(i + 1)) /                             np.cos(np.arctan(slope)) for i in range(0, 50)],axis= 0) 
    gg = netCDF4.Dataset(prep,'r') 
    ggpro = netCDF4.Dataset('/cnrm/cen/users/NO_SAVE/cluzetb/vortex/s2m/12/2016_50_global/mb0001/pro/PRO_2017041410_2017042010.nc', 'r')
    
    dep_tot = gg.variables['DEP_TOT'][:]
    dep_tot_pro = ggpro.variables['DSN_T_ISBA'][:][0,:]
    ttime = ggpro.variables['time'][:]
    diff = sum_tot - dep_tot
    plt.figure()
    plt.scatter(slope,np.squeeze(diff))
    plt.figure()
    plt.plot(np.squeeze(diff))
    plt.plot(np.squeeze(sum_tot-dep_tot_pro))
    plt.title(prep[-13:-3])
    plt.show()


# In[23]:


units = ggpro.variables['time'].getncattr('units')


# In[24]:


gg3 =netCDF4.num2date(ttime, units =units)


# In[25]:


gg3


# In[ ]:




