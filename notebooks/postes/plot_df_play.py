
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib inline')
from postes.explore_obs import find_common_massif_obsbase
import matplotlib.pyplot as plt
import seaborn as sns

plt.close('all')
obsPath = '/home/cluzetb/assim_postes/obs_postes_1983_2019.csv'
obsPath = '/home/cluzetb/assim_postes/OBS_2010080100_2020021023_types1346789.csv'

df, listp, _, _ = find_common_massif_obsbase(obsPath, massifs=[8,9,12,13,15,16])
# D. finalize (rename cols etc.)
df = df.rename(columns = {int(p.find('number').text): p.find('name').text for p in listp})
print('cols', df.columns)
print('number of observed posts:', len(df.columns))


#  E. Play.
corr = df.corr()


plt.figure(figsize = (7,15))
sns.heatmap(corr, vmin = -1, vmax = 1, center = 0, cmap = 'RdBu', xticklabels=True, yticklabels=True)


for col in df.columns:
    plt.figure()
    plt.scatter(df.index, df[col],)
    plt.title(col)
    plt.gcf().autofmt_xdate()
    plt.show()
    
    


# In[2]:

df.columns


# In[3]:

df


# In[ ]:



