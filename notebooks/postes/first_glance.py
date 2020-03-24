
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
import datetime
import numpy as np
from utilpp import RMSE,spread
from CrocOpp import CrocOpp
from crocO import set_options
import matplotlib.pyplot as plt
from postes.explore_metadata import find_name_station

# ########## PARAMS ##########################
run = 'global'  # 'ol' for openloop
nens = 40
neff = 7
assimvars = 'DEP'
year = 2014
assimilate_every = 7
vconf = 'postes_8_9_12_13_15_16_csv'
selMassif = [12]

runs = dict()
for sensor in ['12', '12_X38191400']:
    args = [
        '/home/cluzetb/snowtools_git/assim/crocO.py',
        '--xpid', '{0}_{1}_{2}_{3}_{4}_{5}'.format(year, run, nens, sensor, neff, assimilate_every),
        '--xpidol', '2014_ol_40',
        '--sensor', sensor,
        '--vconf', 'postes_8_9_12_13_15_16_csv',
        '-d', 'all',
        '--vars', 'DEP',
        '--ppvars', 'DEP',
        '-o', 'beaufix',
        # '--readprep',
        '--readobs',
        '--readaux',
        '--readoper',
    ]
    options, conf = set_options(args)

    runs[sensor] = CrocOpp(options, conf)


# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
from pyproj import Proj, transform
plt.figure(figsize = (10,10))
for i in [0,1]:
    plt.subplot(121 + i)
    run = runs['12'] if i==0 else runs['12_X38191400']

    mobs = [i for i, t in enumerate(run.obsTs['time']) 
            if (t.hour ==6 and ((t.year == 2013 and t.month>9)  or (t.month <7 and t.year ==2015)))]
    mens = [i for i, t in enumerate(run.ensProAn['time']) if (t.hour ==6 and (t.month>9 or t.month <7))]
    moper = [i for i, t in enumerate(run.oper['time']) if (t.hour ==6 and (t.month>9 or t.month <7))]
    obs = np.ma.masked_invalid(run.obsTs['DEP'][mobs, :])
    ensAn = run.ensProAn['DEP'][mens, :,:]
    ensOl = run.ensProOl['DEP'][mens, :,:]
    oper = run.oper['DEP'][moper,:,:]
    rmseAn = RMSE(ensAn, obs, aggrTime = True, aggrDomain = False)
    rmseOl = RMSE(ensOl, obs, aggrTime = True, aggrDomain = False)
    rmseOper = RMSE(oper, obs, aggrTime = True, aggrDomain = False)
    spAn = spread(ensAn, aggrTime = True, aggrDomain = False)
    spOl = spread(ensOl, aggrTime = True, aggrDomain = False)
    # CRPS 
    CRPSAn = np.array([EnsembleScores(list(range(ensAn.shape[0])),
                                                 list(range(ensAn.shape[0])),
                                                 obs[:,cl],
                                                 ensAn[:,cl,:].T,
                                                ).CRPS_decomp() for cl in range(run.pgd.npts)])
    CRPSOl = np.array([EnsembleScores(list(range(ensOl.shape[0])),
                                                 list(range(ensOl.shape[0])),
                                                 obs[:,cl],
                                                 ensOl[:,cl,:].T,
                                                ).CRPS_decomp() for cl in range(run.pgd.npts)])
    CRPSS = np.ma.masked_invalid(np.array([1 - CRPSAn[i,1]/CRPSOl[i,1] if CRPSOl[i,1]>0. else np.nan 
                                           for i in range(CRPSAn.shape[0])]))
    print(CRPSS)
    plt.scatter(run.pgd.lon, run.pgd.lat, c = CRPSS, cmap = 'RdBu', vmin = -1, vmax = 1, zorder = 100)
    #plt.colorbar()
    for i, st in enumerate(run.pgd.station):
        if (not(np.ma.is_masked(rmseOl[i])) and not(np.ma.is_masked(rmseAn[i]))):
            pass
            # plt.gca().annotate(st, xy = (run.pgd.lon[i], run.pgd.lat[i]), zorder = 1, fontsize = 8)
    inProj = Proj(init = 'epsg:2154')
    outProj = Proj(init = 'epsg:4326')
    sf = shp.Reader("/home/cluzetb/snowtools_git/DATA/massifs_alpes_L93E.shp")

    for shape in sf.shapeRecords():
        x = [i[0] for i in shape.shape.points[:]]
        y = [i[1] for i in shape.shape.points[:]]
        x2,y2 = transform(inProj,outProj,x,y)
        plt.plot(x2,y2, color = 'k', zorder = 0)
    plt.axis('equal')
    plt.xlim([5.5, 7])
    plt.ylim([44.7, 45.7])
#plt.colorbar()

plt.show()


# In[ ]:





# In[ ]:


"""
%matplotlib notebook
from utilpp import set_itimes
from scores.ensemble import EnsembleScores
mobs = [i for i, t in enumerate(run.obsTs['time']) 
        if (t.hour ==6 and ((t.year == 2013 and t.month>9)  or (t.month <7 and t.year ==2015)))]
mens = [i for i, t in enumerate(run.ensProAn['time']) if (t.hour ==0 and (t.month>9 or t.month <7))]
obs = np.ma.masked_invalid(run.obsTs['DEP'][mobs, :])
ensAn = run.ensProAn['DEP'][mens, :,:]
ensOl = run.ensProOl['DEP'][mens, :,:]

# compute scores
rmseAn = RMSE(ensAn, obs, aggrTime = True, aggrDomain = False)
rmseOl = RMSE(ensOl, obs, aggrTime = True, aggrDomain = False)
spAn = spread(ensAn, aggrTime = True, aggrDomain = False)
spOl = spread(ensOl, aggrTime = True, aggrDomain = False)
""""""


# In[ ]:


"""
# plot
plt.figure()
plt.plot([0,1], [0,1], color = 'k')
plt.scatter(rmseOl,rmseAn)
for i, st in enumerate(run.pgd.station):
    if (not(np.ma.is_masked(rmseOl[i])) and not(np.ma.is_masked(rmseAn[i]))):
        plt.annotate(st, xy = (rmseOl[i], rmseAn[i]))
plt.axis('equal')
plt.xlim([0,1])
plt.ylim([0,1])
plt.figure()
plt.plot(rmseOl, color = 'k')
plt.plot(rmseAn, color = 'b')
plt.figure()
plt.plot(spOl, color = 'k', ls = '--')
plt.plot(spAn, color = 'b', ls = '--')
"""


# In[ ]:





# In[ ]:





# In[5]:


# spaghetti in all classes
get_ipython().run_line_magic('matplotlib', 'inline')
run = runs["12"]
for pt in range(len(run.pgd.station)):
#for pt in range(1):
    plt.figure() 
    plt.plot(run.ensProOl['time'],run.ensProOl['DEP'][:,pt,:], color='0.7', alpha = 0.4) 
    plt.plot(run.ensProAn['time'],run.ensProAn['DEP'][:,pt,:], color='g', alpha = 0.4) 
    plt.plot(run.ensProAn['time'], np.median(run.ensProOl['DEP'][:,pt,:], axis = 1),
             color='0.2',lw = 2)    
    plt.plot(run.ensProAn['time'], np.median(run.ensProAn['DEP'][:,pt,:], axis = 1),
             color='k',lw = 2)
    plt.plot(run.oper['time'], run.oper['DEP'][:,pt,0],
             color='b',lw = 2)
    plt.plot(run.obsTs['time'],run.obsTs['DEP'][:,pt], color = 'r', ls = '--')
    for dd in run.options.dates: 
        plt.scatter(datetime.datetime.strptime(dd, '%Y%m%d%H'), run.obsReal[dd].data['DEP'][pt], color = 'r', zorder = 100.) 
        plt.title(find_name_station(run.pgd.station[pt])) 
    plt.gca().set_xlim([run.ensProAn['time'][0], run.ensProAn['time'][-1]])
    plt.gcf().autofmt_xdate()
    plt.show()


# In[ ]:




