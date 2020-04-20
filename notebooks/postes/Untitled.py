
# coding: utf-8

# In[112]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
import datetime
import numpy as np
from utilpp import RMSE,spread
from CramponPp import CramponPp
from crampon import set_options
import matplotlib.pyplot as plt
from postes.explore_metadata import find_name_station
from utilcrampon import Pgd

# ########## PARAMS ##########################
run = 'global'  # 'ol' for openloop
nens = 40
neff = 7
assimvars = 'DEP'
year = 2016
assimilate_every = 7
vconf = 'postes_8_9_12_13_15_16_csv'
selMassif = [12]
##############################################

# preparation
sensorBase = '_'.join(list(map(str, selMassif)))
pgdPath = '/cnrm/cen/users/NO_SAVE/cluzetb/vortex/s2m/' + vconf + '/spinup/pgd/PGD_' + vconf + '.nc'
pgd = Pgd(pgdPath)

# generate runs on beaufix
classesId = np.array([i for i, m in enumerate(pgd.massif) if m in selMassif])
# sort classes per altitude:
classesId = classesId[np.argsort(pgd.elev[classesId])]

xps = {str(stat): sensorBase + '_X' + str(stat) for stat in pgd.station[classesId]}
xps[str(sensorBase)] = str(sensorBase)
print(xps)

runs = dict()
for station in list(xps.keys()):
    print('sensor', xps[station])
    args = [
        '/home/cluzetb/snowtools_git/assim/crampon.py',
        '--xpid', '{0}_{1}_{2}_{3}_{4}_{5}'.format(year, run, nens, xps[station], neff, assimilate_every),
        '--xpidol', '{0}_ol_40'.format(year),
        '--sensor', xps[station],
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

    runs[station] = CramponPp(options, conf)


# In[113]:


# plot scores only for the considered massif (12)
from scores.ensemble import EnsembleScores
from utilpp import RMSE
CRPSS = np.empty((len(runs.keys()), len(runs.keys())-1))
ReliS = np.empty((len(runs.keys()), len(runs.keys())-1))
for ii, station in enumerate(list(map(str,pgd.station[classesId])) + list(map(str,selMassif))):
    print(find_name_station(str(station)))
    run = runs[station]
    mobs = [i for i, t in enumerate(run.obsTs['time']) 
            if (t.hour ==6 and ((t.year == year and t.month>9)  or (t.month <7 and t.year ==year+1)))]
    mensAn = [i for i, t in enumerate(run.ensProAn['time']) if (t.hour ==6 and (t.month>9 or t.month <7))]
    mensOl = [i for i, t in enumerate(run.ensProOl['time']) if (t.hour ==6 and (t.month>9 or t.month <7))]
    obs = np.ma.masked_invalid(run.obsTs['DEP'][mobs, :][:,classesId])
    ensAn = run.ensProAn['DEP'][mensAn, :,:][:,classesId,:]
    ensOl = run.ensProOl['DEP'][mensOl, :,:][:,classesId,:]
    CRPSAn = np.array([EnsembleScores(list(range(ensAn.shape[0])),
                                                 list(range(ensAn.shape[0])),
                                                 obs[:,cl],
                                                 ensAn[:,cl,:].T,
                                                ).CRPS_decomp() for cl in range(len(classesId))])
    CRPSOl = np.array([EnsembleScores(list(range(ensOl.shape[0])),
                                                 list(range(ensOl.shape[0])),
                                                 obs[:,cl],
                                                 ensOl[:,cl,:].T,
                                                ).CRPS_decomp() for cl in range(len(classesId))])
    CRPSS[ii,:] = np.ma.masked_invalid(np.array([1 - CRPSAn[i,0]/CRPSOl[i,0] if CRPSOl[i,0]>0. else np.nan 
                                           for i in range(CRPSAn.shape[0])]))
    ReliS[ii,:] = np.ma.masked_invalid(np.array([1 - CRPSAn[i,1]/CRPSOl[i,1] if CRPSOl[i,0]>0. else np.nan 
                                           for i in range(CRPSAn.shape[0])]))
    '''
    CRPSAn = RMSE(ensAn, obs, aggrDomain=False, aggrTime = True)
    CRPSOl = RMSE(ensOl, obs, aggrDomain=False, aggrTime = True)

    CRPSS[ii,:] = np.ma.masked_invalid(np.array([1 - CRPSAn[i]/CRPSOl[i] if CRPSOl[i]>0. else np.nan 
                                           for i in range(CRPSAn.shape[0])]))
    ReliS[ii,:] = np.ma.masked_invalid(np.array([1 - CRPSAn[i]/CRPSOl[i] if CRPSOl[i]>0. else np.nan 
                                           for i in range(CRPSAn.shape[0])]))
    '''


# In[114]:


get_ipython().run_line_magic('matplotlib', 'inline')
from assim.postes.explore_metadata import find_name_station
from math import trunc
fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
plt.subplots_adjust(wspace = 0.01, top = 0.95, left = 0.35, bottom = 0.4, right = 0.9)

axs[0].imshow(np.hstack((CRPSS,np.expand_dims(np.nanmean(CRPSS, axis = 1), axis = 1))),
             cmap = 'RdBu',
             vmin = -1, vmax = 1)
im = axs[1].imshow(np.hstack((ReliS,np.expand_dims(np.nanmean(ReliS, axis = 1), axis = 1))),
             cmap = 'RdBu',
             vmin = -1, vmax = 1)
axs[0].set_yticks(range(0,CRPSS.shape[0]+1))
axs[0].set_yticklabels([find_name_station(pgd.station[c]) + ', ({0:d}m)'.format(trunc(pgd.elev[c])) for c in classesId] + ['assimilating all stations'])
axs[1].set_yticklabels([])
for i in [0,1]:
    axs[i].set_xticks(range(0,CRPSS.shape[0]))
    axs[i].set_xticklabels([find_name_station(s) for s in pgd.station[classesId]] + ['mean'])
    axs[i].tick_params(axis='x', rotation=90)
    axs[i].plot([-0.5, CRPSS.shape[0]-0.5], [CRPSS.shape[0]-1.5,CRPSS.shape[0] -1.5], color = 'k', ls = '--', lw = 2)
axs[0].set_title('CRPSS')
axs[0].annotate('removing', xy=(-1.05,0.55), xytext=(-1.1,0.55), xycoords='axes fraction', 
            fontsize=8, ha='center', va='center',
            bbox=dict(boxstyle='square', fc='white'),
            arrowprops=dict(arrowstyle='-[, widthB=8., lengthB=1.5', lw=2.0))
axs[1].set_title('ReliS')
# add colorbar
cax = fig.add_axes([0.91, 0.405, 0.02, 0.55])
cb = plt.colorbar(im, cax=cax, )


plt.savefig('/home/cluzetb/notebooks/postes/CRPS_12X_{0}.pdf'.format(year), dpi = 200)


# In[136]:


# spaghetti in all classes
get_ipython().run_line_magic('matplotlib', 'inline')
for ii, station in enumerate(list(map(str,pgd.station[classesId]))):
    runX = runs[station]
    run12 = runs['12']
    pt = classesId[ii]
    fig, axs = plt.subplots(ncols = 2, figsize = (10,5), dpi = 100)
    for ax, run in zip(axs, [runX, run12]):
        ax.plot(run.ensProOl['time'],run.ensProOl['DEP'][:,pt,:], color='0.7', alpha = .4) 
        ax.plot(run.ensProAn['time'],run.ensProAn['DEP'][:,pt,:], color='g', alpha = 0.4) 
        
        ax.plot(run.ensProAn['time'], np.median(run.ensProAn['DEP'][:,pt,:], axis = 1),
                color='k',lw = 2)
        ax.plot(run.ensProOl['time'], np.median(run.ensProOl['DEP'][:,pt,:], axis = 1),
                color='0.2',lw = 2)
        ax.plot(run.obsTs['time'],run.obsTs['DEP'][:,pt], color = 'r', ls = '--')
        for dd in run.options.dates: 
            #ax.scatter(datetime.datetime.strptime(dd, '%Y%m%d%H'), run.obsReal[dd].data['DEP'][pt],
            #           color = 'r', zorder = 100.) 
            ax.set_title(find_name_station(run.pgd.station[pt])) 
        ax.set_xlim([run.ensProAn['time'][0], run.ensProAn['time'][-1]])
        ax.set_ylim([0,3])
        fig.autofmt_xdate()
    plt.show()


# In[53]:


'''# study on the location of posts
plt.figure(figsize =(5,5))
for i, iid in enumerate(classesId):
    plt.scatter(pgd.elev[classesId],CRPSS[i,:],)
    st = pgd.station[classesId[i]]
    if not np.isnan(CRPSS[-5,i]):
        plt.annotate(find_name_station(str(st)), xy = (1200, CRPSS[i, -3]), ha = 'right')
    plt.xlim([900, 3000])
plt.show()
'''


# In[131]:


for xp in runs.keys():
    fig = plt.figure()
    plt.plot([datetime.datetime.strptime(dd,'%Y%m%d%H') for dd in runs[xp].alpha.keys()],[runs[xp].alpha[dd] for dd in runs[xp].alpha.keys()], label = xp)
    plt.legend()
    fig.autofmt_xdate()


# In[ ]:




