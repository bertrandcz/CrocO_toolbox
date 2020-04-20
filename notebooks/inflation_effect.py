
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
from plotcrampon import prepare_spaghetti_Marie
from utilcrampon import cm2inch,setSubsetclasses


# In[2]:


get_ipython().run_line_magic('run', "/home/cluzetb/assim/scripts/load_multi_xp.py 'inflation_effect'")


# In[3]:


focusCl = setSubsetclasses( RUN[xp].pgd, options.classesE, options.classesA, options.classesS)[0]
nfocusCl = [p for p in range(RUN[xp].pgd.npts) if p not in focusCl ]


# In[4]:


def colors_legs(xp):
    if xp.endswith('_10'):
        col = 'm'
        lab = '$N_{eff}^* = 10$'
    elif xp.endswith('_1'):
        col = 'darkorange'
        lab = '$N_{eff}^* = 1$'
    else:
        col = 'royalblue'
        lab = '$N_{eff}^* = 7$'
    return col, lab


# In[6]:


fig, axs = plt.subplots(nrows=3, figsize = cm2inch(8.3, 15))
fig.subplots_adjust(left=0.2,right = 0.99, top = 0.99, bottom = 0.01)
legEl = {}
for i,xp in enumerate(RUN.keys()):
    col, lab = colors_legs(xp)
    print('xp : ', xp)
    legEl = prepare_spaghetti_Marie(RUN[xp], 'SWE', 77, ax = axs[0], color =col, label=lab, density = True, title = False, legEl = legEl)
    legEl = prepare_spaghetti_Marie(RUN[xp], 'SWE', 77, ax = axs[1], color =col, label=lab, score = '$\sigma$', title = False, legEl = legEl)
    legEl = prepare_spaghetti_Marie(RUN[xp], 'SWE', 77, ax = axs[2], color =col, label=lab, score = 'AE', title = False, legEl = legEl)

legName = []
legItem = []
for key, val in legEl.items():
    legItem.append(val)
    legName.append(key)
fig.legend(legItem, legName, bbox_to_anchor=[0.5, 0.07], loc = 'center',ncol = 2)
for i, label in enumerate(['(a)', '(b)','(c)']):
    ax = axs[i]
    ax.text(0.04, 0.95, label, transform=ax.transAxes,
      fontsize=14, fontweight='bold', va='top')


#toggle comments for saving.
fig.savefig('/home/cluzetb/notebooks/articleGMD/fig/res/infl.pdf')
#plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




