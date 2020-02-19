
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib.pyplot as plt
from utilpp import read_BG
from plotcrocO import Pie
from SemiDistributed import PrepAbs
from copy import copy
from Operators import PrepEnsOperator
from utilcrocO import setSubsetclasses,cm2inch


# In[2]:


get_ipython().run_line_magic('run', "/home/cluzetb/assim/scripts/load_multi_xp.py 'plot_bg'")


# In[3]:


pt=77
dates = RUN[list(RUN.keys())[0]].options.dates
dates = ['2015022010']


# In[4]:


'''
for date in dates:
    print(date)
    for xp in RUN.keys():
        bg = read_BG(RUN[xp].options)
        focusCl = setSubsetclasses(RUN[xp].pgd, RUN[xp].options.classesE,
                                      RUN[xp].options.classesA,RUN[xp].options.classesS)[0]
        focusColors = [0] * len(focusCl)
        focusCl.append(pt)
        focusColors.append(-0.6)
        opts = copy(RUN[xp].options)
        opts.vars = RUN[xp].options.vars
        opts.ppvars = RUN[xp].options.vars
        opts.classesS = ['0', '20', '40']
        bgData = PrepAbs(date, opts, ptinom = 'bg')
        bgData.data = dict()
        for var in RUN[xp].options.vars:
            if var != 'SWE':
                bgData.data[var] = bg[date][var][pt,:]
        bgData.isloaded = True
        pieData = Pie(bgData, focusCl = focusCl, focusColors = focusColors)
        pieData.plot(cmap = 'RdBu', clims = [-1,1], title = '', savefig = True, colortitle = '')
'''


# In[5]:


for date in dates:
    print(date)
    #f = plt.figure(figsize = cm2inch(8.3,15))
    #gs0 = gridspec.GridSpec(3, 1, figure=f)
    for xp in RUN.keys():
        focusCl = setSubsetclasses(RUN[xp].pgd, RUN[xp].options.classesE,
                                      RUN[xp].options.classesA,['0', '20'])[0]
        focusColors = [0] * len(focusCl)
        focusCl.append(pt)
        focusColors.append(-0.6)
        RUN[xp].options.ppvars = RUN[xp].options.vars
        RUN[xp].options.classesS = ['0','20','40']
        po = PrepEnsOperator(RUN[xp].ensOl[date])
        bgData = po.bgcovariance(clA=pt)
        pieData = Pie(bgData, focusCl = focusCl, focusColors = focusColors, maketitle = True)
        # if 'DEP' in xp:
        #     sg = gs0[2]
        # else:
        #     sg = gs0[0:1]
        #    print(sg)
        pieData.plot(cmap = 'RdBu', clims = [-1,1], title = '', savefig = True, colortitle = '$R$')
    plt.show()


# In[6]:


pwd


# In[ ]:




