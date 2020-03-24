
# coding: utf-8

# In[7]:


get_ipython().run_line_magic('matplotlib', 'notebook')
from osgeo import ogr, osr


# In[14]:


import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
from pyproj import Proj, transform
inProj = Proj(init = 'epsg:2154')
outProj = Proj(init = 'epsg:4326')
sf = shp.Reader("/home/cluzetb/snowtools_git/DATA/massifs_alpes_L93E.shp")

plt.figure()
for shape in sf.shapeRecords():
    x = [i[0] for i in shape.shape.points[:]]
    y = [i[1] for i in shape.shape.points[:]]
    x2,y2 = transform(inProj,outProj,x,y)
    plt.plot(x2,y2)
plt.show()


# In[ ]:




