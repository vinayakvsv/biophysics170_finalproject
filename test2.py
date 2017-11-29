
# coding: utf-8

# In[1]:


print("hello world")


# In[13]:


#See https://github.com/mirnylab/cooler-binder/blob/master/cooler_api.ipynb for tutorial
import matplotlib.pyplot as plt
import numpy as np
import pandas
import h5py
import sys,os,re

import cooler


# In[3]:


get_ipython().magic('matplotlib inline')


# In[ ]:


get_ipython().system('curl -O ftp://cooler.csail.mit.edu/coolers/hg19/Rao2014-GM12878-MboI-allreps-filtered.5kb.cool')
get_ipython().system('mv Rao2014-GM12878-MboI-allreps-filtered.5kb.cool ./hic_data/')


# In[5]:


filepath = './hic_data/Rao2014-GM12878-MboI-allreps-filtered.5kb.cool'


# In[6]:


c = cooler.Cooler(filepath)
print(c.info)
print(c.chroms())


# In[14]:


mat = c.matrix(balance=False, sparse=True)[1000:1200, 1000:1200]
mat
arr = mat.toarray()
arr
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
im = ax.matshow(np.log10(arr), cmap='YlOrRd')
fig.colorbar(im)

