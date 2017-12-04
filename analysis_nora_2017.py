
# coding: utf-8

# In[1]:


from common import *
import gzip
import io
import matplotlib.pyplot as plt


# In[61]:


#print(datafiles)
cooldir = os.path.join(DATA_DIR,'coolfiles')
datafiles = [os.path.join(cooldir,i) for i in os.listdir(cooldir) if ".cool" in i]
print(datafiles)
nora_coolerfiles = [i for i in datafiles if '.cool' in i]


# In[64]:


#get metadata
#import re
filenames = [i for i in os.listdir(cooldir) if ".cool" in i]
files_meta = pandas.DataFrame([re.split('_|-|[.]|.cool',i)[:-1] + [j] for i,j in zip(filenames,datafiles)])
files_meta.columns = ["sample","treatment","read number","baseres","file"]
print(files_meta)


# In[110]:


#let's create cooler files for all objects
cooler1 = cooler.Cooler(datafiles[0])#[cooler.Cooler(i) for i in files_meta.loc[:,"file"]]
coolers = [cooler.Cooler(i) for i in files_meta.loc[:,"file"]]

#higher resolution files
higrescool_meta = files_meta[files_meta.baseres == '20000']
print(higrescool_meta)

higrescool = [cooler.Cooler(i) for i in higrescool_meta.loc[:,"file"]]
print(higrescool)


# In[119]:


for i in higrescool:
    imat = i.matrix(balance=False, sparse=True)[500:1200, 500:1200]
    arr = imat.toarray()
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    im = ax.matshow(np.log10(arr), cmap='YlOrRd')
    fig.colorbar(im)
    plt.show()


# # Call TADs in this dataset

# In[ ]:


# Questions for Leonid
1. Nora data: 

