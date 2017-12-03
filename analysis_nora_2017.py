
# coding: utf-8

# In[1]:


from common import *
import gzip
import io


# In[3]:


#print(datafiles)
datafiles = [os.path.join(DATA_DIR,'coolfiles',i) for i in os.listdir(DATA_DIR)]
nora_coolerfiles = [i for i in datafiles if '.cool' in i]


# In[4]:


coolertest = nora_coolerfiles[0]
print(coolertest)


# In[6]:


#decompress the file
#compressedFile = io.StringIO(response.read())
#decompressedFile = gzip.GzipFile(fileobj=coolertest)
#f = gzip.open(coolertest, 'rb')
cooler1 = cooler.Cooler(coolertest)


# In[7]:


cooler1

