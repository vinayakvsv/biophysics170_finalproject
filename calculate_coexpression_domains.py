
# coding: utf-8

# # Introduction
# 
# In this script, we will calculate co-expression domains from the RNA-seq data in Nora et al Cell 2017. 
# 
# The procedure will go as follows:
# 
# 1. Import RNA-seq bigwig files, and define signal over genes (mm9 Refseq annotation), and store as Cooler files
#     a. visualize in HiGlass with TAD data
# 2. Compute a co-expression matrix of genes.
# 3. Utilize TAD-calling strategies to define blocks of significantly-correlated blocks.
# 4. Export the result as a 

# In[ ]:


# imports
from common import * #see common.py for the packages that were just imported

#contents of this directory
datafiles = [os.path.join(DATA_DIR,i) for i in os.listdir(DATA_DIR)]
rnaseqbw = [i for i in datafiles if "RNA-seq" in i]
print(rnaseqbw)


# In[ ]:


#let's import the genes
mm9genes = pandas.read_table(filepath_or_buffer="./mm-data/mm9.txt")
print(mm9genes)


# In[ ]:


#import pyBigWig


# In[3]:


import pybedtools

