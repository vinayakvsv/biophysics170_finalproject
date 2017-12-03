
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

# In[9]:


# imports
from common import * #see common.py for the packages that were just imported

#contents of this directory
datafiles = [os.path.join(DATA_DIR,i) for i in os.listdir(DATA_DIR)]


# In[13]:


#to test, let's draw the cooler files. 
#How to work with compressed .cool files: https://stackoverflow.com/questions/15352668/download-and-decompress-gzipped-file-in-memory
#c = cooler.Cooler("/Users/vinayakvsv/biophysics170/biophysics170_finalproject/../data/GSM2644945_Untreated-R1.100000.cool.gz")
mm9genes = cooler.read_chromsizes(filepath_or="./mm-data/mm9_refseq.txt")

