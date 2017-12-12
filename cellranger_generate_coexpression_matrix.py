
# coding: utf-8

# In[ ]:


# This script will perform the covariation calculations on Orchestra. 
# This is meant to produce a pairwise contact file that we can utilize for analysis on Cooler


# In[ ]:


#Basics
import os, sys, re
import h5py
import multiprocessing as mp

#Packages for big data
import csv
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import scipy.sparse as sps
import scipy.spatial.distance as ds


#import cooler
#import re
#import seaborn


# In[ ]:


# matrix import function
def import_matrix_mtx_file(filename,aspd=False):
    inmat = scipy.io.mmread(source=filename)
    if (aspd==True):
        return(pd.SparseDataFrame(inmat))
    else:
        return(inmat)


# In[ ]:


#####
# Main

# Set imports; these will be the only three inputs we need
in_matrix_mtx = sys.argv[1] # os.path.join(bcell_matrix_dir, "matrix.mtx")
in_genes = sys.argv[2] # os.path.join(bcell_matrix_dir, "genes.tsv")
in_barcodes = sys.argv[3] # os.path.join(bcell_matrix_dir, "barcodes.tsv")


# In[ ]:


# Import the files
# 1. inmatrix
mp_pool = mp.Pool(5)
infile = mp_pool.map(import_matrix_mtx_file, [inmatrix])
mp_pool.close()

# 2. genes_df
gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]
genes_df = pd.DataFrame()
genes_df = genes_df.assign(genename=pd.Series(gene_names),geneid=pd.Series(gene_ids))

# 3. barcodes_df
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]

