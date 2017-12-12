
# coding: utf-8

# # Introduction
# This notebook loads in analysis of B-cells isolated from PBMC's from a healthy donor and analyzed on the 10X Genomics platform. Data was downloaded from the 10X website (https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/fresh_68k_pbmc_donor_a). Instructions on using the ``matrix.mtx`` files are included here: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices. We will follow this tutorial

# In[1]:


# import necessary packages
import csv
import os
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import numpy as np
import scipy.sparse as sps
import multiprocessing as mp
import h5py
import cooler
import re
import seaborn
import scipy.spatial.distance as ds


# In[2]:


# set paths
genome = "hg19"
celltype = "bcells"
matfolder = "filtered_matrices_mex"
clusterfolder = "analysis_csv"
scrnaseq_dir = "../scrnaseq_10x/"

bcell_matrix_dir = os.path.join(scrnaseq_dir,celltype,matfolder,genome)
bcell_analysis_dir = os.path.join(scrnaseq_dir,celltype,clusterfolder)

pca_dir = os.path.join(bcell_analysis_dir,"pca")
tsne_dir = os.path.join(bcell_analysis_dir,"tsne")
kmeans_dir = os.path.join(bcell_analysis_dir,"kmeans")


# In[3]:


#set import functions

def import_matrix_mtx_file(filename,aspd=False):
    inmat = scipy.io.mmread(source=filename)
    if (aspd==True):
        return(pd.SparseDataFrame(inmat))
    else:
        return(inmat)


# In[4]:


#import the file. Use multiprocessing to make this run smoother

inmatrix = os.path.join(bcell_matrix_dir, "matrix.mtx")
genes_path = os.path.join(bcell_matrix_dir, "genes.tsv")
barcodes_path = os.path.join(bcell_matrix_dir, "barcodes.tsv")


# In[5]:


#impor the matrix file
mp_pool = mp.Pool(5)
infile = mp_pool.map(import_matrix_mtx_file, [inmatrix])
mp_pool.close()


# In[6]:


mat = infile[0]
print(mat.shape)
print(type(mat))
print(mat)
#rows are genes, columns are cells!


# ```
# #attempts at h5py
# 
# #attempt to save as an h5py file
# 
# #mat.to_hdf(os.path.join(bcell_matrix_dir, "matrix.h5"), 'sparse_df')
# 
# # with h5py.File(os.path.join(bcell_matrix_dir, "matrix.h5"), 'w') as hf:
# #     hf.create_dataset("bcell_8k",  data=mat)
# 
# #filepath = os.path.join(bcell_matrix_dir, "matrix.h5")
# # filepath = os.path.join(bcell_matrix_dir, "matrix.h5")
# # print(filepath)
# # store1 = pd.HDFStore(os.path.join(bcell_matrix_dir, "matrix.h5"))
# 
# 
# # import h5py
# # h5 = h5py.File(filepath="../scrnaseq_10x/bcells/", 'r')
# 
# # import cooler
# # cooler.matrix()
# 
# 
# # #attempt to load
# # with h5py.File(os.path.join(bcell_matrix_dir, "matrix.h5"), 'r') as hf:
# #     indata = hf['bcell_8k'][:]
# ```

# # Tiling the transcriptome
# 
# Since we have over 32,000 genes, we would like to condense them into genomic tiles that would allow us to work with a smaller, sparser matrix. 
# 
# 1. We will import the GTF, intersect it with 100-kb genomic tiles (giving us roughly 10,000 regions), and convert the matrix from a ~32,000 x 10,000 matrix to a ~10,000 x 10,000 matrix. 
# 2. The reads for genes in each bin for each cell will be summed up so that the "transcriptional output" of the bin can be compared with other bins. 
# 3. We will then calculate the correlation matrix for the cell-bin matrix such that we correlate the transcriptional outputs of bin 1 with bin 2. We then save this as a "contact map" of sorts and output it to a file
# 4. We then convert the matrix into a Cooler file (using the 100-kb bins from earlier) and work it up in Ccooler or Hi-Glass

# In[7]:


# load the gene IDs and names
gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]


# In[8]:


# load the barcode IDs
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]


# In[9]:


genes_df = pd.DataFrame()
genes_df = genes_df.assign(genename=pd.Series(gene_names),geneid=pd.Series(gene_ids))


# for multi-core jobs
# https://stackoverflow.com/questions/33480675/multiprocessing-on-a-set-number-of-cores

# ## Import ENSEMBL GTF

# In[10]:


#set ensembl file
#scrnaseq_dir = "../scrnaseq_10x/"
hg19gtf = "../hg19ref/refdata-cellranger-hg19-1.2.0/genes/genes.gtf"
hg19gtf_dir = os.path.join(scrnaseq_dir,hg19gtf)

hg19gtf_genes = pd.read_csv(hg19gtf,delimiter="\t",comment="#",header=None)


# In[11]:


hg19gtf_genes.columns = ["chr","database","type","start","end","other","strand","other2","metadata"]


# In[12]:


metadata = hg19gtf_genes.metadata.apply(lambda x: re.split(r"[a-zA-Z]+_[a-zA-Z]+ \"|\";",x)[1])


# In[13]:


hg19gtf_genes['gene_id'] = pd.Series(metadata)
#hg19gtf_genes
hg19gtf_genes_expression = hg19gtf_genes[hg19gtf_genes.gene_id.isin(gene_ids) & hg19gtf_genes.type.isin(["gene"])]
print(hg19gtf_genes_expression.shape)


# # Visualization
# Let's view the cell clusters

# We now load T-SNE metadata for this population. We will show what the populations look like and then isolate the cells responsible to generate co-expression blocks

# In[14]:


#tsne_data = scipy.io.mmread(os.path.join(tsne_dir, "projection.csv")) 
tsne_data = os.path.join(tsne_dir, "projection.csv")
tsne_points = [row for row in csv.reader(open(tsne_data), delimiter=",")]
tsne_meta = pd.DataFrame(tsne_points[1:])
tsne_meta.columns = tsne_points[0]


# In[15]:


#import clustering
cluster_10 = os.path.join(kmeans_dir,"10_clusters")
print(cluster_10)
#get the cluster designations
cluster_10_ids = pd.read_csv(os.path.join(cluster_10,"clusters.csv"))
#cluster_10_ids
#kmeans_10 = os.path.join(pca_dir,"10_clusters") 


# In[16]:


plt.figure(figsize=(10, 10))
plt.scatter(tsne_meta.iloc[:,1],tsne_meta.iloc[:,2],c = cluster_10_ids.iloc[:,1])

xlist = tsne_meta.iloc[:,1].tolist()
ylist = tsne_meta.iloc[:,2].tolist()

plt.xticks([float(min(xlist))-1, float(max(xlist))+1])
plt.yticks([float(min(ylist))-1, float(max(ylist))+1])

# plt.xticks(tsne_meta.iloc[:,1].tolist())
# plt.yticks(tsne_meta.iloc[:,2].tolist())

# plt.xticks(np.arange(float(min(xlist))-1, float(max(xlist))+1,0.01))
# plt.yticks(np.arange(float(min(ylist))-1, float(min(ylist))+1,0.01))

plt.show()


# In[17]:


xlist = tsne_meta.iloc[:,1].tolist()
ylist = tsne_meta.iloc[:,2].tolist()

xticks = [float(min(xlist))-1, float(max(xlist))+1]
yticks = [float(min(ylist))-1, float(max(ylist))+1]

fig, ax = plt.subplots()
plt.scatter(tsne_meta.iloc[:,1],tsne_meta.iloc[:,2],c = cluster_10_ids.iloc[:,1])
plt.title("T-SNE Representation of ~10K Human B-cells")

ylims = ax.get_ylim()
# ylims = tuple([( float(format(j, '.2f'))) for j in ylims])
# print(ylims)

ax.set_yticks(ylims)

xlims = ax.get_xlim()
ax.set_xticks(xlims)

plt.show()


# # Isolating co-expression domains
# 
# Rather than solving a "clustering" problem, we would like to identify contiguous blocks along the diagonal where "squares" exist. This is essentially a Hi-C problem (to call TADs), and we can involve "cooltools" to do so.
# 
# We could think of the problem "algorithmically" as such: Given n points and n^2 relationships to each other, how can we identify m groupings such that the total signal within the...
# 
# 
# The counts are rather sparse for a given set of cells. Most of the genes have 0 counts across all of the cells. See below

# In[18]:


plt.hist(np.log10(np.sum(mat,axis=1)+1),bins=50,density=True)
plt.xlabel(s="log10(total counts across cells + 1)")
plt.ylabel(s="density")
plt.title("Distribution of reads mapping to %d genes\nacross all %d human B-cells"%(mat.shape[0],mat.shape[1]))
plt.show()


# ## Remove zero-sum genes

# First, we will remove genes from the matrix that have zero-sum across all of the cells (since such genes have an "indeterminate" status across the cells

# ## Functions for computing pairwise correlations

# Now, we compute the pairwise correlation of genes

# In[19]:


def filter_genes(inmat,ingenes):
    #first, remove genes that have zero sums all throughout the cells
    mat1_sumxcells = np.sum(inmat,axis=1)
    zerogenes_ind = np.where(mat1_sumxcells == mat1_sumxcells.min())[0].tolist()
    nonzerogenes_ind = np.where(mat1_sumxcells != mat1_sumxcells.min())[0].tolist()

#     mask = ingenes.index.isin(zerogenes_ind)
#     nonzero_genes = genes_df[~mask]

    mask = ingenes.index.isin(nonzerogenes_ind)
    nonzero_genes = genes_df[mask]
    
    return(nonzero_genes)

# print(genes_df.shape)
# print(mat.shape)

def calculate_correlation_matrices(inmat,ingenes):
    #first, remove genes that have zero sums all throughout the cells
    mat1_sumxcells = np.sum(inmat,axis=1)
    nonzerogenes_ind = np.where(mat1_sumxcells != mat1_sumxcells.min())[0].tolist()

    mask = ingenes.index.isin(nonzerogenes_ind)
#     nonzero_genes = genes_df[~mask]
    nonzero_genes = genes_df[mask]
    
    mat1_nonzerogenes = inmat.todense()[nonzerogenes_ind,:]
    mat1_nonzerogenes_sparse = sps.coo_matrix(np.corrcoef(mat1_nonzerogenes))
    
    return(mat1_nonzerogenes_sparse)


# ## Correlation for all cells

# In[21]:


# inputmats = [(mat,genes_df)] # for i in cellgroups_keep]


# mp_pool = mp.Pool(4)
# filtgenes = mp_pool.starmap(filter_genes,inputmats)
# cormats = mp_pool.starmap(calculate_correlation_matrices,inputmats)
# mp_pool.close()

filtgenes = filter_genes(inmat=mat,ingenes=genes_df)
# cormat = calculate_correlation_matrices(ingenes=genes_df,inmat=mat)
print(filtgenes.shape)
# mat


# In[ ]:


mp_pool = mp.Pool(5)
cormat = mp_pool.starmap(calculate_correlation_matrices,[genes_df,mat])
mp_pool.close()


# In[ ]:


cormat.shape


# In[ ]:


cormat.todense()


# We'll use this matrix

# In[ ]:


seaborn.heatmap(cormat.todense()[0:2000,0:2000])
plt.xticks([])
plt.yticks([])
plt.show()


# ## Map entries onto genomic intervals

# Now, we want to convert the correlation matrix to a "list of contacts" so that we can work with it as a cooler file. 

# ### Isolate the GTF genes that are used

# In[ ]:


#filtgenes.geneid


# In[ ]:


hg19gtf_genes_filt = hg19gtf_genes[hg19gtf_genes.gene_id.isin(filtgenes.geneid) & hg19gtf_genes.type.isin(["gene"])]
#filtgenes


# In[ ]:


#re-index the filtgenes to match the indexing of the cormat
filtgenes.index = range(cormat.shape[0])

#re-inded the hg19gtf as well
hg19gtf_genes_filt.index =  range(cormat.shape[0])

#print(filtgenes.shape)
#print(cormat.shape)
#hg19gtf_genes_filt


# In[ ]:


hg19gtf_genes_filt.shape


# ### Get the pair-wise listing of intervals and the correlation

# In[ ]:


#test

for i in hg19gtf_genes_filt.index[0:5]:
    for j in hg19gtf_genes_filt.index[0:5]:
        gene1 = hg19gtf_genes_filt.iloc[i,[0,3,4,6]].tolist()
        gene2 = hg19gtf_genes_filt.iloc[j,[0,3,4,6]].tolist()
        gene1[0] = "chr"+str(gene1[0])
        gene2[0] = "chr"+str(gene2[0])
        corscore = cormat.todense()[i,j]
        outline = [str(n) for n in gene1+gene2+[corscore]]
        print("\t".join(outline))

# for i in range(len(hg19gtf_genes_filt.index[0:5])):
#     for j in range(len(hg19gtf_genes_filt.index[0:5])):
#         gene1 = hg19gtf_genes_filt.iloc[i,[0,3,4,6]].tolist()
#         gene2 = hg19gtf_genes_filt.iloc[j,[0,3,4,6]].tolist()
#         gene1[0] = "chr"+str(gene1[0])
#         gene2[0] = "chr"+str(gene2[0])
#         corscore = cormat.todense()[i,j]
#         outline = [str(n) for n in gene1+gene2+[corscore]]
#         print("\t".join(outline))


# May run faster using https://pandas.pydata.org/pandas-docs/version/0.21/generated/pandas.DataFrame.stack.html. If we want to print out 12,500 x 12,500 = ~ 150 million entries rapidly, we will need to either run this on a computing cluster (Orchestra) or utilize some clever strategy

# In[ ]:


print(cormat.shape)
#we 


# Whatever output we get, we will have to use pairix or tabix to index the coordinates before making a cooler file

# In[ ]:


import time


# In[ ]:


n = 100
gene_chunks = []
for i in range(0,len(hg19gtf_genes_filt.index),n):
    gene_chunks.append(hg19gtf_genes_filt[i:i+n])


# In[ ]:


def create_ouput_line(mat,genesgtf,ind1,ind2):
    gene1 = genesgtf.iloc[ind1,[0,3,6]].tolist()
    gene2 = genesgtf.iloc[ind2,[0,3,6]].tolist()
    gene1[0] = "chr"+str(gene1[0])
    gene2[0] = "chr"+str(gene2[0])
    corscore = mat.todense()[ind1,ind2]
    outline = [str(n) for n in gene1+gene2+[corscore]]
    return("\t".join(outline)+"\n")


# In[ ]:


def save_file(intvls,mat=cormat):
    #since this is a symmetric matrix, we will output only a half-triangle
    output = ""
    for i in intvls.index:
        for j in range(i,intvls.index[-1]):
            output += create_ouput_line(mat=mat,genesgtf=intvls,ind1=i,ind2=j)
    return(output)


# In[ ]:


# mp_pool = mp.pool.ThreadPool(4)

# a = mp_pool.map(save_file,[gene_chunks[0]])
#a = save_file(intvls=gene_chunks[0])


# In[ ]:


mp_pool.close()


# In[ ]:


a.get()


# In[ ]:


from multiprocessing.pool import ThreadPool

def foo(word, number):
    print (word * number)
    r[(word,number)] = number
    return(number)

words = ['hello', 'world', 'test', 'word', 'another test']
numbers = [1,2,3,4,5]
pool = ThreadPool(5)
results = []
for i in range(0, len(words)):
    results.append(pool.apply_async(foo, args=(words[i], numbers[i])))

pool.close()
pool.join()
results = [r.get() for r in results]
print(results)


# In[ ]:


a.get()


# In[ ]:


mp_pool.close()


# In[ ]:


for i in gene_chunks:
    print(i)


# In[ ]:


counter = 1

# for faster writing, https://stackoverflow.com/questions/27384093/fastest-way-to-write-huge-data-in-file

genepairs = []

t0 = time.time()
for i in hg19gtf_genes_filt.index[0:1000]:
    for j in hg19gtf_genes_filt.index[0:1000]:
        gene1 = hg19gtf_genes_filt.iloc[i,[0,3,6]].tolist()
        gene2 = hg19gtf_genes_filt.iloc[j,[0,3,6]].tolist()
        gene1[0] = "chr"+str(gene1[0])
        gene2[0] = "chr"+str(gene2[0])
        corscore = cormat.todense()[i,j]
        outline = [str(n) for n in gene1+gene2+[corscore]]
        genepairs.append(outline)
d = time.time() - t0
print("Took %.2f seconds" %d)


# Perhaps we could split the list into blocks and loop through them in parallel, joining the outputs at the end?

# In[ ]:


with open("./scrnaseq_10x_bcell_correlation_contacts_all.txt",'w') as outfile:
#     for i in hg19gtf_genes_filt.index[0:1000]:
#         for j in hg19gtf_genes_filt.index[0:1000]:
#     for i in range(len(hg19gtf_genes_filt.index[0:1000])):
#         for j in range(len(hg19gtf_genes_filt.index[0:1000])):
    for i in hg19gtf_genes_filt.index[0:5]:
        for j in hg19gtf_genes_filt.index[0:5]:
#             if i != j:
            gene1 = hg19gtf_genes_filt.iloc[i,[0,3,6]].tolist()
            gene2 = hg19gtf_genes_filt.iloc[j,[0,3,6]].tolist()
            gene1[0] = "chr"+str(gene1[0])
            gene2[0] = "chr"+str(gene2[0])
                #print(i,j)
            corscore = cormat.todense()[i,j]
            outline = [str(n) for n in gene1+gene2+[corscore]]
                #print(i,j,outline)
            outfile.write("\t".join(outline)+"\n")
            if counter % 10 == 0 or counter == 1:
                print('%d lines complete'%counter)
            counter += 1


# In[ ]:


#hg19gtf_genes_filt


# In[ ]:


#12700
# hg19gtf_genes_filt.iloc[12700,[0,3,6]].tolist()
# cormat.todense()[0,12700]
cormat.shape


# In[ ]:


# def writeout(filename):
#     with open(filename,'w') as outfile:
#         for i in hg19gtf_genes_filt.index[0:1000]:
#             for j in hg19gtf_genes_filt.index[0:1000]:
#                 gene1 = hg19gtf_genes_filt.iloc[i,[0,3,6]].tolist()
#                 gene2 = hg19gtf_genes_filt.iloc[j,[0,3,6]].tolist()
#                 gene1[0] = "chr"+str(gene1[0])
#                 gene2[0] = "chr"+str(gene2[0])
#                 corscore = cormat.todense()[i,j]
#                 outline = [str(n) for n in gene1+gene2+[corscore]]
#                 outfile.write("\t".join(outline)+"\n")
#     return(filename)


# In[ ]:


# mp_pool = mp.Pool(4)
# filenames = mp_pool.map(writeout,["./scrnaseq_10x_bcell_correlation_contacts_1000.txt"])
mp_pool.close()


# In[ ]:


# for i in range(len(filtgenes.index)):
#     for 


# In[ ]:


#print([i for i in filtgenes.index])


# ## Export matrix

# We want to export the correlation matrix so that we can convert it into a cooler file (and view the results in HiGlass)

# ## Other approaches

# ### Clustering

# The clustering analysis turned out differentially-expressed genes for each cluster, so we could calculate a co-expression matrix using just those genes (rather than those with no counts recorded). Nevertheless, now we can perform matrix operations on the coo_matrix to identify the "blocks" of co-expressed genes.
# 
# This could be a job for "block modelling." https://stats.stackexchange.com/questions/138325/clustering-a-correlation-matrix. However, we want to retain the positional ordering of genes in the matrix

# ### Binning genes

# What we will do is import an ENSEMBL GTF of genes for hg19 and create "blocks" from which we can coarse-grain the gene co-expression matrix (and recalculate correlation between all genes being expressed by getting correlation of reads for each block). We can then perform "cooler" operations on the matrix to isolate co-expression blocks

# # Compare co-expression domains and contact scores from HiC data
# 
# We will perform the following procedure
# 1. For the lower-left diagonal of the matrix, take each diagonal slice parallel to the hypotenuse
#     * Each slice represents some scale of genomic distance between regions. We would like to compute the whether particular contact-coexpression trends stand out over long genomic distances
# 2. Obtain both the contact scores and the expression correlation of all genes in the genomic bins throughout the slice. 
# 3. We ultimately want to calculate some "threshold" at which genomic proximity and gene co-expression are potentially significant. 
# * Weight the signal used to compute correlation by the number of genes present (assign an "uncertainty" to the reads mapping to the genes in the bin)

# In[ ]:


#test_cooler
#ctest_coexp = cooler.Cooler("./cooler_files/test_out.cool")
ctest_coexp = cooler.Cooler("./cooler_files/scrnaseq_10x_bcell_correlation_contacts_1000_5kb_bin.cool")

resolution = ctest_coexp.info['bin-size']
balance_weights = cooler.ice
mat2 = ctest_coexp.matrix(balance=False).fetch('chr1:1-1,000,000')
print(mat2)
#probably need a finer resolution
#mat2 = ctest_coexp.matrix().fetch()
#mat2 = ctest_coexp.matrix()


# In[ ]:


ctest_coexp.info


# In[ ]:


#mat2 = ctest_coexp.matrix(as_pixels=False,balance=False).fetch('chr1:1-15,000,000')


# In[ ]:


plt.matshow(mat2, cmap='YlOrRd')
plt.show()


# In[ ]:


# ctest_hic = cooler.Cooler("../data/coolfiles/GSM2644945_Untreated-R1.20000.cool")
ctest_hic = cooler.Cooler("./hic_data/Rao2014-GM12878-MboI-allreps-filtered.5kb.cool")
print(ctest_hic.info['bin-size'])
#balance_weights = cooler.ice
print(ctest_hic.info)
mat_hic_untreated = ctest_hic.matrix(balance=True).fetch('chr1:1-1,000,000')
plt.matshow(mat_hic_untreated, cmap='YlOrRd')
plt.show()


# Clearly, we cannot compare the two matrices unless we consider regions mutual to both the Hi-C file and the coexpression matrix. But, this would make for an good introductory figure. 

# # OUTDATED Components
# At some point, we will delete the below sections

# # Isolate groups of cells
# We can identify putative groups of cells (from the 10-cluster scheme) and calculate co-expression matrices fo the genes using measurements from each cluster of cells (we will have as any matrices as we have cells). We will calculate these co-expression matrices for the whole set of ~30,000 genes and represent it as a plot. 

# In[ ]:


mat1_sumxcells = np.sum(mat1,axis=1)
zerogenes_ind = np.where(mat1_sumxcells != mat1_sumxcells.min())[0].tolist()

mask = genes_df.index.isin(zerogenes_ind)
nonzero_genes = genes_df[~mask]


# In[ ]:


nonzero_genes.shape


# In[ ]:


def get_cluster(clusterid,celldf,cluster_ids):
    cells = cluster_ids[cluster_ids.Cluster == clusterid]
    cells_ind = cells.index
    return(celldf.tocsr()[:,cells_ind].tocoo())


# In[ ]:


mp_pool = mp.Pool(4)
clusters = [(i,mat,cluster_10_ids) for i in set(cluster_10_ids.iloc[:,1])]
cellgroups = mp_pool.starmap(get_cluster, clusters)
mp_pool.close()


# In[ ]:


mat1 = cellgroups[0]
print(type(mat1))


# ## Calcuate co-ex

# In[ ]:


def calculate_sparsemat_correlation(sparsemat):
    #A = sparsemat + 1e-120 #add an infinitesmal pseudocount #I should probably add random "noise" to preserve the correlations
    A = sparsemat # + We want to input some random noise so that the preponderance of zeroes does not turn up nan-level correlations 
    n = A.shape[1] #gets the number of cells; shape[0[] would give the number of variables
    rowsum = A.sum(1)
    
    #compute covariance matrix
    centering = rowsum.dot(rowsum.T.conjugate()) / n
    C = (A.dot(A.T.conjugate()) - centering) / (n - 1)

    # The correlation coefficients are given by
    # C_{i,j} / sqrt(C_{i} * C_{j})
    d = np.diag(C)
    coeffs = C / np.sqrt(np.outer(d, d))

    return C#coeffs


# In[ ]:


# from scipy.sparse import random
# from scipy import stats
# class CustomRandomState(object):
#     def randint(self, k):
#         i = np.random.randint(k)
#         return i - i % 2
# rs = CustomRandomState()
# rvs = stats.poisson(2, loc=1).rvs
# S = random(1000, 10000, density=0.25, random_state=rs, data_rvs=rvs)
# print(S.A)
# print(np.corrcoef(S.A))


# In[ ]:


# def calculate_corr(pd_df):
#     return pd_df.T.corr(method="pearson") #pearson is default


# In[ ]:


for i in cellgroups:
    print(i.shape)


# In[ ]:


#keey the clusters that have more than 'ncell' cells
ncell = 1
inds = [None for i in range(len(cellgroups))]
for i in range(len(cellgroups)):
    if cellgroups[i].shape[1] > ncell:
        inds[i] = i
inds = [i for i in inds if i is not None]
print(inds)
#corrdfs_actual = [corrdfs[i] for i in inds]
cellgroups_keep = [cellgroups[i] for i in inds]
for i in [j for j in cellgroups_keep]:
    print(i.shape)


# In[ ]:


# https://stackoverflow.com/questions/44553858/compute-a-pairwise-distance-matrix-is-a-scalable-big-data-ready-approach-avail
# https://stackoverflow.com/questions/19231268/correlation-coefficients-for-sparse-matrix-in-python
# https://stackoverflow.com/questions/34872854/how-to-ignore-zeros-when-calculating-correlations-between-columns-for-sparse-mat


# In[ ]:


mat1 = cellgroups_keep[0]
print(mat1.shape)


# ### Outdated: correlations with the full set (TOO LONG!)

# In[ ]:


# convert in a numpy array
m = np.array(mat1.todense())

# create the distance matrix using pdist
m1 = m + 1e-119 #add an infinitesmal to pseudocount 

sparsecorr = calculate_sparsemat_correlation(m1)


# In[ ]:


scipy.sparse.random()


# In[ ]:


d = ds.pdist(m[1:100,1:100], 'correlation')

d= ds.squareform(d)
print(np.amax(np.sum(np.corrcoef(m1[1:100,1:100]),axis=1)))


# In[ ]:


mp_pool = mp.Pool(4)

#cellgroups_keep_T = [i.transpose() for i in cellgroups_keep]
# corrdfs = mp_pool.map(calculate_corr,cellgroups_keep)
corrdfs = mp_pool.map(calculate_sparsemat_covariance,cellgroups_keep)


print(corrdfs[0])
#print(cluster_0)
#print(cluster_0.shape)
#cluster_0.corr()
#cov_cluster_0 = calculate_sparsemat_correlation(cluster_0)
#cluster_0.sum(1).dot(cluster_0.sum(1).T.conjugate()).shape


# In[ ]:


# Now, for associate each position on the correlation matrix with a gene from the 


# In[ ]:


#corrdfs_actual[3]


# In[ ]:


def plot_coo_matrix(m):
#     if not isinstance(m, scipy.sparse.coo_matrix):
#         m = scipy.sparse.coo_matrix(m)
    ax = plt.figure()
#    ax = fig.add_subplot(111, axisbg='black')
#    ax.plot(m.col, m.row, 's', ms=1)
#    ax.set_xlim(0, m.shape[1])
#    ax.set_ylim(0, m.shape[0])
#    ax.set_aspect('equal')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.stet_xticks([])
    ax.set_yticks([])
    return ax


# In[ ]:


a = plot_coo_matrix(corrdfs[0])
a.figure


# In[ ]:


# #from https://stackoverflow.com/questions/19231268/correlation-coefficients-for-sparse-matrix-in-python
# def sparse_corrcoef(A, B=None):

#     if B is not None:
#         A = sparse.vstack((A, B), format='csr')

#     A = A.astype(np.float64)
#     n = A.shape[1]
#     print(n)

#     # Compute the covariance matrix
#     rowsum = A.sum(1)
#     centering = rowsum.dot(rowsum.T.conjugate()) / n
#     C = (A.dot(A.T.conjugate()) - centering) / (n - 1)

#     # The correlation coefficients are given by
#     # C_{i,j} / sqrt(C_{i} * C_{j})
#     d = np.diag(C)
#     coeffs = C / np.sqrt(np.outer(d, d))

#     return coeffs


# In[ ]:


# for i in cellgroups:
#     print(i.shape)
# #cellgroups[4]
# #mat_csr = mat.tocsr()
# #mat_csr[2,cellgroups[0].tolist()]
# #    mat_group = mat.tocsc[0:100] cells_ind


# In[ ]:


# test_corr = sparse_corrcoef(A=cellgroups[0])


# In[ ]:


#import statsmodels


# In[ ]:


test_corr_nona = np.nan_to_num(test_corr) #replace all NA with 0


# In[ ]:


test_corr_sparse = scipy.sparse.coo_matrix(test_corr)
#print(type(test_corr_sparse))

#ax = sparse_corrcoef(test_corr)
#ax.figure


# In[ ]:


print(np.amax(np.sum(a=test_corr_nona,axis=0))) #so this does work
test_corr_nona_sparse = scipy.sparse.coo_matrix(test_corr_nona) #store as a sparce matrix


# In[ ]:


#import cooler
test_corr_nona_sparse.todense()


# ### Computing correlations for clusters of cells -- DEPRECATED!

# In[ ]:


mp_pool = mp.Pool(4)
inputmats = [(i,genes_df) for i in cellgroups_keep]
#print(inputmats)
filtgenes = mp_pool.starmap(filter_genes,inputmats)
#cormats = mp_pool.starmap(calculate_correlation_matrices,inputmats)
mp_pool.close()


# In[ ]:


for i in cellgroups_keep:
    print(i)


# In[ ]:


#mp_pool = mp.Pool(4)
#cormats = mp_pool.starmap_async(calculate_correlation_matrices,inputmats)
mp_pool.close()


# In[ ]:


cormats = []
for i in inputmats:
    cormats.append(calculate_correlation_matrices(i[0],i[1]))


# In[ ]:


for i in cormats:
    print(i.shape)


# In[ ]:


cormats[0].todense()


# In[ ]:


cormats[1].todense()


# We can also perform a distance calculation

# In[ ]:


d = ds.pdist(mat1_nonzerogenes_sparse.todense(), 'correlation')

