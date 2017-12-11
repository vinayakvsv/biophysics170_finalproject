
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


# In[4]:


#set import functions

def import_matrix_mtx_file(filename,aspd=False):
    inmat = scipy.io.mmread(source=filename)
    if (aspd==True):
        return(pd.SparseDataFrame(inmat))
    else:
        return(inmat)

#how can we generate a SparseDataFrame from a Sparse Matrix without creating a dense matrix in memory?

# def import_csv(filename,delimiter="\t"):
#     ids = [row[0] for row in csv.reader(open(filename), delimiter=delimiter)]
#     names = [row[1] for row in csv.reader(open(filename), delimiter=delimiter)]
#     return(ids,names)


# In[5]:


#import the file. Use multiprocessing to make this run smoother

inmatrix = os.path.join(bcell_matrix_dir, "matrix.mtx")
genes_path = os.path.join(bcell_matrix_dir, "genes.tsv")
barcodes_path = os.path.join(bcell_matrix_dir, "barcodes.tsv")


# In[6]:


#impor the matrix file
mp_pool = mp.Pool(5)
infile = mp_pool.map(import_matrix_mtx_file, [inmatrix])
mp_pool.close()


# In[7]:


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

# In[8]:


# load the gene IDs and names
gene_ids = [row[0] for row in csv.reader(open(genes_path), delimiter="\t")]
gene_names = [row[1] for row in csv.reader(open(genes_path), delimiter="\t")]


# In[9]:


# load the barcode IDs
barcodes = [row[0] for row in csv.reader(open(barcodes_path), delimiter="\t")]


# In[10]:


genes_df = pd.DataFrame()
genes_df = genes_df.assign(genename=pd.Series(gene_names),geneid=pd.Series(gene_ids))


# for multi-core jobs
# https://stackoverflow.com/questions/33480675/multiprocessing-on-a-set-number-of-cores

# ## Import ENSEMBL GTF

# In[11]:


#set ensembl file
#scrnaseq_dir = "../scrnaseq_10x/"
hg19gtf = "../hg19ref/refdata-cellranger-hg19-1.2.0/genes/genes.gtf"
hg19gtf_dir = os.path.join(scrnaseq_dir,hg19gtf)

hg19gtf_genes = pd.read_csv(hg19gtf,delimiter="\t",comment="#",header=None)


# In[12]:


hg19gtf_genes.columns = ["chr","database","type","start","end","other","strand","other2","metadata"]


# In[13]:


metadata = hg19gtf_genes.metadata.apply(lambda x: re.split(r"[a-zA-Z]+_[a-zA-Z]+ \"|\";",x)[1])


# In[14]:


hg19gtf_genes['gene_id'] = pd.Series(metadata)
#hg19gtf_genes
hg19gtf_genes_expression = hg19gtf_genes[hg19gtf_genes.gene_id.isin(gene_ids) & hg19gtf_genes.type.isin(["gene"])]
print(hg19gtf_genes_expression.shape)


# In[ ]:


#how do we efficiently identify blocks? Bin the genome into 100 kb blocks and calculate correlations between genes within them?

#first, let's get a distribution of correlation values:
#n, bins, patches = plt.hist(test_corr_nona[0:100][0:100])
#plt.show()


# # Visualization
# Let's view the cell clusters

# We now load T-SNE metadata for this population. We will show what the populations look like and then isolate the cells responsible to generate co-expression blocks

# In[15]:


#tsne_data = scipy.io.mmread(os.path.join(tsne_dir, "projection.csv")) 
tsne_data = os.path.join(tsne_dir, "projection.csv")
tsne_points = [row for row in csv.reader(open(tsne_data), delimiter=",")]
tsne_meta = pd.DataFrame(tsne_points[1:])
tsne_meta.columns = tsne_points[0]


# In[16]:


#import clustering
cluster_10 = os.path.join(kmeans_dir,"10_clusters")
print(cluster_10)
#get the cluster designations
cluster_10_ids = pd.read_csv(os.path.join(cluster_10,"clusters.csv"))
#cluster_10_ids
#kmeans_10 = os.path.join(pca_dir,"10_clusters") 


# In[17]:


plt.figure(figsize=(10, 10))
plt.scatter(tsne_meta.iloc[:,1],tsne_meta.iloc[:,2],c = cluster_10_ids.iloc[:,1])
plt.xticks([])
plt.yticks([])
plt.show()


# # Isolate groups of cells
# Now, we want to identify putative groups of cells (from the 10-cluster scheme) and calculate co-expression matrices fo the genes using measurements from each cluster of cells (we will have as any matrices as we have cells). We will calculate these co-expression matrices for the whole set of ~30,000 genes and represent it as a plot. 

# In[ ]:


# cellgroups = [None for i in clusters]
# for i in clusters:
#     cells = cluster_10_ids[cluster_10_ids.Cluster == i] #get the cells in the cluster "i"
#     cells_ind = cells.index
#     mat_csr = mat.tocsr()
#     cellgroups[i-1] = scipy.sparse.coo_matrix(mat_csr[:,cells_ind.tolist()]) #cells_ind  


# In[19]:


# def get_cluster(clusterid,cellmat,cluster_ids):
#     cells = cluster_ids[cluster_ids.Cluster == clusterid]
#     cells_ind = cells.index
#     mat_csr = mat.tocsr()
#     return(scipy.sparse.coo_matrix(mat_csr[:,cells_ind.tolist()]))

def get_cluster(clusterid,celldf,cluster_ids):
    cells = cluster_ids[cluster_ids.Cluster == clusterid]
    cells_ind = cells.index
#     mat_cells_ind = celldf.iloc[:,cells_ind.tolist()]  
    return(celldf.tocsr()[:,cells_ind].tocoo())
#     return(mat_cells_ind)
#    mat_csr = mat.tocsr()
#    return(scipy.sparse.coo_matrix(mat_csr[:,cells_ind.tolist()]))


# In[20]:


mp_pool = mp.Pool(4)
clusters = [(i,mat,cluster_10_ids) for i in set(cluster_10_ids.iloc[:,1])]
cellgroups = mp_pool.starmap(get_cluster, clusters)
mp_pool.close()


# In[21]:


mat1 = cellgroups[0]
print(type(mat1))


# ## Calcuate co-ex

# In[22]:


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


# In[23]:


for i in cellgroups:
    print(i.shape)


# In[24]:


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


# In[25]:


mat1 = cellgroups_keep[0]
print(mat1.shape)


# # Isolating co-expression domains
# 
# Rather than solving a "clustering" problem, we would like to identify contiguous blocks along the diagonal where "squares" exist. This is essentially a Hi-C problem (to call TADs), and we can involve "cooltools" to do so.
# 
# We could think of the problem "algorithmically" as such: Given n points and n^2 relationships to each other, how can we identify m groupings such that the total signal within the...
# 
# 
# The counts are rather sparse for a given set of cells. Most of the genes have 0 counts across all of the cells. See below

# In[26]:


plt.hist(np.log10(np.sum(mat1,axis=1)+1),bins=1000,density=True)
plt.xlabel(s="log10(total counts across cells + 1)")
plt.ylabel(s="density")
plt.show()


# ## Remove zero-sum genes

# First, we will remove genes from the matrix that have zero-sum across all of the cells (since such genes have an "indeterminate" status across the cells

# In[27]:


mat1_sumxcells = np.sum(mat1,axis=1)
zerogenes_ind = np.where(mat1_sumxcells != mat1_sumxcells.min())[0].tolist()

mask = genes_df.index.isin(zerogenes_ind)
nonzero_genes = genes_df[~mask]


# In[29]:


nonzero_genes.shape


# Now, we compute the pairwise correlation of genes

# In[36]:


def filter_genes(imat,ingenes):
    #first, remove genes that have zero sums all throughout the cells
    mat1_sumxcells = np.sum(mat1,axis=1)
    zerogenes_ind = np.where(mat1_sumxcells != mat1_sumxcells.min())[0].tolist()

    mask = ingenes.index.isin(zerogenes_ind)
    nonzero_genes = genes_df[~mask]

    return(nonzero_genes)

def calculate_correlation_matrices(inmat,ingenes):
    #first, remove genes that have zero sums all throughout the cells
    mat1_sumxcells = np.sum(mat1,axis=1)
    zerogenes_ind = np.where(mat1_sumxcells != mat1_sumxcells.min())[0].tolist()

    mask = ingenes.index.isin(zerogenes_ind)
    nonzero_genes = genes_df[~mask]
    
    mat1_nonzerogenes = mat1.todense()[zerogenes_ind,:]
    mat1_nonzerogenes_sparse = sps.coo_matrix(np.corrcoef(mat1_nonzerogenes))
    
    return(mat1_nonzerogenes_sparse)


# In[ ]:


# #print(type(mat1_sumxcells))
# mat1_nonzerogenes = mat1.todense()[zerogenes_ind,:]
# mat1_nonzerogenes_sparse = sps.coo_matrix(np.corrcoef(mat1_nonzerogenes))


# In[57]:


mp_pool.close()


# In[43]:


mp_pool = mp.Pool(4)
inputmats = [(i,genes_df) for i in cellgroups_keep]
#print(inputmats)
filtgenes = mp_pool.starmap(filter_genes,inputmats)
#cormats = mp_pool.starmap(calculate_correlation_matrices,inputmats)
mp_pool.close()


# In[52]:


for i in cellgroups_keep:
    print(i)


# In[45]:


#mp_pool = mp.Pool(4)
#cormats = mp_pool.starmap_async(calculate_correlation_matrices,inputmats)
mp_pool.close()


# In[47]:


cormats = []
for i in inputmats:
    cormats.append(calculate_correlation_matrices(i[0],i[1]))


# In[53]:


for i in cormats:
    print(i.shape)


# In[51]:


cormats[0].todense()


# In[50]:


cormats[1].todense()


# We can also perform a distance calculation

# In[ ]:


d = ds.pdist(mat1_nonzerogenes_sparse.todense(), 'correlation')


# In[ ]:


print(d)


# In[ ]:


print(mat1_nonzerogenes_sparse.todense())


# In[54]:


# for i in cormats:
# #     seaborn.heatmap(mat1_nonzerogenes_sparse.todense()[0:1000,0:1000])
#     seaborn.heatmap(i.todense())
#     plt.xticks([])
#     plt.yticks([])
#     plt.show()
# #plt.plot(mat1_nonzerogenes_sparse.todense()[0:100][0:100])
# #plt.show()


# ## Correlation for all cells

# In[59]:


# inputmats = [(mat,genes_df)] # for i in cellgroups_keep]


# mp_pool = mp.Pool(4)
# filtgenes = mp_pool.starmap(filter_genes,inputmats)
# cormats = mp_pool.starmap(calculate_correlation_matrices,inputmats)
# mp_pool.close()

filtgenes = filter_genes(imat=mat,ingenes=genes_df)
cormat = calculate_correlation_matrices(ingenes=genes_df,inmat=mat)

# mat


# In[61]:


cormat.shape


# In[75]:


cormat.todense()


# We'll use this matrix

# In[79]:


seaborn.heatmap(cormat.todense()[0:2000,0:2000])
plt.xticks([])
plt.yticks([])
plt.show()


# ## Map entries onto genomic intervals

# Now, we want to convert the correlation matrix to a "list of contacts" so that we can work with it as a cooler file. 

# ### Isolate the GTF genes that are used

# In[88]:


#filtgenes.geneid


# In[89]:


hg19gtf_genes_filt = hg19gtf_genes[hg19gtf_genes.gene_id.isin(filtgenes.geneid) & hg19gtf_genes.type.isin(["gene"])]
#filtgenes


# In[124]:


filtgenes


# In[125]:


hg19gtf_genes_filt.shape


# ### Get the pair-wise listing of intervals and the correlation

# In[126]:


#test
for i in range(len(hg19gtf_genes_filt.index[0:5])):
    for j in range(len(hg19gtf_genes_filt.index[0:5])):
        gene1 = hg19gtf_genes_filt.iloc[i,[0,3,4,6]].tolist()
        gene2 = hg19gtf_genes_filt.iloc[j,[0,3,4,6]].tolist()
        gene1[0] = "chr"+str(gene1[0])
        gene2[0] = "chr"+str(gene2[0])
        corscore = cormat.todense()[i,j]
        outline = [str(n) for n in gene1+gene2+[corscore]]
        print("\t".join(outline))


# May run faster using https://pandas.pydata.org/pandas-docs/version/0.21/generated/pandas.DataFrame.stack.html

# Whatever output we get, we will have to use pairix or tabix to index the coordinates before making a cooler file

# In[130]:


with open("./scrnaseq_10x_bcell_correlation_contacts_100.txt",'w') as outfile:
#     for i in hg19gtf_genes_filt.index[0:1000]:
#         for j in hg19gtf_genes_filt.index[0:1000]:
    for i in range(len(hg19gtf_genes_filt.index[0:10])):
        for j in range(len(hg19gtf_genes_filt.index[0:10])):
            if i != j:
                gene1 = hg19gtf_genes_filt.iloc[i,[0,3,6]].tolist()
                gene2 = hg19gtf_genes_filt.iloc[j,[0,3,6]].tolist()
                gene1[0] = "chr"+str(gene1[0])
                gene2[0] = "chr"+str(gene2[0])
                #print(i,j)
                corscore = cormat.todense()[i,j]
                outline = [str(n) for n in gene1+gene2+[corscore]]
                #print(i,j,outline)
                outfile.write("\t".join(outline)+"\n")


# In[121]:


#hg19gtf_genes_filt


# In[117]:


#12700
# hg19gtf_genes_filt.iloc[12700,[0,3,6]].tolist()
# cormat.todense()[0,12700]
cormat.shape


# In[107]:


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


# In[109]:


# mp_pool = mp.Pool(4)
# filenames = mp_pool.map(writeout,["./scrnaseq_10x_bcell_correlation_contacts_1000.txt"])
mp_pool.close()


# In[83]:


# for i in range(len(filtgenes.index)):
#     for 


# In[84]:


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


# # Compare co-expression domains and contact scores from HiC data
# 
# We will perform the following procedure
# 1. For the lower-left diagonal of the matrix, take each diagonal slice parallel to the hypotenuse
#     * Each slice represents some scale of genomic distance between regions. We would like to compute the whether particular contact-coexpression trends stand out over long genomic distances
# 2. Obtain both the contact scores and the expression correlation of all genes in the genomic bins throughout the slice. 
# 3. We ultimately want to calculate some "threshold" at which genomic proximity and gene co-expression are potentially significant. 
# * Weight the signal used to compute correlation by the number of genes present (assign an "uncertainty" to the reads mapping to the genes in the bin)

# In[244]:


#test_cooler
ctest_coexp = cooler.Cooler("./test_out.cool")
ctest_coexp = cooler.Cooler("./scrnaseq_10x_bcell_correlation_contacts_1000_nonintvl.100kb_bin.cool")

resolution = ctest_coexp.info['bin-size']
balance_weights = cooler.ice
mat2 = ctest_coexp.matrix(balance=False).fetch('chr1:1-1000000')
print(mat2)
#probably need a finer resolution
#mat2 = ctest_coexp.matrix().fetch()
#mat2 = ctest_coexp.matrix()


# In[170]:


#mat2 = ctest_coexp.matrix(as_pixels=False,balance=False).fetch('chr1:1-15,000,000')


# In[245]:


plt.matshow(mat2, cmap='YlOrRd')
plt.show()

