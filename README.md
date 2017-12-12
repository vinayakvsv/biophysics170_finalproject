# Biophysics 170 Fall 2017 Final Project
## Introduction
This repository contains files for the Biophysics 170 Fall 2017 final project at Harvard Medical School. 

## Summary of project
Our project aims to answer the following questions
1. What is the extent of overlap between TADs and spatial domains of co-expression
2. How do overlaps of TADs and co-expression modules (CoEx) change with depletion of CTCF and cohesion?
3. How do co-expressed genes with closeby eQTLs ("variable chromatin modules") overlap with TADs?
4. Do Mendelian disease-associated variants overlap with TAD boundaries or CTCF occupancy regions? Do disease-associated eQTL's change the structure of TAD-CoEx or TAD-VCM overlaps?

### Project descriptions from biophysics 170 website
1. Co-expressed genes and eQTLs tend to be close to each other along a human chromosome, forming VCMs (variable chromatin modules see paper). The project would be to test whether VCMs and/or co-expressed genes are located within TADs (compartments and TADs).  We can get you VCM data or single-cell expression data, e.g. from here https://www.10xgenomics.com/single-cell/ [Vinay]
2. TAD boundaries have been shown to play important functional roles in controlling ranges of enhancers activity, and were found to be associated with heritable diseases (PMID 25959774). The project would be to find GWAS peaks or other clinical variants or Mendelian disease-associated variants that can overlay TAD boundaries and/or CTCF occupancy regions. [Noah]
3. In a recent papers we examined the effect of cohesin loss and CTCF loss on chromosome organization. In both cases we observed mild effects on transcription (~1000 genes changing their expression). That was surprising because CTCF and cohesin a believed to be central for enhancer-promoter interactions. It possible that housekeeping genes are less dependent on such interactions, while developmental genes are more dependent. The project is to examine changes in expression upon cohesin loss (and CTCF loss, if time allows) and to test whether housekeeping are less affected than developmental genes by cohesin loss. This project can in principle can be done without coding, but browsing Hi-C and transcription data in higalss.io browser. [Vinay--done as a complement to Topic 2]


## Repository Structure
*To be added*

# Setup
1. We have a number of python packages to install
2. For visualization, we will be running HiGlass from a docker container. See https://github.com/hms-dbmi/higlass/wiki for instructions
3. We will be using the Dockerfile image from the 4DN Nucleosome consortium for HiC data. 
4. Jupyter setup: https://stackoverflow.com/questions/35254852/how-to-change-jupyter-start-folder and http://jupyter.readthedocs.io/en/latest/install.html



## Notes
1. Remove all existing docker containers:
```
https://zaiste.net/posts/removing_docker_containers/
```
2. Install Cooler: https://github.com/mirnylab/cooler
    + sudo pip install --ignore-installed six --upgrade cooler (if you have El Capitan or higher for Mac)
    + https://github.com/mirnylab/cooler 

# Analysis Strategy
*to be added

# Important links

HiGlass: https://github.com/hms-dbmi/higlass/wiki
Cooler: https://github.com/mirnylab/cooler
Hi-C visualization page: http://promoter.bx.psu.edu/hi-c/view.php

# Parallelization

We use multiprocessing package. However, some of these jobs will have to be run on Orchestra (for example, generating the correlation map of the B-cell co-expression data)

