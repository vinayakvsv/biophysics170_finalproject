
# coding: utf-8

# In[1]:


#import basics
import os
import sys
import re

#import Python bioinformatics packages
import Bio
import matplotlib.pyplot as plt
import numpy as np
import pandas
import h5py

#import Hi-C packages and other lab-made tools
import cooler

#set the data directory
dir_path = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(dir_path,"../data")
print(DATA_DIR)

