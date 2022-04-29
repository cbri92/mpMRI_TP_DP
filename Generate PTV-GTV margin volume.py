# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 15:28:42 2021

@author: cbri3325
"""

#%% Import functions 

import matplotlib.pyplot as plt
import matplotlib as mpl
import SimpleITK as sitk
import numpy as np
import pandas as pd
import datetime
import os
import glob
import gzip
import shutil
import xlsxwriter
import time
from scipy.stats.stats import pearsonr
import dicom2nifti
from multiprocessing.pool import ThreadPool
from functools import partial
from ImageAnalysisFunctions import *

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing individual patients directories' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
    
    rel_imgs = subj_dir+'/nii/Relative images/'
    BET_imgs = subj_dir+'/nii/BET images/'
    PTV = sitk.ReadImage(rel_imgs+'PTV.nii')
    GTV = sitk.ReadImage(rel_imgs+'GTV.nii')
    brain_mask = sitk.ReadImage(BET_imgs+'brain_mask.nii')
    
    Margin_volume=PTV-GTV
    Margin_volume=generate_mask(Margin_volume, brain_mask)
    sitk.WriteImage(Margin_volume, rel_imgs+'Margin_volume.nii')
