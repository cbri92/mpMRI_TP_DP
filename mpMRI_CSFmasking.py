# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 16:16:05 2020

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
from scipy.stats.stats import pearsonr
import dicom2nifti
from multiprocessing.pool import ThreadPool
from functools import partial
from ImageAnalysisFunctions import *

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing individual patients directories' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
n_subj = len(subjs_name) #Total number of subjects


#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
   
    nii_dir = subj_dir +'/nii/' 
   
    #Set paths to subfolders    

    bet_imgs = nii_dir +'BETimages'
    
    #Read images
    
    # T1bet = sitk.ReadImage(bet_imgs +'/T1_bet.nii') #read new T1 bet
    brain_mask = sitk.ReadImage(bet_imgs+'/brain_mask.nii')
    brain_mask = sitk.Cast(brain_mask, sitk.sitkUInt8)
    CSF_3 = sitk.ReadImage(bet_imgs+'/FLAIR_bet_seg_3.nii')
    CSF_0 = sitk.ReadImage(bet_imgs+'/T1CE_bet_seg_0.nii')
    CSF_1 = sitk.ReadImage(bet_imgs+'/T1CE_bet_seg_1.nii')
    CSF_full = CSF_0+CSF_1+CSF_3
    CSF_full = sitk.BinaryThreshold(CSF_full, 1.0, 4.0, 1, 0)
    CSF_full = sitk.Cast(CSF_full, sitk.sitkUInt16)
    Tum = sitk.ReadImage(bet_imgs+'/Tumour_mask.nii')
    GTV = sitk.ReadImage(bet_imgs+'/GTV.nii')
    GTV = sitk.Cast(GTV, sitk.sitkUInt16)
    
    #Generate a binary mask excluding the tumour
    Tum = Tum + GTV
    Tum_excl = sitk.BinaryThreshold(Tum, 1.0, 3.0, 0, 1)
    Tum_excl = sitk.Cast(Tum_excl, sitk.sitkUInt16)
    # sitk.WriteImage(Tum_excl, bet_imgs+'/Tum_excl.nii')
    
    #Generate a binary mask of CSF excluding tumour
    CSF = CSF_full*Tum_excl
    CSF = sitk.BinaryThreshold(CSF, 1.0, 2.0, 0, 1)
    CSF = CSF*brain_mask
    CSF = sitk.BinaryMorphologicalOpening(CSF, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
    sitk.WriteImage(CSF, bet_imgs+'/CSF.nii')
