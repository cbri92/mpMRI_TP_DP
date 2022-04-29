# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:01:24 2022

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
from radiomics import featureextractor
import six, numpy as np
from statistics import pvariance
from statistics import mean
from ImageStatisticsFunctions import *
from scipy.stats import normaltest
import seaborn as sns
import pingouin as pg

def reassign_voxel_intensity(image, mask, low_thr, high_thr, new_value):
        
    '''This function set the intensity of all the voxels with original intensity within the low_thr and high_thr to the new_value'''
    masked_img = generate_mask(image, mask)
    mask0 =  generate_thresholded_roi(masked_img, low_thr, high_thr)
    masked_img = set_mask_value(masked_img, mask0, new_value)
    masked_img = generate_mask(masked_img, mask)
    return masked_img

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing individual patients directories' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
     
    subj_dir1 = data_supradir+current+'/TP1/'
    subj_dir2 = data_supradir+current+'/TP2/'
    
    
    print('Reading images for '+current) 
        
    for f in os.scandir(subj_dir1):
        if f.path.find('Inv_ADC') != -1:
            DPresc1_path = f.path
        elif f.path.find('GTV') != -1:
            GTV1_path = f.path
        elif f.path.find('CTV') != -1:
            CTV1_path = f.path
        elif f.path.find('PTV') != -1:
            PTV1_path = f.path 
            
    for f in os.scandir(subj_dir2):
        if f.path.find('Inv_ADC') != -1:
            DPresc2_path = f.path

    #%%Read images
    
    DPresc1 = sitk.ReadImage(DPresc1_path)
    DPresc1 = 80-DPresc1
    GTV1 = sitk.ReadImage(GTV1_path)
    CTV1 = sitk.ReadImage(CTV1_path)
    PTV1 = sitk.ReadImage(PTV1_path)
    
    DPresc2 = sitk.ReadImage(DPresc2_path)
    DPresc2 = 80-DPresc2

#%% Generate masked dose volume

    GTV1 = GTV1>1
    CTV1 = CTV1>1
    PTV1 = PTV1>1
    
    #Set same origin for targets and dose grids
    GTV1.SetOrigin(DPresc1.GetOrigin())
    CTV1.SetOrigin(DPresc1.GetOrigin())
    PTV1.SetOrigin(DPresc1.GetOrigin())
            
    print('Generating dose masked to targets for '+current)   
    
    GTV1_dose = reassign_voxel_intensity(DPresc1, GTV1,0,0,80)
    CTV1_dose = reassign_voxel_intensity(DPresc1, CTV1,0,0,80)
    PTV1_dose = reassign_voxel_intensity(DPresc1, PTV1,0,0,80)
    
    GTV2_dose = reassign_voxel_intensity(DPresc2, GTV1,0,0,80)
    CTV2_dose = reassign_voxel_intensity(DPresc2, CTV1,0,0,80)
    PTV2_dose = reassign_voxel_intensity(DPresc2, PTV1,0,0,80)   
    
    GTV_sub = GTV1_dose-GTV2_dose
    CTV_sub = CTV1_dose-CTV2_dose
    PTV_sub = PTV1_dose-PTV2_dose
    
    #%%Save images
    
    sitk.WriteImage(GTV1_dose, subj_dir1+'Masked_DPresc_GTV.nii')
    sitk.WriteImage(CTV1_dose, subj_dir1+'Masked_DPresc_CTV.nii')
    sitk.WriteImage(PTV1_dose, subj_dir1+'Masked_DPresc_PTV.nii')
    sitk.WriteImage(GTV2_dose, subj_dir2+'Masked_DPresc_GTV.nii')
    sitk.WriteImage(CTV2_dose, subj_dir2+'Masked_DPresc_CTV.nii')
    sitk.WriteImage(PTV2_dose, subj_dir2+'Masked_DPresc_PTV.nii')
    
    sitk.WriteImage(GTV_sub, data_supradir+current+'/DPresc_Diff_GTV.nii')
    sitk.WriteImage(CTV_sub, data_supradir+current+'/DPresc_Diff_CTV.nii')
    sitk.WriteImage(PTV_sub, data_supradir+current+'/DPresc_Diff_PTV.nii')