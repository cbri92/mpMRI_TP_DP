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

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing individual patients directories' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
     
    subj_dir1 = data_supradir+current+'/TP1/'
    subj_dir2 = data_supradir+current+'/TP2/'
    
    print('Gunzip and rename files for '+current)
    
    #Unzip all nii.gz files and delete original gzipped files
    for file in glob.glob(subj_dir1 +'*.gz'):
        gunzip_shutil(file, file[:-3])
        
    for filename in os.listdir(subj_dir1):
        os.chdir(subj_dir1)
        os.rename(filename, filename[7:])
        
    for file in glob.glob(subj_dir2 +'*.gz'):
        gunzip_shutil(file, file[:-3])
        
    for filename in os.listdir(subj_dir2):
        os.chdir(subj_dir2)
        os.rename(filename, filename[7:])
    
    print('Reading images for '+current) 
        
    for f in os.scandir(subj_dir1):
        if f.path.find('DPBN76') != -1:
            Dose1_path = f.path
        elif f.path.find('GTV') != -1:
            GTV1_path = f.path
        elif f.path.find('CTV') != -1:
            CTV1_path = f.path
        elif f.path.find('PTV') != -1:
            PTV1_path = f.path 
            
    for f in os.scandir(subj_dir2):
        if f.path.find('DPBN76') != -1:
            Dose2_path = f.path
        # elif f.path.find('GTV') != -1:
        #     GTV2_path = f.path
        # elif f.path.find('CTV') != -1:
        #     CTV2_path = f.path
        # elif f.path.find('PTV') != -1:
        #     PTV2_path = f.path 
    
    #%%Read images
    
    Dose1 = sitk.ReadImage(Dose1_path)
    GTV1 = sitk.ReadImage(GTV1_path)
    CTV1 = sitk.ReadImage(CTV1_path)
    PTV1 = sitk.ReadImage(PTV1_path)
    
    Dose2 = sitk.ReadImage(Dose2_path)
    # GTV2 = sitk.ReadImage(GTV2_path)
    # CTV2 = sitk.ReadImage(CTV2_path)
    # PTV2 = sitk.ReadImage(PTV2_path)
    
#%% Generate masked dose volume

    GTV1 = GTV1>1
    # sitk.WriteImage(GTV1, subj_dir1+'/GTV_new.nii')
    CTV1 = CTV1>1
    # sitk.WriteImage(CTV1, subj_dir1+'/CTV_new.nii')
    PTV1 = PTV1>1
    # sitk.WriteImage(PTV1, subj_dir1+'/PTV_new.nii')
    
    # GTV2 = GTV2>1
    # # sitk.WriteImage(GTV2, subj_dir2+'/GTV_new.nii')
    # CTV2 = CTV2>1
    # # sitk.WriteImage(CTV2, subj_dir2+'/CTV_new.nii')
    # PTV2 = PTV2>1
    # # sitk.WriteImage(PTV2, subj_dir2+'/PTV_new.nii')
    
    #Set same origin for targets and dose grids
    GTV1.SetOrigin(Dose1.GetOrigin())
    CTV1.SetOrigin(Dose1.GetOrigin())
    PTV1.SetOrigin(Dose1.GetOrigin())
            
    print('Generating dose masked to targets for '+current)   
    
    GTV1_dose = generate_mask(Dose1, GTV1)
    CTV1_dose = generate_mask(Dose1, CTV1)
    PTV1_dose = generate_mask(Dose1, PTV1)
    
    GTV2_dose = generate_mask(Dose2, GTV1)
    CTV2_dose = generate_mask(Dose2, CTV1)
    PTV2_dose = generate_mask(Dose2, PTV1)
    
    GTV_sub = GTV1_dose-GTV2_dose
    CTV_sub = CTV1_dose-CTV2_dose
    PTV_sub = PTV1_dose-PTV2_dose
    
    #%%Save images
    
    sitk.WriteImage(GTV1_dose, subj_dir1+'Masked_Dose_GTV.nii')
    sitk.WriteImage(CTV1_dose, subj_dir1+'Masked_Dose_CTV.nii')
    sitk.WriteImage(PTV1_dose, subj_dir1+'Masked_Dose_PTV.nii')
    sitk.WriteImage(GTV2_dose, subj_dir2+'Masked_Dose_GTV.nii')
    sitk.WriteImage(CTV2_dose, subj_dir2+'Masked_Dose_CTV.nii')
    sitk.WriteImage(PTV2_dose, subj_dir2+'Masked_Dose_PTV.nii')
    
    sitk.WriteImage(GTV_sub, data_supradir+current+'/Dose_Diff_GTV.nii')
    sitk.WriteImage(CTV_sub, data_supradir+current+'/Dose_Diff_CTV.nii')
    sitk.WriteImage(PTV_sub, data_supradir+current+'/Dose_Diff_PTV.nii')