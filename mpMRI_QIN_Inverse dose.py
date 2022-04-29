# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 16:13:06 2020

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
    
    print('Starting analysis for '+current)
    
    nii_dir = subj_dir +'/nii/' 
   
    #Set paths to subfolders    
    orig_imgs = nii_dir +'Original images'
    reg_imgs = nii_dir +'Registered images'
    rel_imgs = nii_dir +'Relative images'
    bet_imgs = nii_dir +'BET images'
    smooth_imgs = nii_dir +'Smoothed images'
    norm_imgs = nii_dir +'Normalised images'
    prob_imgs = nii_dir +'Probability images'
 
#%% Generate DPBN maps by using the following formula D_(p,i)= 〖D_min+(〖D_max-D〗_min)×〖TP〗_i〗^n 
    
    # os.mkdir(nii_dir +'Dose images')
    dose_imgs = nii_dir +'Dose images'
   
    ADC_rBV_DPBN = sitk.ReadImage(dose_imgs +'/ADC_rBV_DPBN.nii')
    ADC_rBF_DPBN = sitk.ReadImage(dose_imgs +'/ADC_rBF_DPBN.nii')
    
    #Generate inverse dose prescription as Dinv = Dmax - Dpresc, where Dmax = 80 Gy
    Inv_ADC_rBV_DPBN = 80-ADC_rBV_DPBN
    Inv_ADC_rBF_DPBN = 80-ADC_rBF_DPBN
    
    sitk.WriteImage(Inv_ADC_rBV_DPBN, dose_imgs +'/Inv_ADC_rBV_DPBN.nii')
    sitk.WriteImage(Inv_ADC_rBF_DPBN, dose_imgs +'/Inv_ADC_rBF_DPBN.nii')
    
    #Generate mask of inverse dose in GTV
    GTV = sitk.ReadImage(rel_imgs +'/GTV.nii') #read GTV 
    Inv_ADC_rBV_DPBN_inGTV = generate_mask(Inv_ADC_rBV_DPBN, GTV)
    Inv_ADC_rBF_DPBN_inGTV = generate_mask(Inv_ADC_rBF_DPBN, GTV)
    
    sitk.WriteImage(Inv_ADC_rBV_DPBN_inGTV, dose_imgs +'/Inv_ADC_rBV_DPBN_inGTV.nii')
    sitk.WriteImage(Inv_ADC_rBF_DPBN_inGTV, dose_imgs +'/Inv_ADC_rBF_DPBN_inGTV.nii')
    
