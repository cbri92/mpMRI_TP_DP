# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 12:54:42 2020

@author: cbri3325
"""

#%% Import functions 

import SimpleITK as sitk
import os
import glob
import shutil
from ImageAnalysisFunctions import *

import sys
sys.path.append("C:/Users/cbri3325/Anaconda3/lib/site-packages/platipy/__init__.py") # Path containing PlatiPy library
from platipy.dicom.nifti_to_series.convert import convert_nifti_to_dicom_series

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing individual patients directories' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names
subjs_name.remove('DICOMrefFLAIR')
subjs_name.remove('DICOMrefT1CE')

n_subj = len(subjs_name) #Total number of subjects

dicom_refFLAIR = data_supradir+'DICOMrefFLAIR/'
dicom_refT1CE = data_supradir+'DICOMrefT1CE/'

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    subj_dir = data_supradir+current
    subj_name = current
   
    # nii_dir = subj_dir +'/nii/' 
   
    #Set paths to subfolders    
    # orig_imgs = nii_dir +'Original images'
    # bet_imgs = nii_dir +'BET images'

 

    #Read T1CE, FLAIR and T1 bet images
    # T1CE = sitk.ReadImage(orig_imgs+'/t1_axial_post.nii', sitk.sitkFloat32) #read T1CE bet image
    # FLAIR = sitk.ReadImage(orig_imgs+'/flair.nii', sitk.sitkFloat32) #read FLAIR bet image
    # # T1CE = sitk.ReadImage(bet_imgs+'/T1CE_bet.nii', sitk.sitkFloat32) #read T1 bet image
    # FLAIR = sitk.ReadImage(bet_imgs+'/FLAIR_bet.nii', sitk.sitkFloat32) #read flair bet image
    T1CE = sitk.ReadImage(subj_dir+'/T1CE.nii', sitk.sitkFloat32) #read T1CE image
    FLAIR = sitk.ReadImage(subj_dir+'/FLAIR.nii', sitk.sitkFloat32) #read FLAIR image
    
    # T1new = sitk.Cast(T1, sitk.sitkUInt16)
    # sitk.WriteImage(T1new, bet_imgs+'/T1new.nii')
    
    T1CEnew = sitk.Cast(T1CE, sitk.sitkUInt16)
    # sitk.WriteImage(T1CEnew, bet_imgs+'/T1CEnew.nii')
    
    FLAIRnew = sitk.Cast(FLAIR, sitk.sitkUInt16)
    # sitk.WriteImage(FLAIRnew, bet_imgs+'/FLAIRnew.nii')
    
    #Make DICOM directory
    # os.mkdir(orig_imgs +'/DICOMS')
    # DICOM_dir = orig_imgs +'/DICOMS'
    # os.mkdir(bet_imgs +'/DICOMS')
    # DICOM_dir = bet_imgs +'/DICOMS'
    os.mkdir(subj_dir +'/DICOMS')
    DICOM_dir = subj_dir +'/DICOMS'
    
    # os.mkdir(DICOM_dir +'/T1')
    # T1_dir = DICOM_dir +'/T1'
    
    os.mkdir(DICOM_dir +'/T1CE')
    T1CE_dir = DICOM_dir +'/T1CE'
    
    os.mkdir(DICOM_dir +'/FLAIR')
    FLAIR_dir = DICOM_dir +'/FLAIR'
    
    # convert_nifti_to_dicom_series(T1new, dicom_ref, output_directory=T1_dir)
    convert_nifti_to_dicom_series(T1CEnew, dicom_refT1CE, output_directory=T1CE_dir)
    convert_nifti_to_dicom_series(FLAIRnew, dicom_refFLAIR, output_directory=FLAIR_dir)
