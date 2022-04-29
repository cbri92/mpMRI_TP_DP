# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 14:33:01 2021

@author: cbri3325
"""


from ConvertDicomsSeries_ToNiiFile import *
from pathlib import Path
import SimpleITK as sitk
import os.path


data_supradir='path to supradirectory containing individual patients directories' #Set path to working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir = data_supradir+current
    
    subj_name = current
    
    print('Starting DICOM series to Nifti file conversion for '+current)
    
   
#%% Set paths to DICOM series folders

    for f in os.scandir(subj_dir):
        # if f.path.find('rBV') != -1:
        #     Inv_ADC_rBV_DP_path = f.path
        if f.path.find('rBF') != -1:
            Inv_ADC_rBF_DP_path = f.path

    # convert_dicomSeries_to_dicom_RTdosefile(Inv_ADC_rBV_DP_path, out_dir, "Inv_ADC_rBV_DP.dcm")
    convert_dicomSeries_to_Niifile(Inv_ADC_rBF_DP_path, subj_dir, "Inv_ADC_rBF_DP.nii")
