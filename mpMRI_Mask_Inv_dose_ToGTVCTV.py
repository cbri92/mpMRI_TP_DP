# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 14:33:01 2021

@author: cbri3325
"""


from ImageAnalysisFunctions import *
from ConvertNii_ToDoseFiles import *
from pathlib import Path
import SimpleITK as sitk
import os.path
import glob


data_supradir='path to supradirectory containing individual patients directories' #Set path to working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir = data_supradir+current
    
    subj_name = current
    
    print('Starting masking Inverse dose into GTV-CTV margins for '+current)
    
#%% Unzip all nii.gz files and delete original gzipped files

    for filename in glob.glob(subj_dir +'/' +'*.gz'):
        gunzip_shutil(filename, filename[:-3])
        
#%% Set paths to DICOM series folders

    for f in os.scandir(subj_dir):
        if f.path.find('DP') != -1:
            Inv_ADC_rBF_DP_path = f.path
        # elif f.path.find('TP') != -1:
        #     ADC_rBF_TP_path = f.path
        elif f.path.find('GTV') != -1:
            GTV_path = f.path
        elif f.path.find('CTV') != -1:
            CTV_path = f.path
        elif f.path.find('PTV') != -1:
            PTV_path = f.path

#%% Read images

    Inv_dose = sitk.ReadImage(Inv_ADC_rBF_DP_path)
    # TP = sitk.ReadImage(ADC_rBF_TP_path)
    GTV = sitk.ReadImage(GTV_path)
    CTV = sitk.ReadImage(CTV_path)
    PTV = sitk.ReadImage(PTV_path)
    
#%% Generate Mergin volume

    GTV = GTV>1
    sitk.WriteImage(GTV, subj_dir+'/GTV_new.nii')
    CTV = CTV>1
    sitk.WriteImage(CTV, subj_dir+'/CTV_new.nii')
    PTV = PTV>1
    sitk.WriteImage(PTV, subj_dir+'/PTV_new.nii')

    Margin_volume = CTV-GTV
    Margin_volume.SetOrigin(Inv_dose.GetOrigin())
    sitk.WriteImage(Margin_volume, subj_dir+'/RT_margins.nii')

#%% Generate masked inverse dose

    Masked_Inv_dose = generate_mask(Inv_dose, Margin_volume)
    sitk.WriteImage(Masked_Inv_dose, subj_dir+'/Masked_Inv_dose.nii')
    
#%% Convert nii dose image into DICOMS RT Dose file

    for f in os.scandir(subj_dir):
        if f.path.find('Aligned') != -1:
            dcm_ref = f.path
            
    convert_nii_to_dicom_RTdosefile(Masked_Inv_dose, dcm_ref, output_directory=subj_dir, out_filename="Masked_Inv_dose.dcm")
