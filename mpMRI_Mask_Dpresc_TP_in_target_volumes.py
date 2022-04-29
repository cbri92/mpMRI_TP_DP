# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:01:24 2022

@author: cbri3325
"""


#%% Import functions 

import SimpleITK as sitk
import numpy as np
import os
import glob
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing DVH files' #Set path to DVH files directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

Dpresc_dir = 'path to supradirectory containing individual patients directories'

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current
    
    subj_Dpresc_dir = Dpresc_dir+current+'/Dose grids DPBN76/TP1/'
    
    struct_dir = subj_dir+'/RT structuresTP1'
    masked_dir = subj_dir+'/Dose masked to targets'
  
    print('Reading images for '+current) 
    
    for filename in glob.glob(subj_Dpresc_dir+'*inverse_dose_grid.nii.gz'):
        gunzip_shutil(filename, filename[:-3])
    
    for f in os.scandir(subj_Dpresc_dir):
        if f.path.find('inverse_dose') != -1:
            Inv_path = f.path
        
    Inv_dose = sitk.ReadImage(Inv_path) #Read inverse dose prescription image
    
    #Find path to GTV, CTV and PTV
    for f in os.scandir(struct_dir):
        if f.path.find('GTV') != -1:
            GTV_path = f.path
        elif f.path.find('CTV') != -1:
            CTV_path = f.path
        elif f.path.find('PTV') != -1:
            PTV_path = f.path
            
    #Read GTV, CTV and PTV
    GTV = sitk.ReadImage(GTV_path)
    CTV = sitk.ReadImage(CTV_path)
    PTV = sitk.ReadImage(PTV_path)
    
    #Make RT struct as binary labels
    GTV = GTV>1
    CTV = CTV>1
    PTV = PTV>1

    #Set same origin for targets and dose grids
    GTV.SetOrigin(Inv_dose.GetOrigin())
    CTV.SetOrigin(Inv_dose.GetOrigin())
    PTV.SetOrigin(Inv_dose.GetOrigin())
    PTV_margins = PTV-CTV
    
    #%%Generate dose prescriptions and TP from inverse dose prescription
    DP_Presc_dose = 80-Inv_dose
    TP = (DP_Presc_dose-60)/20
    TP = generate_mask(TP, CTV) 

    print('Generating dose prescription and tumour probability masked to targets for '+current)   
    
    #Generate DP prescription dose masked to GTV, CTV and PTV
    DP_GTV_Presc_dose = generate_mask(DP_Presc_dose, GTV)
    DP_CTV_Presc_dose = generate_mask(DP_Presc_dose, CTV)
    DP_PTV_Presc_dose = set_mask_value(DP_CTV_Presc_dose, PTV_margins, 60)
    DP_PTV_Presc_dose = generate_mask(DP_PTV_Presc_dose, PTV)

    #Generate STD prescription dose masked to GTV, CTV and PTV
    STD_Presc_dose = set_mask_value(Inv_dose, PTV, 60)
    
    STD_GTV_Presc_dose = generate_mask(STD_Presc_dose, GTV)
    STD_CTV_Presc_dose = generate_mask(STD_Presc_dose, CTV)
    STD_PTV_Presc_dose = generate_mask(STD_Presc_dose, PTV)
    
    #Generate TP masked to GTV, CTV and PTV    
    GTV_TP = generate_mask(TP, GTV)
    CTV_TP = generate_mask(TP, CTV)
    PTV_TP = generate_mask(TP, PTV)

    
    #%%Save images
    
    #Save masked TP
    # os.mkdir(masked_dir+'/Tumour probability masked')
    TP_masked_dir = masked_dir+'/Tumour probability masked'
    
    sitk.WriteImage(GTV_TP, TP_masked_dir+'/TP_GTV.nii')
    sitk.WriteImage(CTV_TP, TP_masked_dir+'/TP_CTV.nii')
    sitk.WriteImage(PTV_TP, TP_masked_dir+'/TP_PTV.nii')
    
    #Save masked DP prescriptions    
    DP_masked_dir = masked_dir+'/Dose painting masked plans'
    
    sitk.WriteImage(DP_GTV_Presc_dose, DP_masked_dir+'/DP_Presc_GTV.nii')
    sitk.WriteImage(DP_CTV_Presc_dose, DP_masked_dir+'/DP_Presc_CTV.nii')
    sitk.WriteImage(DP_PTV_Presc_dose, DP_masked_dir+'/DP_Presc_PTV.nii')

    #Save masked STD prescriptions
    STD_masked_dir = masked_dir+'/Standard masked plans'
    
    sitk.WriteImage(STD_GTV_Presc_dose, STD_masked_dir+'/STD_Presc_GTV.nii')
    sitk.WriteImage(STD_CTV_Presc_dose, STD_masked_dir+'/STD_Presc_CTV.nii')
    sitk.WriteImage(STD_PTV_Presc_dose, STD_masked_dir+'/STD_Presc_PTV.nii')