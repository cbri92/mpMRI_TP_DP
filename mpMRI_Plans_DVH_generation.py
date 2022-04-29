# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 14:33:01 2021

@author: cbri3325
"""


from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *
from pathlib import Path
import SimpleITK as sitk
import os.path
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings


data_supradir='path to supradirectory containing individual patients directories' #Set path to working directory
if not os.path.isdir(data_supradir):
    warnings.warn('invalid data_supradir supplied; quitting', stacklevel=2)
    sys.exit(1)


subjs_path = [f.path for f in os.scandir(data_supradir) if f.is_dir()] #Create a list of the paths to the subjects directories
subjs_name = [f.name for f in os.scandir(data_supradir) if f.is_dir()] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir = data_supradir+current
    
    subj_name = current
    
    print('Starting DVH generation for ' + current)
    
    #%% Set path to dose and structure folders
    
    dose_dir = subj_dir + '/Dose grids DPBN76/TP1'
    struct_dir = subj_dir + '/RT structuresTP1'
    
    #%% Unzip all nii.gz files and delete original gzipped files

    for filename in glob.glob(dose_dir +'/' +'*.gz'):
        gunzip_shutil(filename, filename[:-3])
        
    for filename in glob.glob(struct_dir +'/' +'*.gz'):
        gunzip_shutil(filename, filename[:-3])        

#%%Delete files you don't need

    for filename in glob.glob(dose_dir +'/' +'*inverse*'):
        os.remove(filename)
        
    for filename in glob.glob(dose_dir +'/' +'*summed*'):
        os.remove(filename)
        
    for filename in glob.glob(struct_dir +'/' +'*ADC*'):
        os.remove(filename)     

 # %%Rename dose files such that to remove prefix
    for filename in os.listdir(dose_dir):
        os.chdir(dose_dir)
        os.rename(filename, filename[13:])
        
#%%Define path to dose file
    for f in os.scandir(dose_dir):
        if f.path.find('DPBN76_planned') != -1:
            DP_dose_file_path = f.path
    DP_dose = sitk.ReadImage(DP_dose_file_path)
    
    for f in os.scandir(dose_dir):
        if f.path.find('STANDARD') != -1:
            STD_dose_file_path = f.path
    STD_dose = sitk.ReadImage(STD_dose_file_path)
    
#%%Resample STD plan to 1.2 isotropic resolution of DP plan
    STD_dose = Resample_image(STD_dose, DP_dose)
    sitk.WriteImage(STD_dose, dose_dir+'/STANDARD_plan_dose_grid.nii')
    
    # %%Rename structure files such that to remove prefix
    for filename in os.listdir(struct_dir):
        os.chdir(struct_dir)
        os.rename(filename, filename[13:])
    
    #%%Generate DVH folder
    if os.path.isdir(subj_dir + '/DVH_files'):
        warnings.warn(f'emptying subj_dir directory {subj_dir}', stacklevel=2)
        shutil.rmtree(subj_dir + '/DVH_files')
    [dum, write_dir_name_STD] = os.path.split(STD_dose_file_path)
    write_dir_name_STD = write_dir_name_STD.replace(' ', '_')
    write_dir_name_STD = write_dir_name_STD.replace('.nii', '')
    
    os.mkdir(subj_dir + '/DVH_files')
    os.mkdir(subj_dir + '/DVH_files' + '/' + write_dir_name_STD)
    STD_DVH_dir = subj_dir + '/DVH_files'  + '/' + write_dir_name_STD
    
    [dum, write_dir_name_DP] = os.path.split(DP_dose_file_path)
    write_dir_name_DP = write_dir_name_DP.replace(' ', '_')
    write_dir_name_DP = write_dir_name_DP.replace('.nii', '')
    os.mkdir(subj_dir + '/DVH_files' + '/' + write_dir_name_DP)
    DP_DVH_dir = subj_dir + '/DVH_files'  + '/' + write_dir_name_DP
    
    #%%Generate DVH excel file for each structure
    
    for f in os.scandir(struct_dir):
        file_write_name = f.name.replace(' ','_') # replcaes spaces with _
        [file_write_name, dum] = os.path.splitext(file_write_name)  # removes old extension
        file_write_name = file_write_name + '.xlsx'  # adds new extension

        struct = sitk.ReadImage(f.path)
        struct = struct > 0
        x = allVoxInt(STD_dose, struct) #This function calculates a 2D flat array of the dose for each voxel within the 3D structure
        counts, bin_edges = np.histogram(x, bins=80, range=(0, x.max()), normed=None, weights=None, density=False)
        histcum = 100*(1 - np.cumsum(counts)/len(x)) #cumulative histogram values: y axis
        
        #Save DVH to excel file
        Histo_results = pd.ExcelWriter(STD_DVH_dir +'/' + file_write_name)
        df_dict = {'Dose [Gy]': bin_edges[:-1], 'Relative Volume [%]': histcum}
        Histo_df = pd.DataFrame(data=df_dict)
        Histo_df.to_excel(Histo_results, header=None, index=False)
        Histo_results.save()
        
        y = allVoxInt(DP_dose, struct) #This function calculates a 2D flat array of the dose for each voxel within the 3D structure
        counts, bin_edges = np.histogram(y, bins=80, range=(0, y.max()), normed=None, weights=None, density=False)
        histcum = 100*(1 - np.cumsum(counts)/len(y)) #cumulative histogram values: y axis
        
        #Save DVH to excel file
        Histo_results = pd.ExcelWriter(DP_DVH_dir +'/' + file_write_name)
        df_dict = {'Dose [Gy]': bin_edges[:-1], 'Relative Volume [%]': histcum}
        Histo_df = pd.DataFrame(data=df_dict)
        Histo_df.to_excel(Histo_results, header=None, index=False)
        Histo_results.save()

    #%%Mask plans in RT structures
    
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
    GTV.SetOrigin(DP_dose.GetOrigin())
    CTV.SetOrigin(DP_dose.GetOrigin())
    PTV.SetOrigin(DP_dose.GetOrigin())
    
    #Generate masked targets
    
    os.mkdir(subj_dir + '/Dose masked to targets/')
    masked_dir = subj_dir + '/Dose masked to targets/'
    
    os.mkdir(masked_dir +'Standard masked plans')
    STD_masked_dir = masked_dir +'Standard masked plans'
    
    os.mkdir(masked_dir +'Dose painting masked plans')
    DP_masked_dir = masked_dir +'Dose painting masked plans'
    
    print('Generating dose masked to targets for '+current)   
    
    GTV_DP_dose = generate_mask(DP_dose, GTV)
    CTV_DP_dose = generate_mask(DP_dose, CTV)
    PTV_DP_dose = generate_mask(DP_dose, PTV)
    
    GTV_STD_dose = generate_mask(STD_dose, GTV)
    CTV_STD_dose = generate_mask(STD_dose, CTV)
    PTV_STD_dose = generate_mask(STD_dose, PTV)
    
    #%%Save images
    
    sitk.WriteImage(GTV_DP_dose, DP_masked_dir+'/DP_Dose_GTV.nii')
    sitk.WriteImage(CTV_DP_dose, DP_masked_dir+'/DP_Dose_CTV.nii')
    sitk.WriteImage(PTV_DP_dose, DP_masked_dir+'/DP_Dose_PTV.nii')
    
    sitk.WriteImage(GTV_STD_dose, STD_masked_dir+'/STD_Dose_GTV.nii')
    sitk.WriteImage(CTV_STD_dose, STD_masked_dir+'/STD_Dose_CTV.nii')
    sitk.WriteImage(PTV_STD_dose, STD_masked_dir+'/STD_Dose_PTV.nii')