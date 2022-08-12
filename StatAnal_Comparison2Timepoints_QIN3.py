# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:52:43 2021

@author: Caterina Brighi

This script calculates different repeatability metrics, including the following:
Voxel-wise statistics: 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'.
ROI-wise statistics: 'WMS', 'BMS','wSD', 'bSD', 'tSD', 'RC', 'RC Upper', 'RC Lower', 'wCV', 'ICC'
It also plots scatter density histograms of the voxels intensity distribution within a ROI,
Bolt Altman plots and intensity correlation plotsbetween two timepoints, both for voxel-wise analysis and ROI-wise analysis.
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


#%% Set Working directory
        
data_supradir = 'Path to data directory' #Set working directory

tp1_dir = data_supradir+'Timepoint1/' #Set path to directory containing images from Timepoint 1
tp2_dir = data_supradir+'Timepoint2/' #Set path to directory containing images from Timepoint 2

os.mkdir(data_supradir + 'Intensity comparison at two timepoints_new/')
results_dir = data_supradir + 'Intensity comparison at two timepoints_new/'

subjs_path = [ f.path for f in os.scandir(tp2_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(tp2_dir) if f.is_dir() ] #Create a list of subjects names
# subjs_name.remove('Visualization Results Images')

n_subj = len(subjs_name) #Total number of subjects

#Create an excel spreadsheet with the stats results
Group_stats_results = pd.ExcelWriter(results_dir +'Group_stats_results.xlsx')

#Create empty dataframes to populate as going through the loop

ADCdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
rBVdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
rBFdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
ADC_rBV_TPdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
ADC_rBF_TPdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
ADC_rBV_DPBNdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])
ADC_rBF_DPBNdf = pd.DataFrame(columns=['Subject_ID', 'Mean intensity in ROI at tp1', 'Mean intensity in ROI at tp2', 'Between voxels variance', 'Within voxels variance', 'ROI ICC', 'Within voxel CoV', 'Repeatability coefficient', 'RC Upper', 'RC Lower'])

Repeatibility_statsdf = pd.DataFrame(columns=['Image type', 'WMS', 'BMS','wSD', 'bSD', 'tSD', 'RC', 'RC Upper', 'RC Lower', 'wCV', 'ICC'])

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir1 = tp1_dir+current
    subj_dir2 = tp2_dir+current
    
    subj_name = current
    
    print('Starting feature analysis for '+current)
    
    os.mkdir(results_dir+current)
    out_dir = results_dir+current
   
#%% Set paths to subfolders    
   
    bet_imgs1 = subj_dir1 +'/nii/BET images'
    rel_imgs1 = subj_dir1 +'/nii/Relative images'
    norm_imgs1 = subj_dir1 +'/nii/Normalised images'
    prob_imgs1 = subj_dir1 +'/nii/Probability images'
    dose_imgs1 = subj_dir1 +'/nii/Dose images'
    
    bet_imgs2 = subj_dir2 +'/nii/BET images'
    rel_imgs2 = subj_dir2 +'/nii/Relative images'
    norm_imgs2 = subj_dir2 +'/nii/Normalised images'
    prob_imgs2 = subj_dir2 +'/nii/Probability images'
    dose_imgs2 = subj_dir2 +'/nii/Dose images'
    
#%%Read images from two timepoints
    
    #From Timepoint 1
    ADC_1 = sitk.ReadImage(bet_imgs1 +'/ADC_bet.nii')
    rBV_1 = sitk.ReadImage(rel_imgs1 +'/rBV_rel.nii')
    rBF_1 = sitk.ReadImage(rel_imgs1 +'/rBF_rel.nii')
    ADC_rBV_TP1 = sitk.ReadImage(prob_imgs1 +'/ADC_rBV_def.nii')
    ADC_rBF_TP1 = sitk.ReadImage(prob_imgs1 +'/ADC_rBF_def.nii')
    ADC_rBV_DPBN1 = sitk.ReadImage(dose_imgs1 +'/ADC_rBV_DPBN.nii')
    ADC_rBF_DPBN1 = sitk.ReadImage(dose_imgs1 +'/ADC_rBF_DPBN.nii')
    
    
    #From Timepoint 2
    ADC_2 = sitk.ReadImage(bet_imgs2 +'/ADC_bet.nii')
    rBV_2 = sitk.ReadImage(rel_imgs2 +'/rBV_rel.nii')
    rBF_2 = sitk.ReadImage(rel_imgs2 +'/rBF_rel.nii')
    ADC_rBV_TP2 = sitk.ReadImage(prob_imgs2 +'/ADC_rBV_def.nii')
    ADC_rBF_TP2 = sitk.ReadImage(prob_imgs2 +'/ADC_rBF_def.nii')
    ADC_rBV_DPBN2 = sitk.ReadImage(dose_imgs2 +'/ADC_rBV_DPBN.nii')
    ADC_rBF_DPBN2 = sitk.ReadImage(dose_imgs2 +'/ADC_rBF_DPBN.nii')
    
    #Labels
    Tumour_mask = sitk.ReadImage(rel_imgs1 +'/Margin_volume.nii')
    CTRL_VOI = sitk.ReadImage(rel_imgs1 +'/CTRL_VOI.nii')
    # GTV = sitk.ReadImage(rel_imgs1 +'/GTV.nii')
    
#%%Extract mean image intensity in ROI
    
    #From Timepoint 1
    
    ADC_1_vxls = allVoxInt(ADC_1, Tumour_mask)
    rBV_1_vxls = allVoxInt(rBV_1, Tumour_mask)
    rBF_1_vxls = allVoxInt(rBF_1, Tumour_mask)
    # ADC_rBV_TP1_vxls = allVoxInt(ADC_rBV_TP1, Tumour_mask)
    # ADC_rBF_TP1_vxls = allVoxInt(ADC_rBF_TP1, Tumour_mask)
    # ADC_rBV_DPBN1_vxls = allVoxInt(ADC_rBV_DPBN1, Tumour_mask)
    # ADC_rBF_DPBN1_vxls = allVoxInt(ADC_rBF_DPBN1, Tumour_mask)
    ADC_rBV_TP1_vxls = allVoxInt_threshold(ADC_rBV_TP1, Tumour_mask, 0)
    ADC_rBF_TP1_vxls = allVoxInt_threshold(ADC_rBF_TP1, Tumour_mask, 0)
    ADC_rBV_DPBN1_vxls = allVoxInt_threshold(ADC_rBV_DPBN1, Tumour_mask, 60)
    ADC_rBF_DPBN1_vxls = allVoxInt_threshold(ADC_rBF_DPBN1, Tumour_mask, 60)
    
      
    ADC_1_Mean = ADC_1_vxls.mean()
    rBV_1_Mean = rBV_1_vxls.mean()
    rBF_1_Mean = rBF_1_vxls.mean()
    ADC_rBV_TP1_Mean = ADC_rBV_TP1_vxls.mean()
    ADC_rBF_TP1_Mean = ADC_rBF_TP1_vxls.mean()
    ADC_rBV_DPBN1_Mean = ADC_rBV_DPBN1_vxls.mean()
    ADC_rBF_DPBN1_Mean = ADC_rBF_DPBN1_vxls.mean()
    
     
    #From Timepoint 2
    ADC_2_vxls = allVoxInt(ADC_2, Tumour_mask)
    rBV_2_vxls = allVoxInt(rBV_2, Tumour_mask)
    rBF_2_vxls = allVoxInt(rBF_2, Tumour_mask)
    # ADC_rBV_TP2_vxls = allVoxInt(ADC_rBV_TP2, Tumour_mask)
    # ADC_rBF_TP2_vxls = allVoxInt(ADC_rBF_TP2, Tumour_mask)
    # ADC_rBV_DPBN2_vxls = allVoxInt(ADC_rBV_DPBN2, Tumour_mask)
    # ADC_rBF_DPBN2_vxls = allVoxInt(ADC_rBF_DPBN2, Tumour_mask)
    ADC_rBV_TP2_vxls = allVoxInt_threshold(ADC_rBV_TP2, Tumour_mask, 0)
    ADC_rBF_TP2_vxls = allVoxInt_threshold(ADC_rBF_TP2, Tumour_mask, 0)
    ADC_rBV_DPBN2_vxls = allVoxInt_threshold(ADC_rBV_DPBN2, Tumour_mask, 60)
    ADC_rBF_DPBN2_vxls = allVoxInt_threshold(ADC_rBF_DPBN2, Tumour_mask, 60)
    
    
    ADC_2_Mean = ADC_2_vxls.mean()
    rBV_2_Mean = rBV_2_vxls.mean()
    rBF_2_Mean = rBF_2_vxls.mean()
    ADC_rBV_TP2_Mean = ADC_rBV_TP2_vxls.mean()
    ADC_rBF_TP2_Mean = ADC_rBF_TP2_vxls.mean()
    ADC_rBV_DPBN2_Mean = ADC_rBV_DPBN2_vxls.mean()
    ADC_rBF_DPBN2_Mean = ADC_rBF_DPBN2_vxls.mean()
    

#%% Plot scatter density histogram
    
    h=sns.jointplot(x=ADC_rBV_TP1_vxls, y=ADC_rBV_DPBN1_vxls, kind='scatter', s=5, color='g')
    h.set_axis_labels('Tumour probability', 'Dose prescription [Gy]', fontsize=16)
    plt.tight_layout()
    plt.savefig(out_dir+'/DPTPscatter_ADC_rBV_TP1.png')     
    plt.close()
    
    h=sns.jointplot(x=ADC_rBV_TP2_vxls, y=ADC_rBV_DPBN2_vxls, kind='scatter', s=5, color='g')
    h.set_axis_labels('Tumour probability', 'Dose prescription [Gy]', fontsize=16)
    plt.tight_layout()
    plt.savefig(out_dir+'/DPTPscatter_ADC_rBV_TP2.png')     
    plt.close()

    h=sns.jointplot(x=ADC_rBF_TP1_vxls, y=ADC_rBF_DPBN1_vxls, kind='scatter', s=5, color='r')
    h.set_axis_labels('Tumour probability', 'Dose prescription [Gy]', fontsize=16)
    plt.tight_layout()
    plt.savefig(out_dir+'/DPTPscatter_ADC_rBF_TP1.png')     
    plt.close()
    
    h=sns.jointplot(x=ADC_rBF_TP2_vxls, y=ADC_rBF_DPBN2_vxls, kind='scatter', s=5, color='r')
    h.set_axis_labels('Tumour probability', 'Dose prescription [Gy]', fontsize=16)
    plt.tight_layout()
    plt.savefig(out_dir+'/DPTPscatter_ADC_rBF_TP2.png')     
    plt.close()
    
#%% Calculate ROI-based group repeatibility statistics

    #Calculate ICC within ROI

    ADC_ICC = Repeatibility_metrics(2, ADC_1_vxls, ADC_2_vxls)
    rBV_ICC = Repeatibility_metrics(2, rBV_1_vxls, rBV_2_vxls)
    rBF_ICC = Repeatibility_metrics(2, rBF_1_vxls, rBF_2_vxls)
    ADC_rBV_TP_ICC = Repeatibility_metrics(2, ADC_rBV_TP1_vxls, ADC_rBV_TP2_vxls)
    ADC_rBF_TP_ICC = Repeatibility_metrics(2, ADC_rBF_TP1_vxls, ADC_rBF_TP2_vxls)
    ADC_rBV_DPBN_ICC = Repeatibility_metrics(2, ADC_rBV_DPBN1_vxls, ADC_rBV_DPBN2_vxls)
    ADC_rBF_DPBN_ICC = Repeatibility_metrics(2, ADC_rBF_DPBN1_vxls, ADC_rBF_DPBN2_vxls)
    
    #Plot Bland-Altman plots
    Bland_Altman_plt(ADC_1_vxls, ADC_2_vxls, 'r', 5, 'ADC, GTV-CTV margin volume', out_dir+'/ADC_BAplot.png')
    Bland_Altman_plt(rBV_1_vxls, rBV_2_vxls, 'b', 5, 'rBV, GTV-CTV margin volume', out_dir+'/rBV_BAplot.png')
    Bland_Altman_plt(rBF_1_vxls, rBF_2_vxls, 'g', 5, 'rBF, GTV-CTV margin volume', out_dir+'/rBF_BAplot.png')    
    Bland_Altman_plt(ADC_rBV_TP1_vxls, ADC_rBV_TP2_vxls, 'y', 5, 'ADC_rBV TP, GTV-CTV margin volume', out_dir+'/ADC_rBV_TP_BAplot.png')
    Bland_Altman_plt(ADC_rBF_TP1_vxls, ADC_rBF_TP2_vxls, 'b', 5, 'ADC_rBF TP, GTV-CTV margin volume', out_dir+'/ADC_rBF_TP_BAplot.png')
    Bland_Altman_plt(ADC_rBV_DPBN1_vxls, ADC_rBV_DPBN2_vxls, 'r', 5, 'ADC_rBV DPBN, GTV-CTV margin volume', out_dir+'/ADC_rBV_DPBN_BAplot.png')
    Bland_Altman_plt(ADC_rBF_DPBN1_vxls, ADC_rBF_DPBN2_vxls, 'g', 5, 'ADC_rBF DPBN, GTV-CTV margin volume', out_dir+'/ADC_rBF_DPBN_BAplot.png')    
    
    #Correlation plots
    Vox_Corrleation_plt(ADC_1_vxls, ADC_2_vxls, 'r', 5, 'ADC, GTV-CTV margin volume', out_dir+'/ADC_Corrplot.png')
    Vox_Corrleation_plt(rBV_1_vxls, rBV_2_vxls, 'b', 5, 'rBV, GTV-CTV margin volume', out_dir+'/rBV_Corrplot.png')
    Vox_Corrleation_plt(rBF_1_vxls, rBF_2_vxls, 'g', 5, 'rBF, GTV-CTV margin volume', out_dir+'/rBF_Corrplot.png')    
    Vox_Corrleation_plt(ADC_rBV_TP1_vxls, ADC_rBV_TP2_vxls, 'y', 5, 'ADC_rBV TP, GTV-CTV margin volume', out_dir+'/ADC_rBV_TP_Corrplot.png')
    Vox_Corrleation_plt(ADC_rBF_TP1_vxls, ADC_rBF_TP2_vxls, 'b', 5, 'ADC_rBF TP, GTV-CTV margin volume', out_dir+'/ADC_rBF_TP_Corrplot.png')
    Vox_Corrleation_plt(ADC_rBV_DPBN1_vxls, ADC_rBV_DPBN2_vxls, 'r', 5, 'ADC_rBV DPBN, GTV-CTV margin volume', out_dir+'/ADC_rBV_DPBN_Corrplot.png')
    Vox_Corrleation_plt(ADC_rBF_DPBN1_vxls, ADC_rBF_DPBN2_vxls, 'g', 5, 'ADC_rBF DPBN, GTV-CTV margin volume', out_dir+'/ADC_rBF_DPBN_Corrplot.png')    
    
    
#%% Export subjects stats into excel spreadsheet

    ADC_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':ADC_1_Mean, 'Mean intensity in ROI at tp2':ADC_2_Mean, 'Between voxels variance':ADC_ICC.get('Between subject variance'), 'Within voxels variance':ADC_ICC.get('Within subject variance'), 'ROI ICC':ADC_ICC.get('ICC'), 'Within voxel CoV':ADC_ICC.get('wCV'), 'Repeatability coefficient':ADC_ICC.get('RC'), 'RC Upper':ADC_ICC.get('RC Upper'), 'RC Lower':ADC_ICC.get('RC Lower')}
    rBV_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':rBV_1_Mean, 'Mean intensity in ROI at tp2':rBV_2_Mean, 'Between voxels variance':rBV_ICC.get('Between subject variance'), 'Within voxels variance':rBV_ICC.get('Within subject variance'), 'ROI ICC':rBV_ICC.get('ICC'), 'Within voxel CoV':rBV_ICC.get('wCV'), 'Repeatability coefficient':rBV_ICC.get('RC'), 'RC Upper':rBV_ICC.get('RC Upper'), 'RC Lower':rBV_ICC.get('RC Lower')}
    rBF_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':rBF_1_Mean, 'Mean intensity in ROI at tp2':rBF_2_Mean, 'Between voxels variance':rBF_ICC.get('Between subject variance'), 'Within voxels variance':rBF_ICC.get('Within subject variance'), 'ROI ICC':rBF_ICC.get('ICC'), 'Within voxel CoV':rBF_ICC.get('wCV'), 'Repeatability coefficient':rBF_ICC.get('RC'), 'RC Upper':rBF_ICC.get('RC Upper'), 'RC Lower':rBF_ICC.get('RC Lower')}
    ADC_rBV_TP_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':ADC_rBV_TP1_Mean, 'Mean intensity in ROI at tp2':ADC_rBV_TP2_Mean, 'Between voxels variance':ADC_rBV_TP_ICC.get('Between subject variance'), 'Within voxels variance':ADC_rBV_TP_ICC.get('Within subject variance'), 'ROI ICC':ADC_rBV_TP_ICC.get('ICC'), 'Within voxel CoV':ADC_rBV_TP_ICC.get('wCV'), 'Repeatability coefficient':ADC_rBV_TP_ICC.get('RC'), 'RC Upper':ADC_rBV_TP_ICC.get('RC Upper'), 'RC Lower':ADC_rBV_TP_ICC.get('RC Lower')}
    ADC_rBF_TP_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':ADC_rBF_TP1_Mean, 'Mean intensity in ROI at tp2':ADC_rBF_TP2_Mean, 'Between voxels variance':ADC_rBF_TP_ICC.get('Between subject variance'), 'Within voxels variance':ADC_rBF_TP_ICC.get('Within subject variance'), 'ROI ICC':ADC_rBF_TP_ICC.get('ICC'), 'Within voxel CoV':ADC_rBF_TP_ICC.get('wCV'), 'Repeatability coefficient':ADC_rBF_TP_ICC.get('RC'), 'RC Upper':ADC_rBF_TP_ICC.get('RC Upper'), 'RC Lower':ADC_rBF_TP_ICC.get('RC Lower')}
    ADC_rBV_DPBN_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':ADC_rBV_DPBN1_Mean, 'Mean intensity in ROI at tp2':ADC_rBV_DPBN2_Mean, 'Between voxels variance':ADC_rBV_DPBN_ICC.get('Between subject variance'), 'Within voxels variance':ADC_rBV_DPBN_ICC.get('Within subject variance'), 'ROI ICC':ADC_rBV_DPBN_ICC.get('ICC'), 'Within voxel CoV':ADC_rBV_DPBN_ICC.get('wCV'), 'Repeatability coefficient':ADC_rBV_DPBN_ICC.get('RC'), 'RC Upper':ADC_rBV_DPBN_ICC.get('RC Upper'), 'RC Lower':ADC_rBV_DPBN_ICC.get('RC Lower')}
    ADC_rBF_DPBN_stats = {'Subject_ID':current, 'Mean intensity in ROI at tp1':ADC_rBF_DPBN1_Mean, 'Mean intensity in ROI at tp2':ADC_rBF_DPBN2_Mean, 'Between voxels variance':ADC_rBF_DPBN_ICC.get('Between subject variance'), 'Within voxels variance':ADC_rBF_DPBN_ICC.get('Within subject variance'), 'ROI ICC':ADC_rBF_DPBN_ICC.get('ICC'), 'Within voxel CoV':ADC_rBF_DPBN_ICC.get('wCV'), 'Repeatability coefficient':ADC_rBF_DPBN_ICC.get('RC'), 'RC Upper':ADC_rBF_DPBN_ICC.get('RC Upper'), 'RC Lower':ADC_rBF_DPBN_ICC.get('RC Lower')}

    ADCdf = ADCdf.append(ADC_stats, ignore_index=True)
    rBVdf = rBVdf.append(rBV_stats, ignore_index=True)
    rBFdf = rBFdf.append(rBF_stats, ignore_index=True)
    ADC_rBV_TPdf = ADC_rBV_TPdf.append(ADC_rBV_TP_stats, ignore_index=True)
    ADC_rBF_TPdf = ADC_rBF_TPdf.append(ADC_rBF_TP_stats, ignore_index=True)
    ADC_rBV_DPBNdf = ADC_rBV_DPBNdf.append(ADC_rBV_DPBN_stats, ignore_index=True)
    ADC_rBF_DPBNdf = ADC_rBF_DPBNdf.append(ADC_rBF_DPBN_stats, ignore_index=True)
    
    
    
#%%Calculate group arrays

ADC_meanROI1 = ADCdf['Mean intensity in ROI at tp1'].to_numpy()
ADC_meanROI2 = ADCdf['Mean intensity in ROI at tp2'].to_numpy()

rBV_meanROI1 = rBVdf['Mean intensity in ROI at tp1'].to_numpy()
rBV_meanROI2 = rBVdf['Mean intensity in ROI at tp2'].to_numpy()

rBF_meanROI1 = rBFdf['Mean intensity in ROI at tp1'].to_numpy()
rBF_meanROI2 = rBFdf['Mean intensity in ROI at tp2'].to_numpy()

ADC_rBV_TP_meanROI1 = ADC_rBV_TPdf['Mean intensity in ROI at tp1'].to_numpy()
ADC_rBV_TP_meanROI2 = ADC_rBV_TPdf['Mean intensity in ROI at tp2'].to_numpy()

ADC_rBF_TP_meanROI1 = ADC_rBF_TPdf['Mean intensity in ROI at tp1'].to_numpy()
ADC_rBF_TP_meanROI2 = ADC_rBF_TPdf['Mean intensity in ROI at tp2'].to_numpy()

ADC_rBV_DPBN_meanROI1 = ADC_rBV_DPBNdf['Mean intensity in ROI at tp1'].to_numpy()
ADC_rBV_DPBN_meanROI2 = ADC_rBV_DPBNdf['Mean intensity in ROI at tp2'].to_numpy()

ADC_rBF_DPBN_meanROI1 = ADC_rBF_DPBNdf['Mean intensity in ROI at tp1'].to_numpy()
ADC_rBF_DPBN_meanROI2 = ADC_rBF_DPBNdf['Mean intensity in ROI at tp2'].to_numpy()


#%%Calculate rebeatability metrics

ADC_ROIrepeatibility=Repeatibility_metrics(2, ADC_meanROI1, ADC_meanROI2)
rBV_ROIrepeatibility=Repeatibility_metrics(2, rBV_meanROI1, rBV_meanROI2)
rBF_ROIrepeatibility=Repeatibility_metrics(2, rBF_meanROI1, rBF_meanROI2)
ADC_rBV_TP_ROIrepeatibility=Repeatibility_metrics(2, ADC_rBV_TP_meanROI1, ADC_rBV_TP_meanROI2)
ADC_rBF_TP_ROIrepeatibility=Repeatibility_metrics(2, ADC_rBF_TP_meanROI1, ADC_rBF_TP_meanROI2)
ADC_rBV_DPBN_ROIrepeatibility=Repeatibility_metrics(2, ADC_rBV_DPBN_meanROI1, ADC_rBV_DPBN_meanROI2)
ADC_rBF_DPBN_ROIrepeatibility=Repeatibility_metrics(2, ADC_rBF_DPBN_meanROI1, ADC_rBF_DPBN_meanROI2)


#Add image modality to dictionary
ADC_ROIrepeatibility['Image type']= 'ADC ROI'
rBV_ROIrepeatibility['Image type']= 'rBV ROI'
rBF_ROIrepeatibility['Image type']= 'rBF ROI'
ADC_rBV_TP_ROIrepeatibility['Image type']= 'ADC_rBV_TP ROI'
ADC_rBF_TP_ROIrepeatibility['Image type']= 'ADC_rBF_TP ROI'
ADC_rBV_DPBN_ROIrepeatibility['Image type']= 'ADC_rBV_DPBN ROI'
ADC_rBF_DPBN_ROIrepeatibility['Image type']= 'ADC_rBF_DPBN ROI'


#Append metrics to repeatibility sheet
Repeatibility_statsdf = Repeatibility_statsdf.append(ADC_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(rBV_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(rBF_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(ADC_rBV_TP_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(ADC_rBF_TP_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(ADC_rBV_DPBN_ROIrepeatibility, ignore_index=True)

Repeatibility_statsdf = Repeatibility_statsdf.append(ADC_rBF_DPBN_ROIrepeatibility, ignore_index=True)


#%%Save all dataframes to excel files here
print('Save all results to excel files')

ADCdf.to_excel(Group_stats_results, sheet_name='ADC', index=False)
rBVdf.to_excel(Group_stats_results, sheet_name='rBV', index=False)
rBFdf.to_excel(Group_stats_results, sheet_name='rBF', index=False)
ADC_rBV_TPdf.to_excel(Group_stats_results, sheet_name='ADC_rBV_TP', index=False)
ADC_rBF_TPdf.to_excel(Group_stats_results, sheet_name='ADC_rBF_TP', index=False)
ADC_rBV_DPBNdf.to_excel(Group_stats_results, sheet_name='ADC_rBV_DPBN', index=False)
ADC_rBF_DPBNdf.to_excel(Group_stats_results, sheet_name='ADC_rBF_DPBN', index=False)
Repeatibility_statsdf.to_excel(Group_stats_results, sheet_name='Group Repeatibility metrics', index=False)

Group_stats_results.save()   

#%%Plot results of ICC

ADC_ICC = ADCdf['ROI ICC'].to_numpy()
rBV_ICC = rBVdf['ROI ICC'].to_numpy()
rBF_ICC = rBFdf['ROI ICC'].to_numpy()
ADC_rBV_TP_ICC = ADC_rBV_TPdf['ROI ICC'].to_numpy()
ADC_rBF_TP_ICC = ADC_rBF_TPdf['ROI ICC'].to_numpy()
ADC_rBV_DPBN_ICC = ADC_rBV_DPBNdf['ROI ICC'].to_numpy()
ADC_rBF_DPBN_ICC = ADC_rBF_DPBNdf['ROI ICC'].to_numpy()

ICC_data = [ADC_ICC, rBV_ICC, rBF_ICC, ADC_rBV_TP_ICC, ADC_rBF_TP_ICC, ADC_rBV_DPBN_ICC, ADC_rBF_DPBN_ICC]
ICC_labels = ['ADC', 'rBV', 'rBF', 'ADC_rBV_TP', 'ADC_rBF_TP', 'ADC_rBV_DPBN', 'ADC_rBF_DPBN']

fig, ax = plt.subplots(figsize=(10, 6))
ax.boxplot(ICC_data, labels=ICC_labels)
ax.set_ylim(0, 1)
plt.xlabel('Image modality', fontsize=20)
plt.ylabel('ICC', fontsize=20)
plt.setp(ax.get_xticklabels(), rotation = 45, fontsize=12)
plt.setp(ax.get_yticklabels(), fontsize=16)
fig.tight_layout(pad=2.0)
plt.savefig(results_dir+'ICC.png')     
plt.close()

#%%Plot results of group analysis

#Plot Bland-Altman plots
Bland_Altman_plt(ADC_meanROI1, ADC_meanROI2, 'r', 40, 'ADC, GTV-CTV margin volume', results_dir+'/ADC_BAplot.png')
Bland_Altman_plt(rBV_meanROI1, rBV_meanROI2, 'b', 40,'rBV, GTV-CTV margin volume', results_dir+'/rBV_BAplot.png')
Bland_Altman_plt(rBF_meanROI1, rBF_meanROI2, 'g', 40,'rBF, GTV-CTV margin volume', results_dir+'/rBF_BAplot.png')    
Bland_Altman_plt(ADC_rBV_TP_meanROI1, ADC_rBV_TP_meanROI2, 'y', 40,'ADC_rBV TP,GTV-CTV margin volume', results_dir+'/ADC_rBV_TP_BAplot.png')
Bland_Altman_plt(ADC_rBF_TP_meanROI1, ADC_rBF_TP_meanROI2, 'b', 40,'ADC_rBF TP, GTV-CTV margin volume', results_dir+'/ADC_rBF_TP_BAplot.png')
Bland_Altman_plt(ADC_rBV_DPBN_meanROI1, ADC_rBV_DPBN_meanROI2, 'r', 40,'ADC_rBV DPBN, GTV-CTV margin volume', results_dir+'/ADC_rBV_DPBN_BAplot.png')
Bland_Altman_plt(ADC_rBF_DPBN_meanROI1, ADC_rBF_DPBN_meanROI2, 'g', 40,'ADC_rBF DPBN, GTV-CTV margin volume', results_dir+'/ADC_rBF_DPBN_BAplot.png')    
    
#Correlation plots
Vox_Corrleation_plt(ADC_meanROI1, ADC_meanROI2, 'r', 40,'ADC, GTV-CTV margin volume', results_dir+'/ADC_Corrplot.png')
Vox_Corrleation_plt(rBV_meanROI1, rBV_meanROI2, 'b', 40,'rBV, GTV-CTV margin volume', results_dir+'/rBV_Corrplot.png')
Vox_Corrleation_plt(rBF_meanROI1, rBF_meanROI2, 'g', 40,'rBF, GTV-CTV margin volume', results_dir+'/rBF_Corrplot.png')    
Vox_Corrleation_plt(ADC_rBV_TP_meanROI1, ADC_rBV_TP_meanROI2, 'y', 40,'ADC_rBV TP, GTV-CTV margin volume', results_dir+'/ADC_rBV_TP_Corrplot.png')
Vox_Corrleation_plt(ADC_rBF_TP_meanROI1, ADC_rBF_TP_meanROI2, 'b', 40,'ADC_rBF TP, GTV-CTV margin volume', results_dir+'/ADC_rBF_TP_Corrplot.png')
Vox_Corrleation_plt(ADC_rBV_DPBN_meanROI1, ADC_rBV_DPBN_meanROI2, 'r', 40,'ADC_rBV DPBN, GTV-CTV margin volume', results_dir+'/ADC_rBV_DPBN_Corrplot.png')
Vox_Corrleation_plt(ADC_rBF_DPBN_meanROI1, ADC_rBF_DPBN_meanROI2, 'g', 40,'ADC_rBF DPBN, GTV-CTV margin volume', results_dir+'/ADC_rBF_DPBN_Corrplot.png')    
    

os.mkdir(results_dir + '/Results outputs/')
output_results = results_dir + '/Results outputs'

for file in glob.glob(results_dir +'/' +'*.png'):
    shutil.move(os.path.join(results_dir, file), output_results)
    
shutil.move(results_dir+'Group_stats_results.xlsx', output_results)
