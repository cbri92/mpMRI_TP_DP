# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 14:52:43 2021

@author: Caterina Brighi

This script calculates ICC (2,1) in a ROI between images acquired at two timepoints. 
According to McGraw and Wong (1996) Convention, ICC (2,1) quantifies the variability 
between voxels relative to the measurement error and is calculated 
with a two-way random effects, absolute agreement, single rater/measurement model.
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
        
data_supradir = 'Path to data directory' #Set working directory

tp1_dir = data_supradir+'Timepoint1/' #Set path to directory containing images from Timepoint 1
tp2_dir = data_supradir+'Timepoint2/' #Set path to directory containing images from Timepoint 2

os.mkdir(data_supradir + 'ICC analysis at two timepoints/')
results_dir = data_supradir + 'ICC analysis at two timepoints/'

subjs_path = [ f.path for f in os.scandir(tp2_dir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(tp2_dir) if f.is_dir() ] #Create a list of subjects names

n_subj = len(subjs_name) #Total number of subjects

Results_ICC = {'ADC': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'rBV': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'rBF': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'ADC-rBV TP': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'ADC-rBF TP': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'ADC-rBV DP': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%']), 'ADC-rBF DP': pd.DataFrame(columns=['Subject_ID','Type', 'Description','ICC', 'F', 'df1', 'df2', 'pval', 'CI95%'])}
IntraWriter = pd.ExcelWriter(data_supradir +'ICCResults.xlsx', engine='xlsxwriter')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:
    
    subj_dir1 = tp1_dir+current
    subj_dir2 = tp2_dir+current
    
    subj_name = current
    
    print('Starting ICC analysis for '+current)
    
    # os.mkdir(results_dir+current)
    # out_dir = results_dir+current
    
    # # #Create an excel spreadsheet with the stats results
    # ICC_input_table = pd.ExcelWriter(out_dir +'/ICC_input_table.xlsx', engine='xlsxwriter')
   
#%% Set paths to subfolders    
   
    bet_imgs1 = subj_dir1 +'/nii/BET images'
    rel_imgs1 = subj_dir1 +'/nii/Relative images'
    prob_imgs1 = subj_dir1 +'/nii/Probability images'
    dose_imgs1 = subj_dir1 +'/nii/Dose images'
    
    bet_imgs2 = subj_dir2 +'/nii/BET images'
    rel_imgs2 = subj_dir2 +'/nii/Relative images'
    prob_imgs2 = subj_dir2 +'/nii/Probability images'
    dose_imgs2 = subj_dir2 +'/nii/Dose images'
    
#%%Read images from two timepoints
    
    #From Timepoint 1
    ADC_1 = sitk.ReadImage(bet_imgs1 +'/ADC_bet.nii')
    rBV_1 = sitk.ReadImage(rel_imgs1 +'/rBV_rel.nii')
    rBF_1 = sitk.ReadImage(rel_imgs1 +'/rBF_rel.nii')
    ADC_rBV_TP1 = sitk.ReadImage(prob_imgs1 +'/ADC_rBV_prob.nii')
    ADC_rBF_TP1 = sitk.ReadImage(prob_imgs1 +'/ADC_rBF_prob.nii')
    ADC_rBV_DPBN1 = sitk.ReadImage(dose_imgs1 +'/ADC_rBV_DP_full.nii')
    ADC_rBF_DPBN1 = sitk.ReadImage(dose_imgs1 +'/ADC_rBF_DP_full.nii')
    
    
    #From Timepoint 2
    ADC_2 = sitk.ReadImage(bet_imgs2 +'/ADC_bet.nii')
    rBV_2 = sitk.ReadImage(rel_imgs2 +'/rBV_rel.nii')
    rBF_2 = sitk.ReadImage(rel_imgs2 +'/rBF_rel.nii')
    ADC_rBV_TP2 = sitk.ReadImage(prob_imgs2 +'/ADC_rBV_prob.nii')
    ADC_rBF_TP2 = sitk.ReadImage(prob_imgs2 +'/ADC_rBF_prob.nii')
    ADC_rBV_DPBN2 = sitk.ReadImage(dose_imgs2 +'/ADC_rBV_DP_full.nii')
    ADC_rBF_DPBN2 = sitk.ReadImage(dose_imgs2 +'/ADC_rBF_DP_full.nii')
    
    #Labels
    Tumour_mask = sitk.ReadImage(rel_imgs1 +'/Margin_volume.nii')
    
#%%Extract mean image intensity in ROI
    
    print('Extracting voxel values...')
    #From Timepoint 1
    ADC_1_vxls = allVoxInt(ADC_1, Tumour_mask)
    rBV_1_vxls = allVoxInt(rBV_1, Tumour_mask)
    rBF_1_vxls = allVoxInt(rBF_1, Tumour_mask)
    ADC_rBV_TP1_vxls = allVoxInt(ADC_rBV_TP1, Tumour_mask)
    ADC_rBF_TP1_vxls = allVoxInt(ADC_rBF_TP1, Tumour_mask)
    ADC_rBV_DPBN1_vxls = allVoxInt(ADC_rBV_DPBN1, Tumour_mask)
    ADC_rBF_DPBN1_vxls = allVoxInt(ADC_rBF_DPBN1, Tumour_mask)
    
    #From Timepoint 2
    ADC_2_vxls = allVoxInt(ADC_2, Tumour_mask)
    rBV_2_vxls = allVoxInt(rBV_2, Tumour_mask)
    rBF_2_vxls = allVoxInt(rBF_2, Tumour_mask)
    ADC_rBV_TP2_vxls = allVoxInt(ADC_rBV_TP2, Tumour_mask)
    ADC_rBF_TP2_vxls = allVoxInt(ADC_rBF_TP2, Tumour_mask)
    ADC_rBV_DPBN2_vxls = allVoxInt(ADC_rBV_DPBN2, Tumour_mask)
    ADC_rBF_DPBN2_vxls = allVoxInt(ADC_rBF_DPBN2, Tumour_mask)    
   
#%%Create empty dataframes to populate as going through the loop
    
    ADC_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'ADC'], dtype=float)
    rBV_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'rBV'], dtype=float)
    rBF_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'rBF'], dtype=float)
    ADCrBV_TP_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'ADC-rBV TP'], dtype=float)
    ADCrBF_TP_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'ADC-rBF TP'], dtype=float)
    ADCrBV_DP_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'ADC-rBV DP'], dtype=float)
    ADCrBF_DP_input_df = pd.DataFrame(columns=['Voxel_ID', 'Measurement', 'ADC-rBF DP'], dtype=float)
    
#%% Create a patient voxels intensity input dataframe and append to general

    print('Creating individual dataframes...')    

    ADC_d1 = {'Voxel_ID':[i for i in range(len(ADC_1_vxls))], 'Measurement':[1 for i in range(len(ADC_1_vxls))], 'ADC': [i[0] for i in ADC_1_vxls.tolist()]}
    ADC_df1 = pd.DataFrame(data=ADC_d1)
    
    ADC_d2 = {'Voxel_ID':[i for i in range(len(ADC_2_vxls))], 'Measurement':[2 for i in range(len(ADC_2_vxls))], 'ADC': [i[0] for i in ADC_2_vxls.tolist()]}
    ADC_df2 = pd.DataFrame(data=ADC_d2)

    ADC_input_df=ADC_input_df.append(ADC_df1)
    ADC_input_df=ADC_input_df.append(ADC_df2)
    
    
    rBV_d1 = {'Voxel_ID':[i for i in range(len(rBV_1_vxls))], 'Measurement':[1 for i in range(len(rBV_1_vxls))], 'rBV': [i[0] for i in rBV_1_vxls.tolist()]}
    rBV_df1 = pd.DataFrame(data=rBV_d1)
    
    rBV_d2 = {'Voxel_ID':[i for i in range(len(rBV_2_vxls))], 'Measurement':[2 for i in range(len(rBV_2_vxls))], 'rBV': [i[0] for i in rBV_2_vxls.tolist()]}
    rBV_df2 = pd.DataFrame(data=rBV_d2)

    rBV_input_df=rBV_input_df.append(rBV_df1)
    rBV_input_df=rBV_input_df.append(rBV_df2)
    
    
    rBF_d1 = {'Voxel_ID':[i for i in range(len(rBF_1_vxls))], 'Measurement':[1 for i in range(len(rBF_1_vxls))], 'rBF': [i[0] for i in rBF_1_vxls.tolist()]}
    rBF_df1 = pd.DataFrame(data=rBF_d1)
    
    rBF_d2 = {'Voxel_ID':[i for i in range(len(rBF_2_vxls))], 'Measurement':[2 for i in range(len(rBF_2_vxls))], 'rBF': [i[0] for i in rBF_2_vxls.tolist()]}
    rBF_df2 = pd.DataFrame(data=rBF_d2)

    rBF_input_df=rBF_input_df.append(rBF_df1)
    rBF_input_df=rBF_input_df.append(rBF_df2)


    ADCrBV_TP_d1 = {'Voxel_ID':[i for i in range(len(ADC_rBV_TP1_vxls))], 'Measurement':[1 for i in range(len(ADC_rBV_TP1_vxls))], 'ADC-rBV TP': [i[0] for i in ADC_rBV_TP1_vxls.tolist()]}
    ADCrBV_TP_df1 = pd.DataFrame(data=ADCrBV_TP_d1)
    
    ADCrBV_TP_d2 = {'Voxel_ID':[i for i in range(len(ADC_rBV_TP2_vxls))], 'Measurement':[2 for i in range(len(ADC_rBV_TP2_vxls))], 'ADC-rBV TP': [i[0] for i in ADC_rBV_TP2_vxls.tolist()]}
    ADCrBV_TP_df2 = pd.DataFrame(data=ADCrBV_TP_d2)

    ADCrBV_TP_input_df=ADCrBV_TP_input_df.append(ADCrBV_TP_df1)
    ADCrBV_TP_input_df=ADCrBV_TP_input_df.append(ADCrBV_TP_df2)


    ADCrBF_TP_d1 = {'Voxel_ID':[i for i in range(len(ADC_rBF_TP1_vxls))], 'Measurement':[1 for i in range(len(ADC_rBF_TP1_vxls))], 'ADC-rBF TP': [i[0] for i in ADC_rBF_TP1_vxls.tolist()]}
    ADCrBF_TP_df1 = pd.DataFrame(data=ADCrBF_TP_d1)
    
    ADCrBF_TP_d2 = {'Voxel_ID':[i for i in range(len(ADC_rBF_TP2_vxls))], 'Measurement':[2 for i in range(len(ADC_rBF_TP2_vxls))], 'ADC-rBF TP': [i[0] for i in ADC_rBF_TP2_vxls.tolist()]}
    ADCrBF_TP_df2 = pd.DataFrame(data=ADCrBF_TP_d2)

    ADCrBF_TP_input_df=ADCrBF_TP_input_df.append(ADCrBF_TP_df1)
    ADCrBF_TP_input_df=ADCrBF_TP_input_df.append(ADCrBF_TP_df2)


    ADCrBV_DP_d1 = {'Voxel_ID':[i for i in range(len(ADC_rBV_DP1_vxls))], 'Measurement':[1 for i in range(len(ADC_rBV_DP1_vxls))], 'ADC-rBV DP': [i[0] for i in ADC_rBV_DP1_vxls.tolist()]}
    ADCrBV_DP_df1 = pd.DataFrame(data=ADCrBV_DP_d1)
    
    ADCrBV_DP_d2 = {'Voxel_ID':[i for i in range(len(ADC_rBV_DP2_vxls))], 'Measurement':[2 for i in range(len(ADC_rBV_DP2_vxls))], 'ADC-rBV DP': [i[0] for i in ADC_rBV_DP2_vxls.tolist()]}
    ADCrBV_DP_df2 = pd.DataFrame(data=ADCrBV_DP_d2)

    ADCrBV_DP_input_df=ADCrBV_DP_input_df.append(ADCrBV_DP_df1)
    ADCrBV_DP_input_df=ADCrBV_DP_input_df.append(ADCrBV_DP_df2)


    ADCrBF_DP_d1 = {'Voxel_ID':[i for i in range(len(ADC_rBF_DP1_vxls))], 'Measurement':[1 for i in range(len(ADC_rBF_DP1_vxls))], 'ADC-rBF DP': [i[0] for i in ADC_rBF_DP1_vxls.tolist()]}
    ADCrBF_DP_df1 = pd.DataFrame(data=ADCrBF_DP_d1)
    
    ADCrBF_DP_d2 = {'Voxel_ID':[i for i in range(len(ADC_rBF_DP2_vxls))], 'Measurement':[2 for i in range(len(ADC_rBF_DP2_vxls))], 'ADC-rBF DP': [i[0] for i in ADC_rBF_DP2_vxls.tolist()]}
    ADCrBF_DP_df2 = pd.DataFrame(data=ADCrBF_DP_d2)

    ADCrBF_DP_input_df=ADCrBF_DP_input_df.append(ADCrBF_DP_df1)
    ADCrBF_DP_input_df=ADCrBF_DP_input_df.append(ADCrBF_DP_df2)
        

                 
#%% Calculate ROI-based group repeatibility statistics


    #Calculate ICC
 
    print('Calculating ADC ICC...')
    
    ICC_ADC = pg.intraclass_corr(data=ADC_input_df, targets='Voxel_ID', raters='Measurement', ratings='ADC').round(3)  
    Results_ICC['ADC'] = Results_ICC['ADC'].append({'Subject_ID': current,
                                                          'Type':ICC_ADC.loc[1, 'Type'],
                                                          'Description':ICC_ADC.loc[1, 'Description'], 
                                                          'ICC':ICC_ADC.loc[1, 'ICC'], 
                                                          'F':ICC_ADC.loc[1, 'F'],
                                                          'df1':ICC_ADC.loc[1, 'df1'],
                                                          'df2':ICC_ADC.loc[1, 'df2'],
                                                          'pval':ICC_ADC.loc[1, 'pval'],
                                                          'CI95%':ICC_ADC.loc[1, 'CI95%']
                                                          }, ignore_index=True)
    
    print('Calculating rBV ICC...')
    
    ICC_rBV = pg.intraclass_corr(data=rBV_input_df, targets='Voxel_ID', raters='Measurement', ratings='rBV').round(3)  
    Results_ICC['rBV'] = Results_ICC['rBV'].append({'Subject_ID': current,
                                                          'Type':ICC_rBV.loc[1, 'Type'],
                                                          'Description':ICC_rBV.loc[1, 'Description'], 
                                                          'ICC':ICC_rBV.loc[1, 'ICC'], 
                                                          'F':ICC_rBV.loc[1, 'F'],
                                                          'df1':ICC_rBV.loc[1, 'df1'],
                                                          'df2':ICC_rBV.loc[1, 'df2'],
                                                          'pval':ICC_rBV.loc[1, 'pval'],
                                                          'CI95%':ICC_rBV.loc[1, 'CI95%']
                                                          }, ignore_index=True)
    
    print('Calculating rBF ICC...')
    
    ICC_rBF = pg.intraclass_corr(data=rBF_input_df, targets='Voxel_ID', raters='Measurement', ratings='rBF').round(3)  
    Results_ICC['rBF'] = Results_ICC['rBF'].append({'Subject_ID': current,
                                                          'Type':ICC_rBF.loc[1, 'Type'],
                                                          'Description':ICC_rBF.loc[1, 'Description'], 
                                                          'ICC':ICC_rBF.loc[1, 'ICC'], 
                                                          'F':ICC_rBF.loc[1, 'F'],
                                                          'df1':ICC_rBF.loc[1, 'df1'],
                                                          'df2':ICC_rBF.loc[1, 'df2'],
                                                          'pval':ICC_rBF.loc[1, 'pval'],
                                                          'CI95%':ICC_rBF.loc[1, 'CI95%']
                                                          }, ignore_index=True)
       
 
    print('Calculating ADC-rBV TP ICC...')
    
    ICC_ADCrBV_TP = pg.intraclass_corr(data=ADCrBV_TP_input_df, targets='Voxel_ID', raters='Measurement', ratings='ADC-rBV TP').round(3)  
    Results_ICC['ADC-rBV TP'] = Results_ICC['ADC-rBV TP'].append({'Subject_ID': current,
                                                          'Type':ICC_ADCrBV_TP.loc[1, 'Type'],
                                                          'Description':ICC_ADCrBV_TP.loc[1, 'Description'], 
                                                          'ICC':ICC_ADCrBV_TP.loc[1, 'ICC'], 
                                                          'F':ICC_ADCrBV_TP.loc[1, 'F'],
                                                          'df1':ICC_ADCrBV_TP.loc[1, 'df1'],
                                                          'df2':ICC_ADCrBV_TP.loc[1, 'df2'],
                                                          'pval':ICC_ADCrBV_TP.loc[1, 'pval'],
                                                          'CI95%':ICC_ADCrBV_TP.loc[1, 'CI95%']
                                                          }, ignore_index=True)
    
    print('Calculating ADC-rBF TP ICC...')
    
    ICC_ADCrBF_TP = pg.intraclass_corr(data=ADCrBF_TP_input_df, targets='Voxel_ID', raters='Measurement', ratings='ADC-rBF TP').round(3)  
    Results_ICC['ADC-rBF TP'] = Results_ICC['ADC-rBF TP'].append({'Subject_ID': current,
                                                          'Type':ICC_ADCrBF_TP.loc[1, 'Type'],
                                                          'Description':ICC_ADCrBF_TP.loc[1, 'Description'], 
                                                          'ICC':ICC_ADCrBF_TP.loc[1, 'ICC'], 
                                                          'F':ICC_ADCrBF_TP.loc[1, 'F'],
                                                          'df1':ICC_ADCrBF_TP.loc[1, 'df1'],
                                                          'df2':ICC_ADCrBF_TP.loc[1, 'df2'],
                                                          'pval':ICC_ADCrBF_TP.loc[1, 'pval'],
                                                          'CI95%':ICC_ADCrBF_TP.loc[1, 'CI95%']
                                                          }, ignore_index=True)
    
    print('Calculating ADC-rBV DP ICC...')
    
    ICC_ADCrBV_DP = pg.intraclass_corr(data=ADCrBV_DP_input_df, targets='Voxel_ID', raters='Measurement', ratings='ADC-rBV DP').round(3)  
    Results_ICC['ADC-rBV DP'] = Results_ICC['ADC-rBV DP'].append({'Subject_ID': current,
                                                          'Type':ICC_ADCrBV_DP.loc[1, 'Type'],
                                                          'Description':ICC_ADCrBV_DP.loc[1, 'Description'], 
                                                          'ICC':ICC_ADCrBV_DP.loc[1, 'ICC'], 
                                                          'F':ICC_ADCrBV_DP.loc[1, 'F'],
                                                          'df1':ICC_ADCrBV_DP.loc[1, 'df1'],
                                                          'df2':ICC_ADCrBV_DP.loc[1, 'df2'],
                                                          'pval':ICC_ADCrBV_DP.loc[1, 'pval'],
                                                          'CI95%':ICC_ADCrBV_DP.loc[1, 'CI95%']
                                                          }, ignore_index=True)
    
    print('Calculating ADC-rBF DP ICC...')
    
    ICC_ADCrBF_DP = pg.intraclass_corr(data=ADCrBF_DP_input_df, targets='Voxel_ID', raters='Measurement', ratings='ADC-rBF DP').round(3)  
    Results_ICC['ADC-rBF DP'] = Results_ICC['ADC-rBF DP'].append({'Subject_ID': current,
                                                          'Type':ICC_ADCrBF_DP.loc[1, 'Type'],
                                                          'Description':ICC_ADCrBF_DP.loc[1, 'Description'], 
                                                          'ICC':ICC_ADCrBF_DP.loc[1, 'ICC'], 
                                                          'F':ICC_ADCrBF_DP.loc[1, 'F'],
                                                          'df1':ICC_ADCrBF_DP.loc[1, 'df1'],
                                                          'df2':ICC_ADCrBF_DP.loc[1, 'df2'],
                                                          'pval':ICC_ADCrBF_DP.loc[1, 'pval'],
                                                          'CI95%':ICC_ADCrBF_DP.loc[1, 'CI95%']
                                                          }, ignore_index=True)
 
#%%Save all dataframes to excel files here
print('Save all results to excel files')

for name, df in Results_ICC.items():
    df.to_excel(IntraWriter, sheet_name=name, index=False)
IntraWriter.save()



