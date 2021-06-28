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
        
data_supradir = 'C:/Users/cbri3325/Dropbox (Sydney Uni)/Caterina Brighi/Data/test/' #Set working directory

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
 
      
    # # #Find paths to T1CE, ADC, rBF and rBV images
    # T1CEpath = glob.glob(orig_imgs + '/' + '*t1_axial_post*')
    # T1CEpath = str(T1CEpath[0])
    # T1CEpath = os.path.abspath(T1CEpath)
    
    # FLAIRpath = glob.glob(orig_imgs + '/' + '*flair*')
    # FLAIRpath = str(FLAIRpath[0])
    # FLAIRpath = os.path.abspath(FLAIRpath)
    
    # ADCpath = glob.glob(orig_imgs + '/' + '*highres_adc*')
    # ADCpath = str(ADCpath[0])
    # ADCpath = os.path.abspath(ADCpath)
    
    # rBFpath = glob.glob(orig_imgs + '/' + '*rBF_ge_noAIF*')
    # rBFpath = str(rBFpath[0])
    # rBFpath = os.path.abspath(rBFpath)
    
    # rBVpath = glob.glob(orig_imgs + '/' + '*rBV_ge_noAIF*')
    # rBVpath = str(rBVpath[0])
    # rBVpath = os.path.abspath(rBVpath)
    
    # T1betpath = glob.glob(orig_imgs + '/' + '*t1_bet*')
    # T1betpath = str(T1betpath[0])
    # T1betpath = os.path.abspath(T1betpath)
    
    # DSCpath = glob.glob(orig_imgs + '/' + '*dsc_ge*')
    # DSCpath = str(DSCpath[0])
    # DSCpath = os.path.abspath(DSCpath)
    
    
    # #Read T1CE, FLAIR, ADC, rBF, rBV and T1 bet images
    # T1CE = sitk.ReadImage(T1CEpath, sitk.sitkFloat32) #read T1CE image
    # FLAIR = sitk.ReadImage(FLAIRpath, sitk.sitkFloat32) #read FLAIR image
    # ADC = sitk.ReadImage(ADCpath, sitk.sitkFloat32) #read ADC image
    # rBF = sitk.ReadImage(rBFpath, sitk.sitkFloat32) #read rBF image
    # rBV = sitk.ReadImage(rBVpath, sitk.sitkFloat32) #read rBV image
    # T1bet = sitk.ReadImage(T1betpath, sitk.sitkFloat32) #read T1 bet image
    # DSC = sitk.ReadImage(DSCpath, sitk.sitkFloat32) #read DSC image
 
    #Resample T1CE to isotropic 1.2 mm resolution
    # T1CE = Resample_image_to_specified_resolution(T1CE, 1.2, 1.2, 1.2)
    # sitk.WriteImage(T1CE, reg_imgs +'/T1CE_Resampled.nii')
    T1CE = sitk.ReadImage(reg_imgs +'/T1CE_Resampled.nii')
    
#%% Resample FLAIR, ADC, rCB and rCF to T1CE
  
    # fixed_image = T1CE
    # moving_images = [FLAIR, ADC, rBV, rBF]
    
    # for moving_image in moving_images:
    #     moving_resampled = Resample_image(moving_image, fixed_image)
        
    #     #Save image and transform to file
    #     if moving_image is FLAIR:
    #         sitk.WriteImage(moving_resampled, reg_imgs +'/'+ 'FLAIR_inT1CE.nii')
    #         FLAIR_inT1CE = sitk.ReadImage(reg_imgs +'/'+ 'FLAIR_inT1CE.nii')
    #         # sitk.WriteTransform(final_transform, reg_mtrx +'/'+ 'FLAIRF_2_T1CE.tfm')
    #     elif moving_image is ADC:
    #         sitk.WriteImage(moving_resampled, reg_imgs +'/'+ 'ADC_inT1CE.nii')
    #         ADC_inT1CE = sitk.ReadImage(reg_imgs +'/'+ 'ADC_inT1CE.nii')
    #         # sitk.WriteTransform(final_transform, reg_mtrx +'/'+ 'ADC_2_T1CE.tfm')
    #     elif moving_image is rBF:
    #         sitk.WriteImage(moving_resampled, reg_imgs +'/'+ 'rBF_inT1CE.nii')
    #         rBF_inT1CE = sitk.ReadImage(reg_imgs +'/'+ 'rBF_inT1CE.nii')
    #         # sitk.WriteTransform(final_transform, reg_mtrx +'/''rBF_2_T1CE.tfm')
    #     elif moving_image is rBV:
    #         sitk.WriteImage(moving_resampled, reg_imgs +'/'+ 'rBV_inT1CE.nii')
    #         rBV_inT1CE = sitk.ReadImage(reg_imgs +'/'+ 'rBV_inT1CE.nii')
    #         # sitk.WriteTransform(final_transform, reg_mtrx +'/'+ 'rBV_2_T1CE.tfm')   
    
   
#%% #%% Brain extraction
   
    brain_mask = sitk.ReadImage(bet_imgs +'/brain_mask.nii')
    # brain_mask = sitk.BinaryMedian(brain_mask, (2,2,2), 1.0, 0.0)
    # # brain_mask = Resample_roi(brain_mask, T1CE)
    # sitk.WriteImage(brain_mask, bet_imgs +'/brain_mask.nii')
    
    
    # check2 = input("Check and edit the brain_mask.nii in the Bet images folder. Once done editing type ENTER.")
    # print(check2)
    # brain_mask = sitk.ReadImage(bet_imgs +'/brain_mask.nii') 
    
    # # #Generate brain mask from dsc image 
    # DSC_mask = sitk.OtsuThreshold(rBV_inT1CE, 0, 1) #Generate binary mask of DSC_3D
    # DSC_mask = sitk.BinaryFillhole(DSC_mask)
    # DSC_mask = sitk.BinaryMorphologicalClosing(DSC_mask, (20,20,20), sitk.sitkBall, 1.0) #fill holes in initial VOI
    # DSC_mask = sitk.BinaryMorphologicalOpening(DSC_mask, (1,1,1), sitk.sitkBall, 0.0, 1.0) #remove small structures in initial VOI
    # DSC_mask = generate_mask(DSC_mask, brain_mask)
    # sitk.WriteImage(DSC_mask, bet_imgs +'/DSC_brain_mask.nii')
    # check1 = input("Check and edit the DSC_brain_mask.nii in the Bet images folder. Once done editing type ENTER.")
    # print(check1)
    DSC_mask_inT1CE = sitk.ReadImage(bet_imgs +'/DSC_brain_mask.nii')
    # DSC_mask_inT1CE = sitk.BinaryMedian(DSC_mask_inT1CE, (2,2,2), 1.0, 0.0)
    # # DSC_mask_inT1CE = Resample_roi(DSC_mask_inT1CE, T1CE)
    # sitk.WriteImage(DSC_mask_inT1CE, bet_imgs +'/DSC_brain_mask.nii')
 
       
    # #Brain extract FLAIR, ADC, rBV and rBF
    # T1bet = generate_mask(T1CE, brain_mask) #Refine brain extraction in T1CE
    # FLAIR_bet = generate_mask(FLAIR_inT1CE, brain_mask) #Brain extract FLAIR
    # ADC_bet = generate_mask(ADC_inT1CE, brain_mask) #Brain extract ADC
    # rBV_bet = generate_mask(rBV_inT1CE, brain_mask) #Brain extract rBV
    # rBF_bet = generate_mask(rBF_inT1CE, brain_mask) #Brain extract rBF

    # #Save brain extracted images
    # sitk.WriteImage(T1bet, bet_imgs +'/T1_bet.nii') #save new T1 bet
    # sitk.WriteImage(FLAIR_bet, bet_imgs +'/FLAIR_bet.nii') #save FLAIR bet
    # sitk.WriteImage(ADC_bet, bet_imgs +'/ADC_bet.nii') #save ADC bet
    # sitk.WriteImage(rBV_bet, bet_imgs +'/rBV_bet.nii') #save rBV bet
    # sitk.WriteImage(rBF_bet, bet_imgs +'/rBF_bet.nii') #save rBF bet
    
    #Read brain extracted images
    T1bet = sitk.ReadImage(bet_imgs +'/T1CE_bet.nii') #read new T1 bet
    FLAIR_bet = sitk.ReadImage(bet_imgs +'/FLAIR_bet.nii') #read FLAIR bet
    ADC_bet = sitk.ReadImage(bet_imgs +'/ADC_bet.nii') #read ADC bet
    rBV_bet = sitk.ReadImage(bet_imgs +'/rBV_bet.nii') #read rBV bet
    rBF_bet = sitk.ReadImage(bet_imgs +'/rBF_bet.nii') #read rBF bet


#%% Obtain relative images by divinding intensity by mean intensity of a 2D ROI in contralateral healthy brain
    

    #Read CTRL_ROI file
    # mkROI = input('Generate CTRL_VOI.nii on the T1CE_bet image and save it in the Relative images folder. Once done type ENTER.')
    # print(mkROI)
    CTRL_ROI = sitk.ReadImage(rel_imgs +'/CTRL_VOI.nii') #read CTRL ROI
    # CTRL_ROI = Resample_roi(CTRL_ROI, T1CE)
    # sitk.WriteImage(CTRL_ROI, rel_imgs +'/CTRL_VOI2.nii')
    

    #Calculate stats from CTRL_ROI applied to each image
    CTRL_ROI_FLAIRstats = getStatsRoi(CTRL_ROI, FLAIR_bet) #Calculate stats for CTRL_ROI applied onto FLAIR_bet
    CTRL_ROI_ADCstats = getStatsRoi(CTRL_ROI, ADC_bet) #Calculate stats for CTRL_ROI applied onto ADC_bet
    CTRL_ROI_rBVstats = getStatsRoi(CTRL_ROI, rBV_bet) #Calculate stats for CTRL_ROI applied onto rBV_bet
    CTRL_ROI_rBFstats = getStatsRoi(CTRL_ROI, rBF_bet) #Calculate stats for CTRL_ROI applied onto rBF_bet

    #Calculate mean intensity value from CTRL_ROI applied to each image
    CTRL_ROI_FLAIRmean = CTRL_ROI_FLAIRstats.get('Mean intensity [SUV]') #Calculate mean intensity value of CTRL_ROI applied onto FLAIR_bet
    CTRL_ROI_ADCmean = CTRL_ROI_ADCstats.get('Mean intensity [SUV]') #Calculate mean intensity value of CTRL_ROI applied onto ADC_bet
    CTRL_ROI_rBVmean = CTRL_ROI_rBVstats.get('Mean intensity [SUV]') #Calculate mean intensity value of CTRL_ROI applied onto rBV_bet
    CTRL_ROI_rBFmean = CTRL_ROI_rBFstats.get('Mean intensity [SUV]') #Calculate mean intensity value of CTRL_ROI applied onto rBF_bet
    
    #Generate relative images
    FLAIR_rel = FLAIR_bet/CTRL_ROI_FLAIRmean #Divide intensity in each voxel in FLAIR image by mean intensity from CTRL_ROI
    ADC_rel = ADC_bet/CTRL_ROI_ADCmean #Divide intensity in each voxel in ADC image by mean intensity from CTRL_ROI
    rBV_rel = rBV_bet/CTRL_ROI_rBVmean #Divide intensity in each voxel in rBV image by mean intensity from CTRL_ROI
    rBF_rel = rBF_bet/CTRL_ROI_rBFmean #Divide intensity in each voxel in rBF image by mean intensity from CTRL_ROI
    # rBV_rel = rBV_bet
    # rBF_rel = rBF_bet

    #Save the relative images
    sitk.WriteImage(FLAIR_rel, rel_imgs +'/FLAIR_rel.nii') #save FLAIR rel
    sitk.WriteImage(ADC_rel, rel_imgs +'/ADC_rel.nii') #save ADC rel
    sitk.WriteImage(rBV_rel, rel_imgs +'/rBV_rel.nii') #save rBV rel
    sitk.WriteImage(rBF_rel, rel_imgs +'/rBF_rel.nii') #save rBF rel

#%% Smooth relative images to 10x10x10 mm3 box

    # # Convert all voxels to mean of surrounding 1cm3 ROI
    # # FLAIR_rel = sitk.ReadImage(rel_imgs +'/FLAIR_rel.nii')
    # # ADC_rel = sitk.ReadImage(rel_imgs +'/ADC_rel.nii')
    # # rBV_rel = sitk.ReadImage(rel_imgs +'/rBV_rel.nii')
    # # rBF_rel = sitk.ReadImage(rel_imgs +'/rBF_rel.nii')
    
    # #Smooth over a radius of (23,23,2) voxels which based on images voxel size corresponded to 1cm3
    # # FLAIR_mean = sitk.Mean(FLAIR_rel, (23, 23, 2))
    # ADC_mean = sitk.Mean(ADC_rel, (23, 23, 1))
    # rBV_mean = sitk.Mean(rBV_rel, (23, 23, 1))
    # rBF_mean = sitk.Mean(rBF_rel, (23, 23, 1))
    
    #Smooth images preserving edges
    FLAIR_mean = sitk.Bilateral(FLAIR_rel, 2.0, 5.0, 100)
    ADC_mean = sitk.Bilateral(ADC_rel, 2.0, 5.0, 100)
    rBV_mean = sitk.Bilateral(rBV_rel, 2.0, 5.0, 100)
    rBF_mean = sitk.Bilateral(rBF_rel, 2.0, 5.0, 100)
    
    #Save smoothed images
    sitk.WriteImage(FLAIR_mean, smooth_imgs +'/FLAIR_mean.nii')
    sitk.WriteImage(ADC_mean, smooth_imgs +'/ADC_mean.nii')
    sitk.WriteImage(rBV_mean, smooth_imgs +'/rBV_mean.nii')
    sitk.WriteImage(rBF_mean, smooth_imgs +'/rBF_mean.nii')
    
    # # FLAIR_mean = sitk.ReadImage(smooth_imgs +'/FLAIR_mean.nii')
    # ADC_mean = sitk.ReadImage(smooth_imgs +'/ADC_mean.nii')
    # rBV_mean = sitk.ReadImage(smooth_imgs +'/rBV_mean.nii')
    # rBF_mean = sitk.ReadImage(smooth_imgs +'/rBF_mean.nii')
    
    
    

#%% Zero all voxels outside brain

    FLAIR_zero = generate_mask(FLAIR_mean, DSC_mask_inT1CE)
    ADC_zero = generate_mask(ADC_mean, DSC_mask_inT1CE)
    rBV_zero = generate_mask(rBV_mean, DSC_mask_inT1CE)
    rBF_zero = generate_mask(rBF_mean, DSC_mask_inT1CE)
    
    #Save zero images
    sitk.WriteImage(FLAIR_zero, norm_imgs +'/FLAIR_zero.nii')
    sitk.WriteImage(ADC_zero, norm_imgs +'/ADC_zero.nii')
    sitk.WriteImage(rBV_zero, norm_imgs +'/rBV_zero.nii')
    sitk.WriteImage(rBF_zero, norm_imgs +'/rBF_zero.nii')
    
#%% Normalize images (center and scale)

    #Convert sitk images into numpy arrays
    ADC_zero_nda = sitk.GetArrayFromImage(ADC_zero)
    rBV_zero_nda = sitk.GetArrayFromImage(rBV_zero)
    rBF_zero_nda = sitk.GetArrayFromImage(rBF_zero)

    #Calculate the mean of the zero image
    Mean_ADC = ADC_zero_nda[np.nonzero(ADC_zero_nda)].mean()
    Mean_rBV = rBV_zero_nda[np.nonzero(rBV_zero_nda)].mean()
    Mean_rBF = rBF_zero_nda[np.nonzero(rBF_zero_nda)].mean()
    
    #Generate centred image by subtracting the mean value from the intensity value of each voxel in the zero image
    ADC_center = ADC_zero - Mean_ADC
    rBV_center = rBV_zero - Mean_rBV
    rBF_center = rBF_zero - Mean_rBF

    ADC_center = generate_mask(ADC_center, DSC_mask_inT1CE)
    rBV_center = generate_mask(rBV_center, DSC_mask_inT1CE)
    rBF_center = generate_mask(rBF_center, DSC_mask_inT1CE)    
    
    #Convert sitk images into numpy arrays
    ADC_center_nda = sitk.GetArrayFromImage(ADC_center)
    rBV_center_nda = sitk.GetArrayFromImage(rBV_center)
    rBF_center_nda = sitk.GetArrayFromImage(rBF_center)

    #Calculate the standard deviation of the center image
    STD_ADC = ADC_center_nda[np.nonzero(ADC_center_nda)].std()
    STD_rBV = rBV_center_nda[np.nonzero(rBV_center_nda)].std()
    STD_rBF = rBF_center_nda[np.nonzero(rBF_center_nda)].std()    

    #Generate normalised image by dividing center image by the standard deviation value
    ADC_norm = ADC_center / STD_ADC
    rBV_norm = rBV_center / STD_rBV
    rBF_norm = rBF_center / STD_rBF
    
    ADC_norm = generate_mask(ADC_norm, DSC_mask_inT1CE)
    rBV_norm = generate_mask(rBV_norm, DSC_mask_inT1CE)
    rBF_norm = generate_mask(rBF_norm, DSC_mask_inT1CE)

    #Save the normalised images
    sitk.WriteImage(ADC_norm, norm_imgs +'/ADC_norm.nii')
    sitk.WriteImage(rBV_norm, norm_imgs +'/rBV_norm.nii')
    sitk.WriteImage(rBF_norm, norm_imgs +'/rBF_norm.nii')    
    
#%% Apply formula from regression analysis

    #Formula for combination of ADC and rBV: FORMULA = 1.6679 + 1.8564*ADC + 1.4792*DSC-CBV
    ADC_mul = sitk.Multiply(ADC_norm, 1.8564)
    rBV_mul = sitk.Multiply(rBV_norm, 1.4792)    
    # ADC_rBV = sitk.Add(ADC_mul, rBV_mul)
    # ADC_rBV = sitk.Add(ADC_rBV, 1.3966)
    ADC_rBV = 1.6679 + ADC_mul + rBV_mul #Combined image generated from combination of ADC and rBV images
    sitk.WriteImage(ADC_rBV, prob_imgs +'/ADC_rBV.nii') #Save ADC_rBV combination image    
    
    #Formula for combination of ADC and rBF: FORMULA = 1.6405 + 1.8038*ADC + 1.4332*DSC-CBF
    ADC_mul2 = sitk.Multiply(ADC_norm, 1.8038)
    rBF_mul = sitk.Multiply(rBF_norm, 1.4332)
    # ADC_rBF = sitk.Add(ADC_mul, rBF_mul)
    # ADC_rBF = sitk.Add(ADC_rBF, 1.5916)     
    ADC_rBF = 1.6405 + ADC_mul2 + rBF_mul #Combined image generated from combination of ADC and rBF images
    sitk.WriteImage(ADC_rBF, prob_imgs +'/ADC_rBF.nii') #Save ADC_rBF combination image


#%% Transform logg-odds in probability with formula: Tumour probability = exp(x)/(1+exp(x))

    #Generate tumour probability map for ADC_rBV
    ADC_rBV_exp = sitk.Exp(ADC_rBV)
    ADC_rBV_exp1 = sitk.Add(ADC_rBV_exp, 1)  
    ADC_rBV_prob = sitk.Divide(ADC_rBV_exp, ADC_rBV_exp1)
    sitk.WriteImage(ADC_rBV_prob, prob_imgs +'/ADC_rBV_prob.nii') #Save ADC_rBV probability map
    
    #Generate tumour probability map for ADC_rBF
    ADC_rBF_exp = sitk.Exp(ADC_rBF)
    ADC_rBF_exp1 = sitk.Add(ADC_rBF_exp, 1)  
    ADC_rBF_prob = sitk.Divide(ADC_rBF_exp, ADC_rBF_exp1)
    sitk.WriteImage(ADC_rBF_prob, prob_imgs +'/ADC_rBF_prob.nii') #Save ADC_rBF probability map
    
    
#%% Remove CSF and create defintive ADC/rBV and ADC/rBF probability maps

    #Create binary mask of CSF with FLAIR: threshold value might need adjustement according to voxel values
    # CSF_thr = input('Inspect FLAIR_zero.nii image and enter CSF threshold as a float: ')
    # CSF = FLAIR_zero > float(CSF_thr)
    # sitk.WriteImage(CSF, prob_imgs +'/CSF.nii') #Save CSF mask
    CSF = sitk.ReadImage(prob_imgs +'/CSF.nii')
    
    #Remove CSF from probability map
    ADC_rBV_probnocsf = generate_mask(ADC_rBV_prob, CSF)
    ADC_rBF_probnocsf = generate_mask(ADC_rBF_prob, CSF)

    #Remove voxel values outside brain due to removal of CSF and obtain definitive probability map
    ADC_rBV_def = generate_mask(ADC_rBV_probnocsf, DSC_mask_inT1CE)
    ADC_rBF_def = generate_mask(ADC_rBF_probnocsf, DSC_mask_inT1CE)

    #Save definitive tumour probability maps
    sitk.WriteImage(ADC_rBV_def, prob_imgs +'/ADC_rBV_def.nii') #Save ADC_rBV definitive probability map    
    sitk.WriteImage(ADC_rBF_def, prob_imgs +'/ADC_rBF_def.nii') #Save ADC_rBF definitive probability map
    
    #Read def tumour probability maps
    ADC_rBV_def = sitk.ReadImage(prob_imgs +'/ADC_rBV_def.nii')
    ADC_rBF_def = sitk.ReadImage(prob_imgs +'/ADC_rBF_def.nii')
    
    #Save segmentation of tumour probability map above certain threshold    
    ADC_rBV_90 = ADC_rBV_def >0.85
    ADC_rBF_90 = ADC_rBF_def >0.85
    
    sitk.WriteImage(ADC_rBV_90, prob_imgs +'/ADC_rBV_85.nii') #Save ADC_rBV probability map > 85% 
    sitk.WriteImage(ADC_rBF_90, prob_imgs +'/ADC_rBF_85.nii') #Save ADC_rBF probability map > 85%
    
#%% Generate tumour probability map only within FLAIR-enhancing tumour

    # #Generate Tumour_mask.nii on a FLAIR_inT1CE
    # # mkTUM = input('Generate Tumour_mask.nii on the FLAIR_bet image and save it in the Relative images folder. Once done type ENTER.')
    # # print(mkTUM)
    # tumour_mask = sitk.ReadImage(rel_imgs +'/Tumour_mask.nii') #read Tumour mask 
    # # tumour_mask = sitk.BinaryMedian(tumour_mask, (1,1,1), 1.0, 0.0)
    # # sitk.WriteImage(tumour_mask, rel_imgs +'/tumour_mask.nii')
    
    # #Apply tumour mask onto tumour probability maps
    # ADC_rBV_final = generate_mask(ADC_rBV_def, tumour_mask)
    # ADC_rBF_final = generate_mask(ADC_rBF_def, tumour_mask)
    
    # #Save final tumour probability maps in FLAIR-enhancing tumour
    # sitk.WriteImage(ADC_rBV_final, prob_imgs +'/ADC_rBV_probinFLAIR.nii') #Save ADC_rBV final probability map    
    # sitk.WriteImage(ADC_rBF_final, prob_imgs +'/ADC_rBF_probinFLAIR.nii') #Save ADC_rBF final probability map
    
    # #Save segmentation of tumour probability map above certain threshold    
    # ADC_rBV_95 = ADC_rBV_final >0.95
    # ADC_rBF_95 = ADC_rBF_final >0.95
    
    # sitk.WriteImage(ADC_rBV_95, prob_imgs +'/ADC_rBV_95_inFLAIR.nii') #Save ADC_rBV probability map > 95% 
    # sitk.WriteImage(ADC_rBF_95, prob_imgs +'/ADC_rBF_95_inFLAIR.nii') #Save ADC_rBF probability map > 95%
    
    # ADC_rBV_85 = ADC_rBV_final >0.85
    # ADC_rBF_85 = ADC_rBF_final >0.85
    
    # sitk.WriteImage(ADC_rBV_85, prob_imgs +'/ADC_rBV_85_inFLAIR.nii') #Save ADC_rBV probability map > 85% 
    # sitk.WriteImage(ADC_rBF_85, prob_imgs +'/ADC_rBF_85_inFLAIR.nii') #Save ADC_rBF probability map > 85%
    
    

#%% Generate tumour probability map only within GTV


    GTV = sitk.ReadImage(rel_imgs +'/GTV.nii') #read GTV 
    # PTV = sitk.BinaryDilate(GTV, (20,20,20))
    # sitk.WriteImage(PTV, rel_imgs +'/PTV.nii')  
    PTV = sitk.ReadImage(rel_imgs +'/PTV.nii')
    # Margin_volume=PTV-GTV
    # Margin_volume=generate_mask(Margin_volume, brain_mask)
    # sitk.WriteImage(Margin_volume, rel_imgs+'/Margin_volume.nii')
    Margin_volume = sitk.ReadImage(rel_imgs+'/Margin_volume.nii')
    
    # #Apply PTV mask onto tumour probability maps
    # ADC_rBV_final2 = generate_mask(ADC_rBV_def, PTV)
    # ADC_rBF_final2 = generate_mask(ADC_rBF_def, PTV)
    
    # #Save final tumour probability maps in PTV
    # sitk.WriteImage(ADC_rBV_final2, prob_imgs +'/ADC_rBV_probinPTV.nii') #Save ADC_rBV final2 probability map    
    # sitk.WriteImage(ADC_rBF_final2, prob_imgs +'/ADC_rBF_probinPTV.nii') #Save ADC_rBF final2 probability map
    
    # #Save segmentation of tumour probability map above certain threshold    
    # ADC_rBV_95_PTV = ADC_rBV_final2 >0.95
    # ADC_rBF_95_PTV = ADC_rBF_final2 >0.95
    
    # sitk.WriteImage(ADC_rBV_95_PTV, prob_imgs +'/ADC_rBV_95_inPTV.nii') #Save ADC_rBV probability map > 95% 
    # sitk.WriteImage(ADC_rBF_95_PTV, prob_imgs +'/ADC_rBF_95_inPTV.nii') #Save ADC_rBF probability map > 95%
    
    # ADC_rBV_85_PTV = ADC_rBV_final2 >0.85
    # ADC_rBF_85_PTV = ADC_rBF_final2 >0.85
    
    # sitk.WriteImage(ADC_rBV_85_PTV, prob_imgs +'/ADC_rBV_85_inPTV.nii') #Save ADC_rBV probability map > 85% 
    # sitk.WriteImage(ADC_rBF_85_PTV, prob_imgs +'/ADC_rBF_85_inPTV.nii') #Save ADC_rBF probability map > 85%


#%% Generate DPBN maps by using the following formula D_(p,i)= 〖D_min+(〖D_max-D〗_min)×〖TP〗_i〗^n 
    
    # os.mkdir(nii_dir +'Dose images')
    dose_imgs = nii_dir +'Dose images'
    
    # ADC_rBV_final2 = sitk.ReadImage(prob_imgs +'/ADC_rBV_probinPTV.nii')
    # ADC_rBF_final2 = sitk.ReadImage(prob_imgs +'/ADC_rBF_probinPTV.nii')
    
    # ADC_rBV_DPBN = 60 + 20*ADC_rBV_final2
    # ADC_rBF_DPBN = 60 + 20*ADC_rBF_final2
    
    # sitk.WriteImage(ADC_rBV_DPBN, dose_imgs +'/ADC_rBV_DPBN_inPTV.nii')
    # sitk.WriteImage(ADC_rBF_DPBN, dose_imgs +'/ADC_rBF_DPBN_inPTV.nii')
    
    #Read def tumour probability maps
    ADC_rBV_def = sitk.ReadImage(prob_imgs +'/ADC_rBV_def.nii')
    ADC_rBF_def = sitk.ReadImage(prob_imgs +'/ADC_rBF_def.nii')
    
    ADC_rBV_DPBN = 60 + 20*ADC_rBV_def
    ADC_rBF_DPBN = 60 + 20*ADC_rBF_def
    
    sitk.WriteImage(ADC_rBV_DPBN, dose_imgs +'/ADC_rBV_DPBN.nii')
    sitk.WriteImage(ADC_rBF_DPBN, dose_imgs +'/ADC_rBF_DPBN.nii')
    

    # ADC_rBV_DPBN = sitk.ReadImage(dose_imgs +'/ADC_rBV_DPBN.nii')
    # ADC_rBF_DPBN = sitk.ReadImage(dose_imgs +'/ADC_rBF_DPBN.nii')
    
    #Generate inverse dose prescription as Dinv = Dlow + Dbase - Dpresc, where Dlow = std dose to CTV, i.e. 60 Gy, and Dbase = 50 Gy
    Inv_ADC_rBV_DPBN = 60+50-ADC_rBV_DPBN
    Inv_ADC_rBF_DPBN = 60+50-ADC_rBF_DPBN
    
    sitk.WriteImage(Inv_ADC_rBV_DPBN, dose_imgs +'/Inv_ADC_rBV_DPBN.nii')
    sitk.WriteImage(Inv_ADC_rBF_DPBN, dose_imgs +'/Inv_ADC_rBF_DPBN.nii')
    
    #Generate mask of inverse dose in GTV
    Inv_ADC_rBV_DPBN_inGTV = generate_mask(Inv_ADC_rBV_DPBN, PTV)
    Inv_ADC_rBF_DPBN_inGTV = generate_mask(Inv_ADC_rBF_DPBN, PTV)
    
    sitk.WriteImage(Inv_ADC_rBV_DPBN_inGTV, dose_imgs +'/Inv_ADC_rBV_DPBN_inPTV.nii')
    sitk.WriteImage(Inv_ADC_rBF_DPBN_inGTV, dose_imgs +'/Inv_ADC_rBF_DPBN_inPTV.nii')






