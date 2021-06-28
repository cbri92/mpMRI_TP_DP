# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:54:41 2020

@author: s4451992
"""

import numpy as np
import pandas as pd
import SimpleITK as sitk
import os
from ImageResultsVisualizationFunctions import *

data_supradir = 'C:/Users/cbri3325/Dropbox (Sydney Uni)/Caterina Brighi/Data/test/' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a lisdt of subjects names

#Make sure that the data_supradir contains an excel file with the number of z_slice with tumour in CT_FET space for each patient

#Make a directory for the images
os.mkdir(data_supradir +'/Visualization Results Images')
imgs_dir = data_supradir +'/Visualization Results Images' #Set paths to plot dir

#Import z_slices values into a dataframe
z_slices = pd.read_excel(data_supradir + 'z_slices.xlsx', sheet_name='Sheet1', index_col='Subject_ID')

for current in subjs_name:
    
    subj_dir = data_supradir+current
    subj_name = current
    
    print('Starting generating results figures for '+current)
    
    nii_dir = subj_dir +'/nii/' 
   
    #Set paths to subfolders    
    bet_imgs = nii_dir +'BET images'
    prob_imgs = nii_dir +'Probability images'
    rel_imgs = nii_dir +'Relative images'
    dose_imgs = nii_dir +'Dose images'
        
    #Read resampled bet images   
    T1CE_bet = sitk.ReadImage(bet_imgs +'/T1_bet.nii') #read T1CE_bet
    FLAIR_bet = sitk.ReadImage(bet_imgs +'/FLAIR_bet.nii') #read T1CE_bet
    
    #Read def tumour probability maps
    ADC_rBV_def = sitk.ReadImage(prob_imgs +'/ADC_rBV_def.nii')
    ADC_rBF_def = sitk.ReadImage(prob_imgs +'/ADC_rBF_def.nii')
    
    #Read final tumour probability maps (masked in FLAIR)
    ADC_rBV_final = sitk.ReadImage(prob_imgs +'/ADC_rBV_probinFLAIR.nii')
    ADC_rBF_final = sitk.ReadImage(prob_imgs +'/ADC_rBF_probinFLAIR.nii')
    
    #Read segmentations of >85% and >95% tumour probability in FLAIR
    ADC_rBV_85 = sitk.ReadImage(prob_imgs +'/ADC_rBV_85_inFLAIR.nii')
    ADC_rBF_85 = sitk.ReadImage(prob_imgs +'/ADC_rBF_85_inFLAIR.nii')
    
    ADC_rBV_95 = sitk.ReadImage(prob_imgs +'/ADC_rBV_95_inFLAIR.nii')
    ADC_rBF_95 = sitk.ReadImage(prob_imgs +'/ADC_rBF_95_inFLAIR.nii')
    
    #Read final tumour probability maps (masked in PTV)
    ADC_rBV_final2 = sitk.ReadImage(prob_imgs +'/ADC_rBV_probinPTV.nii')
    ADC_rBF_final2 = sitk.ReadImage(prob_imgs +'/ADC_rBF_probinPTV.nii')
    
    #Read segmentations of >85% and >95% tumour probability
    ADC_rBV_85_PTV = sitk.ReadImage(prob_imgs +'/ADC_rBV_85_inPTV.nii')
    ADC_rBF_85_PTV = sitk.ReadImage(prob_imgs +'/ADC_rBF_85_inPTV.nii')
    
    ADC_rBV_95_PTV = sitk.ReadImage(prob_imgs +'/ADC_rBV_95_inPTV.nii')
    ADC_rBF_95_PTV = sitk.ReadImage(prob_imgs +'/ADC_rBF_95_inPTV.nii')
    
    #Read DBPN maps (masked in PTV)
    ADC_rBV_DPBN = sitk.ReadImage(dose_imgs +'/ADC_rBV_DPBN_inPTV.nii')
    ADC_rBF_DPBN = sitk.ReadImage(dose_imgs +'/ADC_rBF_DPBN_inPTV.nii')
    
    #Read brain mask and tumour mask
    brain_mask = sitk.ReadImage(bet_imgs +'/brain_mask.nii')
    tumour_mask = sitk.ReadImage(rel_imgs+'/Tumour_mask.nii')
    PTV = sitk.ReadImage(rel_imgs+'/PTV.nii')
        
        
    #Make a folder of current in visualisation image folder
    os.mkdir(imgs_dir +'/'+ current)
    img_current = imgs_dir +'/'+ current
        
    #Generate images + respective overlays
    slice_n = int(z_slices.loc[current,'z_slice'])
    
    brain_mask_coronal = sitk.Cast(brain_mask, sitk.sitkUInt8)
    brain_mask_coronal = brain_mask_coronal[:,:,slice_n]
    # sitk.WriteImage(brain_mask_coronal, img_current +'/'+ subj_name +'_brain_mask.png')

    tumour_mask_coronal = sitk.Cast(tumour_mask, sitk.sitkUInt8)
    tumour_mask_coronal = tumour_mask_coronal[:,:,slice_n]
    # sitk.WriteImage(tumour_mask_coronal, img_current +'/'+ subj_name +'_tumour_mask.png')  
    
    PTV_coronal = sitk.Cast(PTV, sitk.sitkUInt8)
    PTV_coronal = PTV_coronal[:,:,slice_n]
    # sitk.WriteImage(PTV_coronal, img_current +'/'+ subj_name +'_PTV.png')
    
    
    # T1CE_255 = window_img(T1CE_bet)
    # T1CE_coronal = T1CE_255[:,:,slice_n]
    T1CE_coronal = grayScale_to_ColorMap(T1CE_bet, slice_n, 'Grey')
    T1CE_coronal = mask_image_multiply(brain_mask_coronal, T1CE_coronal)
    sitk.WriteImage(T1CE_coronal, img_current +'/'+ subj_name +'_T1CE.png')
    # change_brightness(img_current +'/'+ subj_name +'_T1CE.png', 1.5)
    
    # FLAIR_255 = window_img(FLAIR_bet)
    # FLAIR_coronal = FLAIR_255[:,:,slice_n]
    FLAIR_coronal = grayScale_to_ColorMap(FLAIR_bet, slice_n, 'Grey')
    FLAIR_coronal = mask_image_multiply(brain_mask_coronal, FLAIR_coronal)
    sitk.WriteImage(FLAIR_coronal, img_current +'/'+ subj_name +'_FLAIR.png')
    # change_brightness(img_current +'/'+ subj_name +'_FLAIR.png', 2.5)
    
    ADC_rBV_def_coronal = grayScale_to_ColorMap(ADC_rBV_def, slice_n, 'Jet')
    ADC_rBV_def_coronal = mask_image_multiply(brain_mask_coronal, ADC_rBV_def_coronal)
    sitk.WriteImage(ADC_rBV_def_coronal, img_current +'/'+ subj_name +'_ADC_rBV_def.png')
    
    ADC_rBF_def_coronal = grayScale_to_ColorMap(ADC_rBF_def, slice_n, 'Jet')
    ADC_rBF_def_coronal = mask_image_multiply(brain_mask_coronal, ADC_rBF_def_coronal)
    sitk.WriteImage(ADC_rBF_def_coronal, img_current +'/'+ subj_name +'_ADC_rBF_def.png')
    
    #Images in FLAIR enhancing mask
    
    ADC_rBV_final_coronal = grayScale_to_ColorMap(ADC_rBV_final, slice_n, 'Jet')
    ADC_rBV_final_coronal = mask_image_multiply(tumour_mask_coronal, ADC_rBV_final_coronal)
    sitk.WriteImage(ADC_rBV_final_coronal, img_current +'/'+ subj_name +'_ADC_rBV_probinFLAIR.png')
    
    ADC_rBF_final_coronal = grayScale_to_ColorMap(ADC_rBF_final, slice_n, 'Jet')
    ADC_rBF_final_coronal = mask_image_multiply(tumour_mask_coronal, ADC_rBF_final_coronal)
    sitk.WriteImage(ADC_rBF_final_coronal, img_current +'/'+ subj_name +'_ADC_rBF_probinFLAIR.png')
    
    ADC_rBV_85_coronal = grayScale_to_ColorMap(ADC_rBV_85, slice_n, 'Copper')
    ADC_rBV_85_coronal = mask_image_multiply(tumour_mask_coronal, ADC_rBV_85_coronal)
    sitk.WriteImage(ADC_rBV_85_coronal, img_current +'/'+ subj_name +'_ADC_rBV_85_inFLAIR.png')
    
    ADC_rBF_85_coronal = grayScale_to_ColorMap(ADC_rBF_85, slice_n, 'Copper')
    ADC_rBF_85_coronal = mask_image_multiply(tumour_mask_coronal, ADC_rBF_85_coronal)
    sitk.WriteImage(ADC_rBF_85_coronal, img_current +'/'+ subj_name +'_ADC_rBF_85_inFLAIR.png')
    
    ADC_rBV_95_coronal = grayScale_to_ColorMap(ADC_rBV_95, slice_n, 'Red')
    ADC_rBV_95_coronal = mask_image_multiply(tumour_mask_coronal, ADC_rBV_95_coronal)
    sitk.WriteImage(ADC_rBV_95_coronal, img_current +'/'+ subj_name +'_ADC_rBV_95_inFLAIR.png')
    
    ADC_rBF_95_coronal = grayScale_to_ColorMap(ADC_rBF_95, slice_n, 'Red')
    ADC_rBF_95_coronal = mask_image_multiply(tumour_mask_coronal, ADC_rBF_95_coronal)
    sitk.WriteImage(ADC_rBF_95_coronal, img_current +'/'+ subj_name +'_ADC_rBF_95_inFLAIR.png')
    
    #Images in PTV mask
    
    ADC_rBV_final2_coronal = grayScale_to_ColorMap(ADC_rBV_final2, slice_n, 'Jet')
    ADC_rBV_final2_coronal = mask_image_multiply(PTV_coronal, ADC_rBV_final2_coronal)
    sitk.WriteImage(ADC_rBV_final2_coronal, img_current +'/'+ subj_name +'_ADC_rBV_probinPTV.png')
    
    ADC_rBF_final2_coronal = grayScale_to_ColorMap(ADC_rBF_final2, slice_n, 'Jet')
    ADC_rBF_final2_coronal = mask_image_multiply(PTV_coronal, ADC_rBF_final2_coronal)
    sitk.WriteImage(ADC_rBF_final2_coronal, img_current +'/'+ subj_name +'_ADC_rBF_probinPTV.png')
    
    ADC_rBV_85_PTV_coronal = grayScale_to_ColorMap(ADC_rBV_85_PTV, slice_n, 'Copper')
    ADC_rBV_85_PTV_coronal = mask_image_multiply(PTV_coronal, ADC_rBV_85_PTV_coronal)
    sitk.WriteImage(ADC_rBV_85_PTV_coronal, img_current +'/'+ subj_name +'_ADC_rBV_85_inPTV.png')
    
    ADC_rBF_85_PTV_coronal = grayScale_to_ColorMap(ADC_rBF_85_PTV, slice_n, 'Copper')
    ADC_rBF_85_PTV_coronal = mask_image_multiply(PTV_coronal, ADC_rBF_85_PTV_coronal)
    sitk.WriteImage(ADC_rBF_85_PTV_coronal, img_current +'/'+ subj_name +'_ADC_rBF_85_inPTV.png')
    
    ADC_rBV_95_PTV_coronal = grayScale_to_ColorMap(ADC_rBV_95_PTV, slice_n, 'Red')
    ADC_rBV_95_PTV_coronal = mask_image_multiply(PTV_coronal, ADC_rBV_95_PTV_coronal)
    sitk.WriteImage(ADC_rBV_95_PTV_coronal, img_current +'/'+ subj_name +'_ADC_rBV_95_inPTV.png')
    
    ADC_rBF_95_PTV_coronal = grayScale_to_ColorMap(ADC_rBF_95_PTV, slice_n, 'Red')
    ADC_rBF_95_PTV_coronal = mask_image_multiply(PTV_coronal, ADC_rBF_95_PTV_coronal)
    sitk.WriteImage(ADC_rBF_95_PTV_coronal, img_current +'/'+ subj_name +'_ADC_rBF_95_inPTV.png')
    
    ADC_rBV_DPBN_coronal = grayScale_to_ColorMap(ADC_rBV_DPBN, slice_n, 'Jet')
    ADC_rBV_DPBN_coronal = mask_image_multiply(PTV_coronal, ADC_rBV_DPBN_coronal)
    sitk.WriteImage(ADC_rBV_DPBN_coronal, img_current +'/'+ subj_name +'_ADC_rBV_DPBNinPTV.png')
    
    ADC_rBF_DPBN_coronal = grayScale_to_ColorMap(ADC_rBF_DPBN, slice_n, 'Jet')
    ADC_rBF_DPBN_coronal = mask_image_multiply(PTV_coronal, ADC_rBF_DPBN_coronal)
    sitk.WriteImage(ADC_rBF_DPBN_coronal, img_current +'/'+ subj_name +'_ADC_rBF_DPBNinPTV.png')
    
    
    
    
    FLAIR_ADC_rBV_def_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_def.png', 0.4)
    FLAIR_ADC_rBF_def_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_def.png', 0.4)
    
    
    FLAIR_ADC_rBV_final_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_probinFLAIR.png', 0.5)
    FLAIR_ADC_rBF_final_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_probinFLAIR.png', 0.5)
    
    FLAIR_ADC_rBV_85_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_85_inFLAIR.png', 0.4)
    FLAIR_ADC_rBV_95_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_95_inFLAIR.png', 0.4)
    FLAIR_ADC_rBV_8595_ovl = simple_blend(img_current +'/'+ subj_name +'_ADC_rBV_85_inFLAIR_overlay.png', img_current +'/'+ subj_name +'_ADC_rBV_95_inFLAIR.png', 0.3, '85_overlay')
    
    FLAIR_ADC_rBF_85_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_85_inFLAIR.png', 0.4)
    FLAIR_ADC_rBF_95_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_95_inFLAIR.png', 0.4)
    FLAIR_ADC_rBF_8595_ovl = simple_blend(img_current +'/'+ subj_name +'_ADC_rBF_85_inFLAIR_overlay.png', img_current +'/'+ subj_name +'_ADC_rBF_95_inFLAIR.png', 0.3, '85_overlay')
    
    
    
    FLAIR_ADC_rBV_final2_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_probinPTV.png', 0.5)
    FLAIR_ADC_rBF_final2_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_probinPTV.png', 0.5)
    
    FLAIR_ADC_rBV_85_PTV_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_85_inPTV.png', 0.4)
    FLAIR_ADC_rBV_95_PTV_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_95_inPTV.png', 0.4)
    FLAIR_ADC_rBV_8595_PTV_ovl = simple_blend(img_current +'/'+ subj_name +'_ADC_rBV_85_inPTV_overlay.png', img_current +'/'+ subj_name +'_ADC_rBV_95_inPTV.png', 0.3, '85_overlay')
    
    FLAIR_ADC_rBF_85_PTV_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_85_inPTV.png', 0.4)
    FLAIR_ADC_rBF_95_PTV_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_95_inPTV.png', 0.4)
    FLAIR_ADC_rBF_8595_PTV_ovl = simple_blend(img_current +'/'+ subj_name +'_ADC_rBF_85_inPTV_overlay.png', img_current +'/'+ subj_name +'_ADC_rBF_95_inPTV.png', 0.3, '85_overlay')
    
    FLAIR_ADC_rBV_DPBN_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBV_DPBNinPTV.png', 0.5)
    FLAIR_ADC_rBF_DPBN_ovl = simple_blend(img_current +'/'+ subj_name +'_FLAIR.png', img_current +'/'+ subj_name +'_ADC_rBF_DPBNinPTV.png', 0.5)
    

    # FLAIR_ADC_rBV_def_ovl = colormap_imgs_overlay(FLAIR_bet, ADC_rBV_def, slice_n, 'Jet', alpha_value=0.50)    
    # sitk.WriteImage(FLAIR_ADC_rBV_def_ovl, img_current +'/'+ subj_name +'_FLAIR_ADC_rBV_def_ovl.png')
    
    # FLAIR_ADC_rBF_def_ovl = colormap_imgs_overlay(FLAIR_bet, ADC_rBF_def, slice_n, 'Jet', alpha_value=0.50)    
    # sitk.WriteImage(FLAIR_ADC_rBF_def_ovl, img_current +'/'+ subj_name +'_FLAIR_ADC_rBF_def_ovl.png')
    
    # FLAIR_ADC_rBV_final_ovl = colormap_imgs_overlay(FLAIR_bet, ADC_rBV_final, slice_n, 'Jet', alpha_value=0.50)    
    # sitk.WriteImage(FLAIR_ADC_rBV_final_ovl, img_current +'/'+ subj_name +'_FLAIR_ADC_rBV_final_ovl.png')
    
    # FLAIR_ADC_rBF_final_ovl = colormap_imgs_overlay(FLAIR_bet, ADC_rBF_final, slice_n, 'Jet', alpha_value=0.50)    
    # sitk.WriteImage(FLAIR_ADC_rBF_final_ovl, img_current +'/'+ subj_name +'_FLAIR_ADC_rBF_final_ovl.png')
    

