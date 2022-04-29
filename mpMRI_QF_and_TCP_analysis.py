# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 16:01:24 2022

@author: cbri3325
"""


#%% Import functions 

import SimpleITK as sitk
import os
import pandas as pd
from ImageAnalysisFunctions import *
from ImageStatisticsFunctions import *

#%% Set Working directory
        
data_supradir = 'path to supradirectory containing individual patients directories' #Set working directory

subjs_path = [ f.path for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of the paths to the subjects directories
subjs_name = [ f.name for f in os.scandir(data_supradir) if f.is_dir() ] #Create a list of subjects names

Results = {'Quality Factor': pd.DataFrame(columns=['Subject_ID','QF STD GTV', 'QF DP GTV','QF STD CTV', 'QF DP CTV','QF STD PTV', 'QF DP PTV']),'Tumour Control Probability': pd.DataFrame(columns=['Subject_ID','TCP STD GTV', 'TCP DP GTV','TCP STD CTV', 'TCP DP CTV','TCP STD PTV', 'TCP DP PTV']) }

out_dir = 'C:/Users/cbri3325/OneDrive - The University of Sydney (Staff)/Caterina Brighi/Data/QIN/RT plans analysis/'
ResultsWriter = pd.ExcelWriter(out_dir +'QF_TCP_results.xlsx', engine='xlsxwriter')

#%%Create a for loop to perform image analysis on each subject sequentially

for current in subjs_name:

    subj_dir = data_supradir+current
    masked_dir = subj_dir+'/Dose masked to targets'
    
    DP_dir = masked_dir+'/Dose painting masked plans/'
    STD_dir = masked_dir+'/Standard masked plans/'
    TP_dir = masked_dir+'/Tumour probability masked/'
    struct_dir = subj_dir+'/RT structuresTP1/'
  
    print('Reading images for '+current) 
          
    DP_dose = sitk.ReadImage(DP_dir+'DP_Dose_PTV.nii') #Read DP dose in PTV
    DP_presc = sitk.ReadImage(DP_dir+'DP_Presc_PTV.nii') #Read DP prescription in PTV
    DP_presc.SetOrigin(DP_dose.GetOrigin())
    
    STD_dose = sitk.ReadImage(STD_dir+'STD_Dose_PTV.nii') #Read STD dose in PTV
    STD_presc = sitk.ReadImage(STD_dir+'STD_Presc_PTV.nii') #Read STD prescription in PTV
    STD_presc.SetOrigin(STD_dose.GetOrigin())
    
    TP = sitk.ReadImage(TP_dir+'TP_PTV.nii')
    
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
    
    #%%Generate Quality Factor Images and calculating QF
       
    print('Generating QF and QF images for '+current)   
    
    #Calculate QF for DP Plans
    QF_DP_GTV = calculate_QF(sitk.Cast(DP_presc, sitk.sitkFloat64), DP_dose, GTV)
    QF_DP_CTV = calculate_QF(sitk.Cast(DP_presc, sitk.sitkFloat64), DP_dose, CTV)
    QF_DP_PTV = calculate_QF(sitk.Cast(DP_presc, sitk.sitkFloat64), DP_dose, PTV)
    
    #Calculate QFfor STD Plans
    QF_STD_GTV = calculate_QF(sitk.Cast(STD_presc, sitk.sitkFloat64), STD_dose, GTV)
    QF_STD_CTV = calculate_QF(sitk.Cast(STD_presc, sitk.sitkFloat64), STD_dose, CTV)
    QF_STD_PTV = calculate_QF(sitk.Cast(STD_presc, sitk.sitkFloat64), STD_dose, PTV)
    
    #Generate QF images for DP Plans
    QF_DP_GTV_img = generate_QF_map(sitk.Cast(DP_presc, sitk.sitkFloat64), DP_dose, GTV)
    QF_DP_CTV_img = generate_QF_map(sitk.Cast(DP_presc, sitk.sitkFloat64), DP_dose, CTV)
    QF_DP_PTV_img = generate_QF_map(sitk.Cast(DP_presc, sitk.sitkFloat64), DP_dose, PTV)
    
    #Generate QF images for STD Plans
    QF_STD_GTV_img = generate_QF_map(sitk.Cast(STD_presc, sitk.sitkFloat64), STD_dose, GTV)
    QF_STD_CTV_img = generate_QF_map(sitk.Cast(STD_presc, sitk.sitkFloat64), STD_dose, CTV)
    QF_STD_PTV_img = generate_QF_map(sitk.Cast(STD_presc, sitk.sitkFloat64), STD_dose, PTV)
    
    
    #Save images
    
    sitk.WriteImage(QF_DP_GTV_img, DP_dir+'QF_DP_GTV.nii')
    sitk.WriteImage(QF_DP_CTV_img, DP_dir+'QF_DP_CTV.nii')
    sitk.WriteImage(QF_DP_PTV_img, DP_dir+'QF_DP_PTV.nii')
    
    sitk.WriteImage(QF_STD_GTV_img, STD_dir+'QF_STD_GTV.nii')
    sitk.WriteImage(QF_STD_CTV_img, STD_dir+'QF_STD_CTV.nii')
    sitk.WriteImage(QF_STD_PTV_img, STD_dir+'QF_STD_PTV.nii')

    
    #%%Calculate TCP within target volumes
    
    print('Calculating TCP within target volumes for '+current)
    
    #Calculate TCP in GTV
    DP_plan_GTV = allVoxInt(DP_dose, GTV)
    STD_plan_GTV = allVoxInt(STD_dose, GTV)
    TP_GTV = allVoxInt(TP, GTV)

    DP_TCP_GTV = calculate_tcp(DP_plan_GTV, 0.12, 0.03, TP_GTV)
    STD_TCP_GTV = calculate_tcp(STD_plan_GTV, 0.12, 0.03, TP_GTV)
    
    #Calculate TCP in PTV
    DP_plan_CTV = allVoxInt(DP_dose, CTV)
    STD_plan_CTV = allVoxInt(STD_dose, CTV)
    TP_CTV = allVoxInt(TP, CTV)

    DP_TCP_CTV = calculate_tcp(DP_plan_CTV, 0.12, 0.03, TP_CTV)
    STD_TCP_CTV = calculate_tcp(STD_plan_CTV, 0.12, 0.03, TP_CTV)
    
    #Calculate TCP in PTV
    DP_plan_PTV = allVoxInt(DP_dose, PTV)
    STD_plan_PTV = allVoxInt(STD_dose, PTV)
    TP_PTV = allVoxInt(TP, PTV)

    DP_TCP_PTV = calculate_tcp(DP_plan_PTV, 0.12, 0.03, TP_PTV)
    STD_TCP_PTV = calculate_tcp(STD_plan_PTV, 0.12, 0.03, TP_PTV)
    
    #%%Append QF and TCP values to Results dataframe
    
    Results['Quality Factor'] = Results['Quality Factor'].append({'Subject_ID': current,'QF STD GTV':QF_STD_GTV, 'QF DP GTV':QF_DP_GTV,'QF STD CTV':QF_STD_CTV, 'QF DP CTV':QF_DP_CTV,'QF STD PTV':QF_STD_PTV, 'QF DP PTV':QF_DP_PTV}, ignore_index=True)    
    Results['Tumour Control Probability'] = Results['Tumour Control Probability'].append({'Subject_ID': current,'TCP STD GTV':STD_TCP_GTV, 'TCP DP GTV':DP_TCP_GTV,'TCP STD CTV':STD_TCP_CTV, 'TCP DP CTV':DP_TCP_CTV,'TCP STD PTV':STD_TCP_PTV, 'TCP DP PTV':DP_TCP_PTV}, ignore_index=True)
    

    
    
#%%Save all dataframes to excel files here
print('Save all results to excel files')

for name, df in Results.items():
    df.to_excel(ResultsWriter, sheet_name=name, index=False)
ResultsWriter.save()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    