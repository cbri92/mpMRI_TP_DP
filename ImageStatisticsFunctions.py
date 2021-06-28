# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:23:04 2021

@author: cbri3325
"""


#%% Set working environments 

import matplotlib.pyplot as plt
import matplotlib as mpl
import SimpleITK as sitk
import numpy as np
import scipy.ndimage as ndimage
import pandas as pd
import datetime
import os
import glob
import gzip
import shutil
import xlsxwriter
from scipy.stats.stats import pearsonr
from multiprocessing.pool import ThreadPool
from functools import partial
from skimage.feature import peak_local_max
from statistics import pvariance
from statistics import mean
import math


def allVoxInt(image, roi):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.
    The process requires the generation of a mask (hence a maskname) and its writing into and reading from a temp folder (provide fullpath to temp folder).'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image,roi) #Generate a mask of the roi applied to image
    array = sitk.GetArrayFromImage(mask) #convert image to array
    flat = array.flatten() #flatten array to 1D
    list = flat[flat != 0]
    return list

def ICC_ROI(image1_vxls, image2_vxls):
    
    '''This fucntion is used to calculate teh ICC between two images acquired at two timepoints in a test-retest experiment.
    The function returns a dictionary containing the follwoing parameters:
        'Between voxels variance'
        'Within voxels variance'
        'ROI ICC' '''

    
    image_mean_vxls = (image2_vxls+image1_vxls)/2
    image_betVox_Var = pvariance(image_mean_vxls)
    image_diff_vxls = image1_vxls-image2_vxls
    image_withVox_Var = pvariance(image_diff_vxls)
    roi_ICC = image_betVox_Var/(image_betVox_Var+image_withVox_Var)
    
    return {'Between voxels variance':image_betVox_Var, 'Within voxels variance': image_withVox_Var, 'ROI ICC': roi_ICC}
        
def Bland_Altman_plt(image1_vxls, image2_vxls, color, size, title, outFilePath):
    
    '''This function plots the Blant-Altman plot for two repeated measurements of all the voxels in a ROI, with a scatter plot of colour specified'''
    
    # image1_vxls = allVoxInt(image1, roi)
    # image2_vxls = allVoxInt(image2, roi)
    image_mean_vxls = (image2_vxls+image1_vxls)/2
    image_diff_vxls = (image1_vxls-image2_vxls)
    
    plt.scatter(image_mean_vxls, image_diff_vxls, s=size, marker='o', c=color)
    plt.title (title)
    plt.axhline(y = 0.0, color='black', linestyle = '-')
    plt.ylabel('Difference (Timepoint 1 - Timepoint 2)')
    plt.xlabel('Mean (Timepoint 1 and Timepoint 2)')
    plt.savefig(outFilePath)     
    plt.close()
    
def Vox_Corrleation_plt(image1_vxls, image2_vxls, color, size, title, outFilePath):
    
    '''This function plots the voxelwise correlation plot for two repeated measurements of all the voxels in a ROI, with a scatter plot of colour specified'''
    
    # image1_vxls = allVoxInt(image1, roi)
    # image2_vxls = allVoxInt(image2, roi)
    
    plt.scatter(image1_vxls, image2_vxls, s=size, marker='o', c=color)
    plt.title (title)
    plt.plot(image1_vxls, image1_vxls, color='black')
    plt.ylabel('Timepoint 2')
    plt.xlabel('Timepoint 1')
    plt.savefig(outFilePath)     
    plt.close()
    
#%%Repeatibility metrics on a group level

def BMS(K, repeat1, repeat2):
    '''This function calculates the between subjects mean given K number of repeats, and two columns array representing repeat 1 and repeat 2.'''
    
    mean_of_repeats = (repeat1+repeat2)/K
    n=len(repeat1)
    overall_mean = (sum(repeat1)+sum(repeat2))/(n*K)
    x=0
    for i in mean_of_repeats:
        x=x+(((i-overall_mean)**2)/n)
    BMS=K*x
    return BMS
        
def WMS(K, repeat1, repeat2):
    '''This function calculates the within subject mean given K number of repeats, and two columns array representing repeat 1 and repeat 2.'''
    
    mean_of_repeats = (repeat1+repeat2)/K
    n=len(repeat1)
    x=0
    for i,r1,r2 in zip(mean_of_repeats, repeat1, repeat2):
        x=x+(((r1-i)**2)/(n*(K-1))+((r2-i)**2)/(n*(K-1)))
    WMS=x
    return WMS
        
def tSD(K, BMS, WMS):
    '''This function returns the total standard deviation given the number of repeats K, the between subject mean BMS and the within subject mean WMS.'''
    
    tSD = math.sqrt(((BMS+(K-1)*WMS)/K))
    return tSD

def bSD(K, BMS, WMS):
    '''This function returns the between subjects standard deviation given the number of repeats K, the between subject mean BMS and the within subject mean WMS.'''
    
    bSD = math.sqrt(((BMS-WMS)/K))
    return bSD

def wSD(WMS):
    '''This function returns the within subject standard deviation given the within subject mean WMS.'''
    
    wSD = math.sqrt(WMS)
    return wSD

def RC(wSD, WMS, n, K):
    '''This function returns the repeatibility coefficient and the upper and lower confidence intervals given the within subject standard deviation wSD, the within subject mean WMS, the number of subjects n and the number of repeats K.'''
    RC = 2.77*wSD
    RC_Lower = 2.77*math.sqrt((n*(K-1)*WMS)/3.82)
    RC_Upper = 2.77*math.sqrt((n*(K-1)*WMS)/21.92)
    return RC, RC_Lower, RC_Upper

def wCV(wSD, Overall_mean):
    '''This function returns the within subject coefficient of variation given the within subject standard deviation and the overall mean.'''
    wCV=wSD/Overall_mean
    return wCV

def Repeatibility_metrics(K, repeat1, repeat2):
    '''This function returns a dictionary with the following repeatibility metrics:
        'BMS: Between subjects mean'
        'WMS: Within subject mean'
        'tSD: Total standard deviation'
        'bSD: Between subjects standard deviation'
        'wSD: Within subject standard deviation'
        'RC: Repeatibility coefficient'
        'RC Lower: Lower 95% CI'
        'RC Upper: Upper 95% CI'
        'wCV: Within subject coefficient of variation'
        'ICC: intra-correlation coefficient'
        '''
    mean_of_repeats = (repeat1+repeat2)/K
    n=len(repeat1)
    overall_mean = (sum(repeat1)+sum(repeat2))/(n*K)
    x=0
    for i in mean_of_repeats:
        x=x+(((i-overall_mean)**2)/n)
    BMS=K*x
    
    y=0
    for i,r1,r2 in zip(mean_of_repeats, repeat1, repeat2):
        y=y+(((r1-i)**2)/(n*(K-1))+((r2-i)**2)/(n*(K-1)))
    WMS=y
    
    tSD = math.sqrt(((BMS+(K-1)*WMS)/K))
    bSD = math.sqrt(((BMS-WMS)/K))
    wSD = math.sqrt(WMS)
    RC = 2.77*wSD
    RC_Upper = 2.77*math.sqrt((n*(K-1)*WMS)/3.82)
    RC_Lower = 2.77*math.sqrt((n*(K-1)*WMS)/21.92)
    wCV=wSD/overall_mean
    
    BetwSubj_Var = pvariance(mean_of_repeats)
    repeats_diff = repeat1-repeat2
    WithSubj_Var = pvariance(repeats_diff)
    ICC = BetwSubj_Var/(BetwSubj_Var+WithSubj_Var)
    
    return {'BMS':BMS, 'WMS':WMS, 'tSD':tSD, 'bSD':bSD, 'wSD':wSD, 'RC':RC, 'RC Lower':RC_Lower, 'RC Upper':RC_Upper, 'wCV':wCV, 'ICC':ICC}
