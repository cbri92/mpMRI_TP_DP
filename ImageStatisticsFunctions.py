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
from statistics import variance
from statistics import mean
import math
from scipy.stats import gaussian_kde
import matplotlib.cm as cm
import scipy.stats




def allVoxInt(image, roi):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.'''
    
    image = sitk.GetArrayFromImage(image)
    roi = sitk.GetArrayFromImage(roi)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    masked = image[~np.array(mask)]
    return masked

def allVoxInt_threshold(image, roi, threshold):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi with intensity value above the set threshold value applied to an image.'''
    
    image = sitk.GetArrayFromImage(image)
    roi = sitk.GetArrayFromImage(roi)
    mask = np.array(roi,dtype='bool')
    mask = np.invert(mask)
    masked = image[~np.array(mask)]
    list = masked[masked > threshold]
    return list

def ICC_ROI(image1_vxls, image2_vxls):
    
    '''This fucntion is used to calculate teh ICC between two images acquired at two timepoints in a test-retest experiment.
    The function returns a dictionary containing the follwoing parameters:
        'Between voxels variance'
        'Within voxels variance'
        'ROI ICC' '''

    
    image_mean_vxls = np.mean(np.array([image2_vxls,image1_vxls]), axis=0)
    image_betVox_Var = np.var(image_mean_vxls)
    image_diff_vxls = image1_vxls-image2_vxls
    image_withVox_Var = np.var(image_diff_vxls)
    wCV=(np.std(image_diff_vxls)/math.sqrt(2))/np.mean(image_mean_vxls, axis=0)
    RC =RC = 2.77*image_withVox_Var
    roi_ICC = image_betVox_Var/(image_betVox_Var+image_withVox_Var)
    
    return {'Between voxels variance':image_betVox_Var, 'Within voxels variance': image_withVox_Var, 'ROI ICC': roi_ICC, 'Within voxel CoV':wCV, 'Repeatability coefficient':RC}
        
def Bland_Altman_plt(image1_vxls, image2_vxls, color, size, title, outFilePath):
    
    '''This function plots the Blant-Altman plot for two repeated measurements of all the voxels in a ROI, with a scatter plot of colour specified'''
    
    # image1_vxls = allVoxInt(image1, roi)
    # image2_vxls = allVoxInt(image2, roi)
    image_mean_vxls = np.mean(np.array([image2_vxls,image1_vxls]), axis=0)
    image_diff_vxls = (image1_vxls-image2_vxls)
    
    mean = np.mean(image_diff_vxls)
    SD = np.std(image_diff_vxls)
       
    plt.scatter(image_mean_vxls, image_diff_vxls, s=size, marker='o', c=color)
    plt.title (title, size =16)
    plt.axhline(y = 0.0, color='black', linestyle = '-')
    plt.axhline(y = mean, color='black', linestyle = '--')
    plt.axhline(y = mean+1.96*SD, color='red', linestyle = '--')
    plt.axhline(y = mean-1.96*SD, color='red', linestyle = '--')
    plt.ylabel('Difference', size=16)
    plt.xlabel('Mean', size=16)
    # plt.xlim(0,image1_vxls.max())
    # plt.ylim(-image1_vxls.max(),image1_vxls.max())
    plt.tight_layout()
    plt.savefig(outFilePath)     
    plt.close()
    
def Bland_Altman_density_plt(image1_vxls, image2_vxls, size, title, outFilePath):
    
    '''This function plots the Blant-Altman density plot for two repeated measurements of all the voxels in a ROI, with a scatter plot of colour reflecting the density'''
    
    # image1_vxls = allVoxInt(image1, roi)
    # image2_vxls = allVoxInt(image2, roi)
    image_mean_vxls = np.mean(np.array([image2_vxls,image1_vxls]), axis=0)
    image_diff_vxls = (image1_vxls-image2_vxls)
    
    # Calculate the point density
    xy = np.vstack([image_mean_vxls,image_diff_vxls])
    z = 1 - gaussian_kde(xy)(xy) 
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    image_mean_vxls, image_diff_vxls, z = image_mean_vxls[idx], image_diff_vxls[idx], z[idx]

    fig, ax = plt.subplots()
    ax.scatter(image_mean_vxls, image_diff_vxls, c=z, s=size)
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
    plt.title (title, size=16)
    plt.plot(image1_vxls, image1_vxls, color='black')
    plt.ylabel('Timepoint 2', size=16)
    plt.xlabel('Timepoint 1', size=16)
    plt.tight_layout()
    plt.savefig(outFilePath)     
    plt.close()

def Scatter_plt(x, y, color, size, title, xlabel, ylabel, outFilePath):
    
    '''This function plots the scatter plot of colour specified'''
    
    # image1_vxls = allVoxInt(image1, roi)
    # image2_vxls = allVoxInt(image2, roi)
    
    plt.scatter(x, y, s=size, marker='o', c=color)
    plt.title (title, size=16)
    # plt.plot(x, y, color='black')
    plt.ylabel(ylabel, size=16)
    plt.xlabel(xlabel, size=16)
    plt.tight_layout()
    plt.savefig(outFilePath)     
    plt.close()
    
def Density_scatter_plt(image1_vxls, image2_vxls, size, title, xlabel, ylabel, outFilePath):
    
    '''This function generates a density scatter plot, with colors of points representing the density of each value'''
    
    # Calculate the point density
    xy = np.vstack([image1_vxls,image2_vxls])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    image1_vxls, image2_vxls, z = image1_vxls[idx], image2_vxls[idx], z[idx]

    fig, ax = plt.subplots()
    ax_ = ax.scatter(image1_vxls, image2_vxls, c=z, cmap=cm.plasma, s=size)
    plt.colorbar(ax_)
    plt.title (title, size=16)
    # plt.plot(image1_vxls, image1_vxls, color='black')
    # plt.axhline(y = 0.0, color='black', linestyle = '-')
    
    plt.xlim(60,80)
    plt.ylim(60,80)
    plt.ylabel(ylabel, size=16)
    plt.xlabel(xlabel, size=16)
    plt.savefig(outFilePath)     
    plt.close()


def Histogram_plot(x, density, bins, xlabel, ylabel, label, legend_position, fontsize, color, xlim_min, xlim_max, outFilePath):
    '''This function returns a histogram plot of x. If density is set to True the y axis will have the probability density. If set to False, it will have the number of counts.'''
    
    plt.hist(x, density=density, bins=bins, label=label, color=color)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.tight_layout()
    # mn, mx = plt.xlim()
    # plt.xlim(mn, mx)
    plt.xlim(xlim_min, xlim_max)
    plt.legend(loc=legend_position, fontsize=fontsize)
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
    df=n*(K-1)
    #find Chi-Square critical value
    chi2_975=scipy.stats.chi2.ppf(0.975, df=df)
    chi2_025=scipy.stats.chi2.ppf(0.025, df=df)
    RC_Lower = 2.77*math.sqrt((n*(K-1)*WMS)/chi2_975)
    RC_Upper = 2.77*math.sqrt((n*(K-1)*WMS)/chi2_025)
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
    
    bsVar = bSD**2
    wsVar = wSD**2
    
    df=n*(K-1)
    #find Chi-Square critical value
    chi2_975=scipy.stats.chi2.ppf(0.975, df=df)
    chi2_025=scipy.stats.chi2.ppf(0.025, df=df)
    RC = 2.77*wSD
    RC_Lower = 2.77*math.sqrt((n*(K-1)*WMS)/chi2_975)
    RC_Upper = 2.77*math.sqrt((n*(K-1)*WMS)/chi2_025)
    wCV=wSD/overall_mean
    
    ICC = (tSD**2-wSD**2)/tSD**2
    
    return {'BMS':BMS, 'WMS':WMS, 'tSD':tSD, 'bSD':bSD, 'wSD':wSD, 'Within subject variance':wsVar, 'Between subject variance':bsVar, 'RC':RC, 'RC Lower':RC_Lower, 'RC Upper':RC_Upper, 'wCV':wCV, 'ICC':ICC}

#%%Statistical tests

def f_test(x, y):
    x = np.array(x)
    y = np.array(y)
    f = np.var(x, ddof=1)/np.var(y, ddof=1) #calculate F test statistic 
    dfn = x.size-1 #define degrees of freedom numerator 
    dfd = y.size-1 #define degrees of freedom denominator 
    p = 1-scipy.stats.f.cdf(f, dfn, dfd) #find p-value of F test statistic 
    return f, p

#%%RT plans evaluation functions

def calculate_QF(Prescribed_dose_image, Planned_dose_image, roi):
    
    '''This function calculates the Quality Factor between the planned and the prescribed dose within a roi. The closer QF is to 100, the better the quality of the plan wrt to the dose prescription.'''
    
    Rel_dose = abs((Planned_dose_image-Prescribed_dose_image)/Prescribed_dose_image)
    vector = allVoxInt(Rel_dose, roi)
    QF = 100 - mean(vector)*100
    return QF