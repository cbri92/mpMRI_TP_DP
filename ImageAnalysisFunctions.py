# -*- coding: utf-8 -*-
"""
Created on Tue May  5 11:17:17 2020

@author: s4451992
"""

#%% Set working environments 

import matplotlib.pyplot as plt
import matplotlib as mpl
import dicom2nifti
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


#%% dicom to nifti conversion

#Converting a directory with dicom files to nifti files
def dcm2niiDir(dicom_directory, output_folder):
    
    '''Converts a directory with dicom files to nifti files'''
    
    return dicom2nifti.convert_directory(dicom_directory, output_folder)

#Converting a directory with only 1 series to 1 nifti file
def dcm2nii(dicom_directory, output_file):
    
    '''Converts a directory with only 1 series to 1 nifti file'''
    
    return dicom2nifti.dicom_series_to_nifti(dicom_directory, output_file, reorient_nifti=True)

#%% Gunzip all .gz files

def gunzip_shutil(source_filepath, dest_filepath, block_size=65536):
    
    '''Extract all .gz compressed files and delete original .gz files'''
    
    with gzip.open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        shutil.copyfileobj(s_file, d_file, block_size)
    os.remove(source_filepath) #remove original .gz files
    
#%% Rename images in subj_dir with sensible name

def rename_img_files(subj_dir, subj_name):
    
    '''Renames all .nii images to sensible names including subj_name and scan type'''
    
    # for filename in glob.glob(subj_dir +'/' +'*.nii'):
    #     if "mprage" in filename:
    #         os.rename(filename, subj_dir +'/'+ subj_name +'_T1CE.nii')
    #     elif "FET_PET" in filename:
    #         os.rename(filename, subj_dir +'/'+ subj_name +'_PET_FET_raw.nii')
    #     elif "FET_CT" in filename:
    #         os.rename(filename, subj_dir +'/'+ subj_name +'_CT_FET.nii')
    #     elif "PSMA_PET" in filename:
    #         os.rename(filename, subj_dir +'/'+ subj_name +'_PET_TER_raw.nii')
    #     elif "PSMA_CT" in filename:
    #         os.rename(filename, subj_dir +'/'+ subj_name +'_CT_TER.nii')
    
    for filename in glob.glob(subj_dir +'/' +'*.nii'):
        if "FET PET Brain" in filename:
            os.rename(filename, subj_dir +'/'+ subj_name +'_PET_FET_raw.nii')
        elif "FET Fused CT" in filename:
            os.rename(filename, subj_dir +'/'+ subj_name +'_CT_FET.nii')
        elif "PSMA PET" in filename:
            os.rename(filename, subj_dir +'/'+ subj_name +'_PET_TER_raw.nii')
        elif "PSMA Fused CT" in filename:
            os.rename(filename, subj_dir +'/'+ subj_name +'_CT_TER.nii')


#%% PET preprocessing
    
def PET_preprocessing(subj_dir, subj_name):
    
    '''Generates a dictionary containing PET acquisition parameters, exports these data into excel file,
       then uses these data to decay correct the raw PET images and writes them as PET_FET.nii nd PET_TER.nii images in the subj_dir folder'''
    
    subjBodyWeight = float(open(subj_dir +'/'+ 'subject_weight.txt').read()) #import patient's body weight in grams
    
    tracers = ['F', 'Ga']
    param_names = ['Injection dose [Bq]', 'Injection time [yyyy-mm-dd HH:MM:SS]', 'Scan time [yyyy-mm-dd HH:MM:SS]', 'Time passed [minutes]', 'Half-life [minutes]', 'Units conversion factor']
    
    half_life_minutes = {'F': 109.77, 'Ga': 67.71}
    ln2 = 0.69314718056
    
    dictionary = {}
    for element in tracers:
        for filename in glob.glob(subj_dir +'/' +'*.txt'):
            if element+".txt" in filename:
                if "injection_dose" in filename:
                    InjDose = float(open(filename).read())
                elif "injection_time" in filename:
                    InjTimeFile = open(filename, 'r').read()
                    InjTime = datetime.datetime.strptime(InjTimeFile, '%Y %m %d %H %M %S')
                elif "scan_time" in filename: 
                    ScanTimeFile = open(filename, 'r').read()
                    ScanTime = datetime.datetime.strptime(ScanTimeFile, '%Y %m %d %H %M %S')
        
        TimePassed = ScanTime - InjTime
        TimePassed_minutes = TimePassed.total_seconds()/60
        
        decayConstant = ln2/half_life_minutes[element]
        decayFactor = np.exp(decayConstant*TimePassed_minutes)
        unitsFactor = subjBodyWeight/InjDose
        ConversionFactor = decayFactor*unitsFactor
        
        
                
        dictionary[element] = [InjDose, InjTime, ScanTime, TimePassed_minutes,\
                  half_life_minutes[element], ConversionFactor]

#    print(dictionary)

    InfoFile = pd.DataFrame.from_dict(dictionary)       
    InfoFile.index = param_names
    InfoFile.to_excel(subj_dir +'/'+'InfoFile_Subject_'+subj_name+'.xlsx')

    #Multiply PET raw images by conversion factor
    FconvFact = InfoFile.loc["Units conversion factor", "F"]
    GaconvFact = InfoFile.loc["Units conversion factor", "Ga"]
    
    PET_FET_raw = sitk.ReadImage(subj_dir +'/'+ subj_name +'_PET_FET_raw.nii') #read PET FET raw image
    PET_FET_conv = sitk.ShiftScale(PET_FET_raw, shift = 0, scale = FconvFact) #multiply PET FET image by F conversion factor
    sitk.WriteImage(PET_FET_conv, subj_dir +'/'+ subj_name +'_PET_FET.nii') #write converted image as a nifti file PET FET
    
    PET_TER_raw = sitk.ReadImage(subj_dir +'/'+ subj_name +'_PET_TER_raw.nii') #read PET TER image
    PET_TER_conv = sitk.ShiftScale(PET_TER_raw, shift = 0, scale = GaconvFact) #multiply PET TER image by Ga conversion factor
    sitk.WriteImage(PET_TER_conv, subj_dir +'/'+ subj_name +'_PET_TER.nii') #write converted image as a nifti file PET TER
    
#%% Registration and resampling of images
    
def Resample_image(input_image, reference_image):
    
    '''Returns the input image resampled to the reference image space.
       Remember to write the output image into an image file after applying this function'''
       
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(reference_image)   
    output_image = resample.Execute(input_image)
    return output_image

def Resample_image_to_specified_resolution(input_image, x_spacing, y_spacing, z_spacing):
    
    '''Resample the input image to the specified resolution, maintaining the same origin, size and direction.'''
    
    resample = sitk.ResampleImageFilter()
    resample.SetOutputOrigin(input_image.GetOrigin())
    resample.SetOutputDirection(input_image.GetDirection())
    resample.SetOutputSpacing((x_spacing, y_spacing, z_spacing))
    resample.SetSize((int((input_image.GetSize()[0]*input_image.GetSpacing()[0])/x_spacing), int((input_image.GetSize()[1]*input_image.GetSpacing()[1])/y_spacing), int((input_image.GetSize()[2]*input_image.GetSpacing()[2])/z_spacing)))
    output_image = resample.Execute(input_image)
    return output_image



def start_observer():
            global metricvalue_parameters_list 
            metricvalue_parameters_list = []
    
def iteration_observer(registration_method):    
    metricvalue_parameters_list.append((registration_method.GetMetricValue(), registration_method.GetOptimizerPosition()))
    

def multires_registration(fixed_image, moving_image, initial_transform):
    
    '''This is the registration configuration which we use in all cases. The only parameter that we vary is the initial_transform.'''
    
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=100, estimateLearningRate=registration_method.Once)
    registration_method.SetOptimizerScalesFromPhysicalShift() 
    registration_method.SetInitialTransform(initial_transform, inPlace=False)
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas = [2,1,0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    final_transform = registration_method.Execute(fixed_image, moving_image)
    print('Final metric value: {0}'.format(registration_method.GetMetricValue()))
    print('Optimizer\'s stopping condition, {0}'.format(registration_method.GetOptimizerStopConditionDescription()))
    return (final_transform, registration_method.GetMetricValue())

def evaluate_metric(current_rotation, tx, f_image, m_image):
    
    '''This function evaluates the metric value in a thread safe manner'''
    
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)
    registration_method.SetInterpolator(sitk.sitkLinear)
    current_transform = sitk.Euler3DTransform(tx)
    current_transform.SetRotation(*current_rotation)
    registration_method.SetInitialTransform(current_transform)
    res = registration_method.MetricEvaluate(f_image, m_image)
    return res

def RegisterResample_image(fixed_image, moving_image, modality, optimization_method):
    
    '''Returns the moving_image registered and resampled to the fixed_image space and the transformation.
       Parameters selection arguments:
           modality='same' or 'different'
           optimization method='Multires' or 'ExplorExploit' or 'MetricEvaluate' or 'Exhaustive'
           
       Returns a moving_resampled image and a final_transform transformation. Remember to save them both to files after this function is applyed.'''
    
    #Initial alignment
    initial_transform = sitk.CenteredTransformInitializer(fixed_image, 
                                                          moving_image, 
                                                          sitk.Euler3DTransform(), 
                                                          sitk.CenteredTransformInitializerFilter.GEOMETRY)    
    #Registration
    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    if modality == 'same':
        registration_method.SetMetricAsMeanSquares()
    else:
        registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)
    
    #Interpolator settings.
    registration_method.SetInterpolator(sitk.sitkLinear)
    
    if optimization_method == 'Multires':
        # Multi-resolution frmework
        registration_method.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=100, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
        registration_method.SetOptimizerScalesFromPhysicalShift()
        registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1]) # Setup for the multi-resolution framework.
        registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
        registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
        registration_method.SetInitialTransform(initial_transform, inPlace=False)
        final_transform = registration_method.Execute(sitk.Cast(fixed_image, sitk.sitkFloat32), 
                                                      sitk.Cast(moving_image, sitk.sitkFloat32))
    
    elif optimization_method == 'ExplorExploit':
        # Exploration step
        registration_method.SetOptimizerAsExhaustive(numberOfSteps=[0,1,1,0,0,0], stepLength = np.pi)
        registration_method.SetOptimizerScales([1,1,1,1,1,1])
        registration_method.SetInitialTransform(initial_transform, inPlace=True)        
        registration_method.AddCommand(sitk.sitkIterationEvent, lambda: iteration_observer(registration_method))
        registration_method.AddCommand(sitk.sitkStartEvent, start_observer)
        _ = registration_method.Execute(fixed_image, moving_image)
        # Exploitation step.                 
        metricvalue_parameters_list.sort(key=lambda x: x[0]) #Sort our list from most to least promising solutions (low to high metric values).
        k_most_promising = min(3, len(metricvalue_parameters_list)) # We exploit the k_most_promising parameter value settings.
        final_results = []
        for metricvalue, parameters in metricvalue_parameters_list[0:k_most_promising]:
            initial_transform.SetParameters(parameters)
            final_results.append(multires_registration(fixed_image, moving_image, initial_transform))
        final_transform, _ = min(final_results, key=lambda x: x[1])
   
    elif optimization_method == 'MetricEvaluate':
        all_orientations = {'x=0, y=0, z=180': (0.0,0.0,np.pi),
                            'x=0, y=180, z=0': (0.0,np.pi,0.0),
                            'x=0, y=180, z=180': (0.0,np.pi,np.pi)}    
        p = ThreadPool(len(all_orientations)+1)
        orientations_list = [(0,0,0)] + list(all_orientations.values())
        all_metric_values = p.map(partial(evaluate_metric, 
                                          tx = initial_transform, 
                                          f_image = fixed_image,
                                          m_image = moving_image),
                                          orientations_list)
        best_orientation = orientations_list[np.argmin(all_metric_values)]
        initial_transform.SetRotation(*best_orientation)
        final_transform,_ = multires_registration(fixed_image, moving_image, initial_transform)
        
    elif optimization_method == 'Exhaustive':
        registration_method.SetOptimizerAsExhaustive(numberOfSteps=[0,1,1,0,0,0], stepLength = np.pi)
        registration_method.SetOptimizerScales([1,1,1,1,1,1])
        registration_method.SetInitialTransform(initial_transform, inPlace=True) #Perform the registration in-place so that the initial_transform is modified.
        registration_method.Execute(fixed_image, moving_image)
        final_transform, _ = multires_registration(fixed_image, moving_image, initial_transform)
        
    moving_resampled = sitk.Resample(moving_image, fixed_image, final_transform, sitk.sitkLinear, 0.0, moving_image.GetPixelID())
    
    return moving_resampled, final_transform


def apply_tfm_to_image(input_image, reference_image, transform):
    
    '''Applies a predefined transform to an input_image and returns the transformed_image with the same resolution and in the same space as the reference_image'''
    
    transformed_image = sitk.Resample(input_image, reference_image, transform, sitk.sitkLinear, 0.0)
    return transformed_image


#%% Registration and resampling of rois
    
def Resample_roi(input_roi, reference_image):
    
    '''Returns the input roi resampled to the reference image space.
       Uses a NearesrNeighbour interpolator and erodes roi of margin kernel type: Ball, with kernel radius: 1
       Remember to write the output roi into an image file after applying this function'''
    
    resample = sitk.ResampleImageFilter()
    resample.SetReferenceImage(reference_image)
    resample.SetInterpolator(sitk.sitkNearestNeighbor)
    output_roi = resample.Execute(input_roi)
    output_roi = sitk.BinaryErode(output_roi) #erode rois of margin kernel type: Ball, kernel radius: 1
    return output_roi


def apply_tfm_to_roi(input_roi, reference_image, transform):
    
    '''Applies a predefined transform to an input_roi, erodes roi of margin kernel type: Ball, kernel radius: 1 
    and returns the transformed_roi with the same resolution and in the same space as the reference_image'''
    
    transformed_roi = sitk.Resample(input_roi, reference_image, transform, sitk.sitkNearestNeighbor, 0.0)
    transformed_roi = sitk.BinaryErode(transformed_roi) #erode rois of margin kernel type: Ball, kernel radius: 1
    return transformed_roi


#%% Generate additional rois

def flip(roi):
    
    '''Returns roi flipped wrt the x axis. E.g. used to create mirror image roi in contralateral part of the body.
    Remember to save the flipped roi file after applying this function.'''
    
    flipped = sitk.Flip(sitk.Cast(roi, roi.GetPixelID()), (True, False, False), False)
    flipped.SetDirection(roi.GetDirection())
    flipped.SetOrigin(roi.GetOrigin())
    return flipped
    

def generate_3Dcrescent_roi(original_2D_crescent):
    
    '''Returns a crescent-shaped roi by dilating the morphology of a 2D crescent roi with a ball kernel with a radius of 3 pixels.
    Remember to save the dilated roi file after applying this function.'''

    dilate = sitk.DilateObjectMorphologyImageFilter() #Dilate morphology of crescent roi with a ball kernel with a radius of 3 pixels
    dilate.SetKernelType(sitk.sitkBall)
    dilate.SetKernelRadius(3)
    dilated = dilate.Execute(original_2D_crescent)
    return dilated

def generate_mask(image, roi):
    
    '''Returns the masked image of the roi applied to the image.
    Remember to save the masked image file after applying this function.'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image, roi)
    return mask

def generate_thresholded_roi(mask, low_thr, high_thr):
    
    '''Returns a binary roi of all the pixels in the mask within the specified low_thr and high_thr intensity values as thr_roi.
    Remember to save the thresholded roi file after applying this function.'''
    
    thr_roi = sitk.BinaryThreshold(mask, lowerThreshold=low_thr, upperThreshold=high_thr, insideValue=1, outsideValue=0)
    thr_roi = sitk.Cast(thr_roi, sitk.sitkUInt16)
    return thr_roi  
    
def set_mask_value(image, mask, value):
    msk32 = sitk.Cast(mask, sitk.sitkFloat32)
    return sitk.Cast(sitk.Cast(image, sitk.sitkFloat32) *
                     sitk.InvertIntensity(msk32, maximum=1.0) + 
                     msk32*value, image.GetPixelID())

#%% Generate smoothed image by applying a Gaussian smooth image filter

def gauss_smooth(image, sigma):
    
    '''Returns the smoothed image generated from the application of a Gaussian smoothing filter of specified sigma (in units of mm) to the original image.
    Remember to save the smoothed image file after applying this function.'''
    
    smooth = sitk.SmoothingRecursiveGaussianImageFilter()
    smooth.SetSigma(sigma)
    smoothed = smooth.Execute(image)
    smoothed = sitk.Cast(smoothed, image.GetPixelID())
    return smoothed

#%% Statistics on an image

def getImageMean(image):
    
    '''Calculates and returns the mean intensity value of an image'''
    
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    mean = stats.GetMean()
    
    return mean


def getImageSTD(image):
    
    '''Calculates and returns the standard deviation of the mean intensity of an image'''
    
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    std = stats.GetSigma()
    
    return std

def getNonZeroImageMean(image):
    
    '''Calculates and returns the mean intensity value of the non-zero voxels in an image.'''
    
    nda = sitk.GetArrayFromImage(image)
    mean = nda[np.nonzero(nda)].mean()
    
    return mean

def getNonZeroImageSTD(image):
    
    '''Calculates and returns the standard deviation of the mean intensity value of the non-zero voxels in an image.'''
    
    nda = sitk.GetArrayFromImage(image)
    std = nda[np.nonzero(nda)].std()
    
    return std

def getNonZeroStats(image):
    
    '''Calculates general stats of non-zero values from the image and returns them in a dictionary.
    The returned dictionary contains the follwoing parameters:
        'Volume [mm3]'
        'Mean intensity [SUV]'
        'Std of Mean [SUV]'
        'Median intensity [SUV]'
        'Max intensity [SUV]' 
        'Min intensity [SUV]'
        '''
    nda = sitk.GetArrayFromImage(image)
    volume = np.count_nonzero(nda)*image.GetSpacing()[0]*image.GetSpacing()[1]*image.GetSpacing()[2]
    mean = np.mean(nda[np.nonzero(nda)])
    std = np.std(nda[np.nonzero(nda)])
    median = np.median(nda[np.nonzero(nda)])
    maxVal = np.amax(nda[np.nonzero(nda)])
    minVal = np.amin(nda[np.nonzero(nda)])
    return {'Volume [mm3]':volume, 'Mean intensity [SUV]':mean, 'Std of Mean [SUV]': std, 'Median intensity [SUV]': median, 'Max intensity [SUV]':maxVal, 'Min intensity [SUV]':minVal}


def Normalize(image):
    '''Calculates the mean and std of the image and then scales the intensities to mean 0 and std 1'''
    
    nda = sitk.GetArrayFromImage(image)
    mask = np.array(nda>0, dtype=int)
    mask = sitk.GetImageFromArray(mask)
    mean = np.mean(nda[np.nonzero(nda)])
    std = np.std(nda[np.nonzero(nda)])
    centred = nda - mean
    scaled = centred/std
    scaled = sitk.GetImageFromArray(scaled)
    scaled = generate_mask(scaled, mask)
    return scaled, mask

def NanZeros(image):
    
    '''This fucntion replaces with Nan all the zero values in the input image.'''
    
    nda = sitk.GetArrayFromImage(image)
    nda[nda==0]=np.nan
    new = sitk.GetImageFromArray(nda)    
    return new

#%% Get general stats from ROI applied on an image
        
def getStatsRoi(roi, image):
    
    '''Calculates general stats from the roi applied on the image and returns them in a dictionary.
    The returned dictionary contains the follwoing parameters:
        'Volume [mm3]'
        'Mean intensity [SUV]'
        'Std of Mean [SUV]'
        'Median intensity [SUV]'
        'Max intensity [SUV]' 
        'Min intensity [SUV]'
        '''
    
    stats = sitk.LabelIntensityStatisticsImageFilter()
    stats.Execute(roi, image)
    volume = stats.GetPhysicalSize(1)
    mean = stats.GetMean(1)
    std = stats.GetStandardDeviation(1)
    median = stats.GetMedian(1)
    maxVal = stats.GetMaximum(1)
    minVal = stats.GetMinimum(1)
    return {'Volume [mm3]':volume, 'Mean intensity [SUV]':mean, 'Std of Mean [SUV]': std, 'Median intensity [SUV]': median, 'Max intensity [SUV]':maxVal, 'Min intensity [SUV]':minVal}

def getMaxRoi(roi, image):
    
    '''Calculates max intensity value from the roi applied on the image and returns it as a float'''
    
    stats = sitk.LabelIntensityStatisticsImageFilter()
    stats.Execute(roi, image)
    maxVal = stats.GetMaximum(1)
    return maxVal

def getMeanRoi(roi, image):
    
    '''Calculates mean intensity value from the roi applied on the image and returns it as a float'''
    
    stats = sitk.LabelIntensityStatisticsImageFilter()
    stats.Execute(roi, image)
    meanVal = stats.GetMean(1)
    return meanVal

def TBR(dictionary_stats_roi_on_image, dictionary_stats_ctrl_roi_on_image):
    
    '''Returns the TBRmean, stdTBRmean and TBRmax of a roi applied to an image, provided that there exists a dictionary containing the stats
    from the roi applied on the image calculated according to the getStatsRoi function for both the roi and it's cotralteral roi.'''
    
    TBRmean = dictionary_stats_roi_on_image.get('Mean intensity [SUV]')/dictionary_stats_ctrl_roi_on_image.get('Mean intensity [SUV]')
    stdTBRmean = (dictionary_stats_roi_on_image.get('Std of Mean [SUV]')/dictionary_stats_roi_on_image.get('Mean intensity [SUV]'))+(dictionary_stats_ctrl_roi_on_image.get('Std of Mean [SUV]')/dictionary_stats_ctrl_roi_on_image.get('Mean intensity [SUV]'))
    TBRmax = dictionary_stats_roi_on_image.get('Max intensity [SUV]')/dictionary_stats_ctrl_roi_on_image.get('Max intensity [SUV]')
    return TBRmean, stdTBRmean, TBRmax


def rois_volume_ratio(dictionary_stats_roi_on_image1, dictionary_stats_roi_on_image2):
    
    '''Returns the ratio between the volumes of roi 1 and roi 2, provided that there exists a dictionary containing the stats
    from the roi applied on an image calculated according to the getStatsRoi function for both the rois.'''
    
    RoisVolRatio = dictionary_stats_roi_on_image1.get('Volume [mm3]')/dictionary_stats_roi_on_image2.get('Volume [mm3]')
    return RoisVolRatio


#
    
def allPixInt(image, roi, maskName, temp):
    
    '''Returns a flatten array (in 2D) of the intensity values from all pixels in roi applied to an image.
    The process requires the generation of a mask (hence a maskname) and its writing into and reading from a temp folder (provide fullpath to temp folder).'''
    
    masking = sitk.MaskImageFilter() 
    mask = masking.Execute(image,roi) #Generate a mask of the roi applied to image
    sitk.WriteImage(mask, temp +'/'+ maskName +'.nii') #save image
    masked = sitk.ReadImage(temp +'/'+ maskName +'.nii') #read image
    array = sitk.GetArrayFromImage(masked) #convert image to array
    flat = array.flatten() #flatten array to 1D
    list = flat[flat != 0]
#    plt.hist(list)
    return list

#%%Segmentation

#Get intensity and coordinates of voxel with highest intensity in an image

def getMaxVox(image):
    
    '''Returns a float with max value of intensity and a tuple of three int values representing the coordinates x, y, z of the voxel with highest intensity in the image'''
    
    MinMax = sitk.MinimumMaximumImageFilter()
    MinMax.Execute(image)
    maxVal = MinMax.GetMaximum()
    nda = sitk.GetArrayFromImage(image)
    z, y, x = np.where(np.isclose(nda, maxVal)) # when comparing floating-point arrays    
    seed = (int(x[0]+1), int(y[0]+1), int(z[0]+1))
    return maxVal, seed


#Find local maxima in an image and return intensity and coordinates

def getLocalMaxVox(image, img_thresh):
    
    '''Returns a list of [float with max value of intensity and a tuple of three int values representing the coordinates x, y, z] for each local maxima in the image'''
    
    img = sitk.GetArrayFromImage(image)
    
    # Get local maximum values of desired neighborhood
    # I'll be looking in a 50x50x50 area
    img2 = ndimage.maximum_filter(img, size=(20, 20, 20))
    
    # Threshold the image to find locations of interest
    # I'm assuming 6 standard deviations above the mean for the threshold
    # img_thresh = img2.mean() + img2.std() * 6
    
    # Comparison between img2 and img to find the coordinates of local maxima
    coordinates = peak_local_max(img, min_distance=20, threshold_abs=img_thresh)
    
    coords = []
    for i in coordinates:
        z=int(i[0])
        y=int(i[1])
        x=int(i[2])
        coord_vox = (x, y, z)
        coords.append(coord_vox)

    return coords

#%%Transformations

# This function is from https://github.com/rock-learning/pytransform3d/blob/7589e083a50597a75b12d745ebacaa7cc056cfbd/pytransform3d/rotations.py#L302
def matrix_from_axis_angle(a):
    """ Compute rotation matrix from axis-angle.
    This is called exponential map or Rodrigues' formula.
    Parameters
    ----------
    a : array-like, shape (4,)
        Axis of rotation and rotation angle: (x, y, z, angle)
    Returns
    -------
    R : array-like, shape (3, 3)
        Rotation matrix
    """
    ux, uy, uz, theta = a
    c = np.cos(theta)
    s = np.sin(theta)
    ci = 1.0 - c
    R = np.array([[ci * ux * ux + c,
                   ci * ux * uy - uz * s,
                   ci * ux * uz + uy * s],
                  [ci * uy * ux + uz * s,
                   ci * uy * uy + c,
                   ci * uy * uz - ux * s],
                  [ci * uz * ux - uy * s,
                   ci * uz * uy + ux * s,
                   ci * uz * uz + c],
                  ])

    # This is equivalent to
    # R = (np.eye(3) * np.cos(theta) +
    #      (1.0 - np.cos(theta)) * a[:3, np.newaxis].dot(a[np.newaxis, :3]) +
    #      cross_product_matrix(a[:3]) * np.sin(theta))

    return R


def resample(image, transform):
    """
    This function resamples (updates) an image using a specified transform
    :param image: The sitk image we are trying to transform
    :param transform: An sitk transform (ex. resizing, rotation, etc.
    :return: The transformed sitk image
    """
    reference_image = image
    interpolator = sitk.sitkLinear
    default_value = 0
    return sitk.Resample(image, reference_image, transform,
                         interpolator, default_value)


def get_center(img):
    """
    This function returns the physical center point of a 3d sitk image
    :param img: The sitk image we are trying to find the center of
    :return: The physical center point of the image
    """
    width, height, depth = img.GetSize()
    return img.TransformIndexToPhysicalPoint((int(np.ceil(width/2)),
                                              int(np.ceil(height/2)),
                                              int(np.ceil(depth/2))))


def rotation3d(image, theta_z, show=False):
    """
    This function rotates an image across each of the x, y, z axes by theta_x, theta_y, and theta_z degrees
    respectively
    :param image: An sitk MRI image
    :param theta_x: The amount of degrees the user wants the image rotated around the x axis
    :param theta_y: The amount of degrees the user wants the image rotated around the y axis
    :param theta_z: The amount of degrees the user wants the image rotated around the z axis
    :param show: Boolean, whether or not the user wants to see the result of the rotation
    :return: The rotated image
    """
    theta_z = np.deg2rad(theta_z)
    euler_transform = sitk.Euler3DTransform()
    print(euler_transform.GetMatrix())
    image_center = get_center(image)
    euler_transform.SetCenter(image_center)

    direction = image.GetDirection()
    axis_angle = (direction[2], direction[5], direction[8], theta_z)
    np_rot_mat = matrix_from_axis_angle(axis_angle)
    euler_transform.SetMatrix(np_rot_mat.flatten().tolist())
    resampled_image = resample(image, euler_transform)
    if show:
        slice_num = int(input("Enter the index of the slice you would like to see"))
        plt.imshow(sitk.GetArrayFromImage(resampled_image)[slice_num])
        plt.show()
    return resampled_image


def affine_translate(image, dimension, x_translation=0.0, y_translation=0.0, z_translation=0.0):
    new_transform = sitk.AffineTransform(dimension)
    new_transform.SetTranslation((x_translation, y_translation, z_translation))
    resampled = resample(image, new_transform)
    return resampled

def set_mask_value(image, mask, value):    
    '''This function set the intensity of every voxel of the image within the specified mask to a set value.'''
    msk32 = sitk.Cast(mask, sitk.sitkFloat32)
    return sitk.Cast(sitk.Cast(image, sitk.sitkFloat32) *
                     sitk.InvertIntensity(msk32, maximum=1.0) + 
                     msk32*value, image.GetPixelID())

def point2str(point, precision=1):
    """
    Format a point for printing, based on specified precision with trailing zeros. Uniform printing for vector-like data 
    (tuple, numpy array, list).
    
    Args:
        point (vector-like): nD point with floating point coordinates.
        precision (int): Number of digits after the decimal point.
    Return:
        String represntation of the given point "xx.xxx yy.yyy zz.zzz...".
    """
    return ' '.join(format(c, '.{0}f'.format(precision)) for c in point)

# def translate_point(point, translation_vector):
#     '''Translate a point in 3D space given a point as a tuple of x, y, z cooridnates and a translation_vector as a tuple.
#     Returns the translated point coordinates as a tuple'''
#     dimension = 3        
#     offset = translation_vector # offset can be any vector-like data  
#     translation = sitk.TranslationTransform(dimension, offset)
#     transformed_point = (translation.TransformPoint(point))
#     return transformed_point

# def rotation3d_point(point, image, theta, x=0, y=0, z=0):
#     """
#     This function rotates a point in 3D with rotation around x, y, z axis of theta angle as a versor around a fixed center defined by the image.
#     :return: The rotated point
#     """
#     # theta = np.deg2rad(theta)   
#     image_center = (image.GetSize()[0]/2, image.GetSize()[1]/2, image.GetSize()[2]/2)
#     rotation = sitk.VersorTransform([x,y,z,theta], image_center)
#     rotated_point=(rotation.TransformPoint(point))
#     point2=(round(rotated_point[0]), round(rotated_point[1]), round(rotated_point[2]))
#     return point2

def translate_point(point, dimension, x_translation=0.0, y_translation=0.0, z_translation=0.0):
    '''Translate a point in 3D space given a point as a tuple of x, y, z cooridnates and a translation_vector as a tuple.
    Returns the translated point coordinates as a tuple'''
         
    t =(x_translation,y_translation,z_translation) 
    translation = sitk.TranslationTransform(dimension, t)

    # Only need to copy the translational component.
    rigid_euler = sitk.Euler3DTransform()
    rigid_euler.SetTranslation(translation.GetOffset())
    translated_point = (rigid_euler.TransformPoint(point))
    # new_transform = sitk.AffineTransform(dimension)
    # new_transform.SetTranslation((x_translation, y_translation, z_translation))
    # translated_point = (new_transform.TransformPoint(point))
    return translated_point


def rotation3d_point(point, image, theta_z):
    """
    This function rotates a point in 3D with rotation around x, y, z axis of theta angle as a versor around a fixed center defined by the image.
    :return: The rotated point
    """
    # theta = np.deg2rad(theta)   
    theta_z = np.deg2rad(theta_z)
    euler_transform = sitk.Euler3DTransform()
    print(euler_transform.GetMatrix())
    image_center = get_center(image)
    euler_transform.SetCenter(image_center)

    direction = image.GetDirection()
    axis_angle = (direction[2], direction[5], direction[8], theta_z)
    np_rot_mat = matrix_from_axis_angle(axis_angle)
    euler_transform.SetMatrix(np_rot_mat.flatten().tolist())
    rotated_point = (euler_transform.TransformPoint(point))
    point2=(round(rotated_point[0]), round(rotated_point[1]), round(rotated_point[2]))
    return point2

#%%Radiotherapy dose painting analysis functions

def reassign_voxel_intensity(image, mask, low_thr, high_thr, new_value):
        
    '''This function set the intensity of all the voxels belonging to the mask with original intensity within the low_thr and high_thr to the new_value'''
    masked_img = generate_mask(image, mask)
    mask0 =  generate_thresholded_roi(masked_img, low_thr, high_thr)
    masked_img = set_mask_value(masked_img, mask0, new_value)
    masked_img = generate_mask(masked_img, mask)
    return masked_img


def calculate_tcp(weighted_dose, alpha, C, tumour_presence_prob=None):
    
    '''This function returns the tumour control probability in a region of interest, from the dose prescription and the voxel-wise probability of tumour presence.
    C is a constrant (determined from radiobiological parameters) and alpha is a radiobiological parameter.'''
    
    if tumour_presence_prob is None:
        # default assumption is that there is tumour everywhere
        tumour_presence_prob = np.ones(weighted_dose.shape)

    ln_tcp = C * np.sum(tumour_presence_prob * np.exp(-alpha * weighted_dose))
    tcp = np.exp(-1*ln_tcp)
    return tcp

def generate_QF_map(Prescribed_dose_image, Planned_dose_image, roi):
    '''This function generates Quality Factor images between the planned and the prescribed dose within a roi. The closer QF is to 100, the better the quality of the plan wrt to the dose prescription.'''
    
    Rel_dose = abs((Planned_dose_image-Prescribed_dose_image)/Prescribed_dose_image)
    QF_img = 100 - Rel_dose*100
    QF_img = generate_mask(QF_img, roi)
    return QF_img