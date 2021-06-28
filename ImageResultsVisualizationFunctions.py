# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:57:21 2020

@author: s4451992
"""

#%% Set environmental variables

import numpy as np
import SimpleITK as sitk
import os
from PIL import Image, ImageEnhance

#%% Obtain an image intensity range - useful for windowing
    
def get_img_intensity_range(image):
    
    '''Prints the image intensity rang and returns a dictionary with values for 'Window Min' and 'Window Max'.'''
    
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    Max_int = stats.GetMaximum()
    Min_int = stats.GetMinimum()
    print("The image intensity range is: " + str(Min_int) + "-" + str(Max_int))
    return {'Window Max': Max_int, 'Window Min':Min_int}


def window_img(image, outMin=0.0, outMax=255.0):
    
    '''Use the image_range Min and Max intensity values to perform intensity windowing and map the intensity values to [0,255] and cast to 8-bit unsigned int.
    Returns the windowed image.'''
    
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    Max_int = stats.GetMaximum()
    Min_int = stats.GetMinimum()
    windowed = sitk.Cast(sitk.IntensityWindowing(image, windowMinimum=Min_int, windowMaximum=Max_int, 
                                             outputMinimum=outMin, outputMaximum=outMax), sitk.sitkUInt8)
    return windowed

 
def change_brightness(image_path, factor):
    
    '''This function allows to change the brightness of an image by a defined factor. 
    factor = 1 gives the original image
    factor < 1 darkens the image
    factor > 1 brightens the image'''
    
    #Open the image
    img = Image.open(image_path)
    
    #image brightness enhancer
    enhancer = ImageEnhance.Brightness(img)
    
    im_output = enhancer.enhance(factor)
    im_output.save(image_path[:-4]+'_bright.png')
    
    return im_output



def window_probMap(image):
    
    '''Use the image_range Min and Max intensity values to perform intensity windowing and map the intensity values to [0.00001,1] and cast to 8-bit unsigned int.
    Returns the windowed image.'''
    
    stats = sitk.StatisticsImageFilter()
    stats.Execute(image)
    Max_int = stats.GetMaximum()
    Min_int = stats.GetMinimum()
    windowed = sitk.IntensityWindowing(image, windowMinimum=Min_int, windowMaximum=Max_int, 
                                             outputMinimum=0.00001, outputMaximum=Max_int)
    return windowed

#%% Alpha blending
    
def simple_blend(background_path, overlay_path, alpha, outImg_name=None):
    
    '''This function generated a composite image given by the overlay of the overlay onto the background.'''
    
    background = Image.open(background_path)
    overlay = Image.open(overlay_path)

    new_img = Image.blend(background, overlay, alpha)
    
    if not outImg_name:
        new_img.save(overlay_path[:-4]+'_overlay.png',"PNG")
    else:
        new_img.save(overlay_path[:-4]+'_'+outImg_name+'.png',"PNG")
    return new_img



def mask_image_multiply(mask, image):
    components_per_pixel = image.GetNumberOfComponentsPerPixel()
    if  components_per_pixel == 1:
        return mask*image
    else:
        return sitk.Compose([mask*sitk.VectorIndexSelectionCast(image,channel) for channel in range(components_per_pixel)])
    
def alpha_blend(image1, image2, alpha = 0.5, mask1=None,  mask2=None):
    '''
    Alpha blend two images, pixels can be scalars or vectors.
    The region that is alpha blended is controlled by the given masks.
    '''
    
    if not mask1:
        mask1 = sitk.Image(image1.GetSize(), sitk.sitkFloat32) + 1.0
        mask1.CopyInformation(image1)
    else:
        mask1 = sitk.Cast(mask1, sitk.sitkFloat32)
    if not mask2:
        mask2 = sitk.Image(image2.GetSize(),sitk.sitkFloat32) + 1
        mask2.CopyInformation(image2)
    else:        
        mask2 = sitk.Cast(mask2, sitk.sitkFloat32)

    components_per_pixel = image1.GetNumberOfComponentsPerPixel()
    if components_per_pixel>1:
        img1 = sitk.Cast(image1, sitk.sitkVectorFloat32)
        img2 = sitk.Cast(image2, sitk.sitkVectorFloat32)
    else:
        img1 = sitk.Cast(image1, sitk.sitkFloat32)
        img2 = sitk.Cast(image2, sitk.sitkFloat32)
        
    intersection_mask = mask1*mask2
    
    intersection_image = mask_image_multiply(alpha*intersection_mask, img1) + \
                         mask_image_multiply((1-alpha)*intersection_mask, img2)
    return intersection_image + mask_image_multiply(mask2-intersection_mask, img2) + \
           mask_image_multiply(mask1-intersection_mask, img1)
           
def grayScale_to_ColorMap(image, z_slice, colorMap):
    
    
    '''This function returns the image of a coronal section taken at the specified z_slice as 
    a colorMap. The colorMap can be selected amongst the following options:
        Autumn, Blue, Cool, Copper, Green, Grey, HSV, Hot, Jet, OverUnder, Red, Spring, Summer, Winter.'''
    
    img_255 = window_img(image)
    img_coronal = img_255[:,:,z_slice]
    
    if colorMap == 'Jet':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Jet)
    elif colorMap == 'Autumn':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Autumn)
    elif colorMap == 'Blue':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Blue)
    elif colorMap == 'Cool':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Cool)
    elif colorMap == 'Copper':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Copper)
    elif colorMap == 'Green':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Green)
    elif colorMap == 'Grey':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Grey)
    elif colorMap == 'HSV':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.HSV)
    elif colorMap == 'Hot':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Hot)
    elif colorMap == 'OverUnder':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.OverUnder)
    elif colorMap == 'Red':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Red)
    elif colorMap == 'Spring':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Spring)
    elif colorMap == 'Summer':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Summer)
    elif colorMap == 'Winter':
        overlay_color_img = sitk.ScalarToRGBColormap(img_coronal, sitk.ScalarToRGBColormapImageFilter.Winter)
        
    return overlay_color_img
    
    
    
#%%Image overlay and roi overlay display

def colormap_imgs_overlay(image1, image2, z_slice, colorMap, alpha_value=0.5):
    
    '''This function returns the image of a coronal section taken at the specified z_slice of image 2 expressed as 
    a colorMap alpha blended on top of image 1 in grayscale. The colorMap can be selected amongst the following options:
        Autumn, Blue, Cool, Copper, Green, Grey, HSV, Hot, Jet, OverUnder, Red, Spring, Summer, Winter.'''
    
    img1_255 = window_img(image1)
    img2_255 = window_img(image2)
    img1_coronal = img1_255[:,:,z_slice]
    img2_coronal = img2_255[:,:,z_slice]
    
    if colorMap == 'Jet':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Jet)
    elif colorMap == 'Autumn':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Autumn)
    elif colorMap == 'Blue':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Blue)
    elif colorMap == 'Cool':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Cool)
    elif colorMap == 'Copper':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Copper)
    elif colorMap == 'Green':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Green)
    elif colorMap == 'Grey':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Grey)
    elif colorMap == 'HSV':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.HSV)
    elif colorMap == 'Hot':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Hot)
    elif colorMap == 'OverUnder':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.OverUnder)
    elif colorMap == 'Red':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Red)
    elif colorMap == 'Spring':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Spring)
    elif colorMap == 'Summer':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Summer)
    elif colorMap == 'Winter':
        overlay_color_img = sitk.ScalarToRGBColormap(img2_coronal, sitk.ScalarToRGBColormapImageFilter.Winter)
    
    combined_volume = sitk.Cast(alpha_blend(sitk.Compose(img1_coronal, img1_coronal, img1_coronal), overlay_color_img, alpha = alpha_value), sitk.sitkVectorUInt8)
    return combined_volume






def roi_to_img_overlay(image, roi, z_slice, color, opacity=0.2):
    
    '''This function returns the image of a coronal section taken at the specified z_slice of an roi alpha blended 
    at the given opacity value on top of the given image in grayscale. 
    The color of the roi can be selected amongst the following options:
        red, green, blue, pink, magenta, cyan, orange, yellow, aqua marine, sky blue, purple.
        The default value of opacity is 0.2.'''
    
    img_255 = window_img(image)
    roi_255 = window_img(roi)
    img_coronal = img_255[:,:,z_slice]
    roi_coronal = roi_255[:,:,z_slice]
    
    if color == 'red':
        color_roi = [255, 0, 0]
    elif color == 'green':
        color_roi = [0, 255, 0]
    elif color == 'blue':
        color_roi = [0, 0, 255]
    elif color == 'pink':
        color_roi = [255,105,180]
    elif color == 'gold':
        color_roi = [255,215,0]
    elif color == 'magenta':
        color_roi = [255,0,255]
    elif color == 'cyan':
        color_roi = [0,255,255]
    elif color == 'orange':
        color_roi = [255,165,0]
    elif color == 'yellow':
        color_roi = [255,255,0]
    elif color == 'aqua marine':
        color_roi = [127,255,212]
    elif color == 'sky blue':
        color_roi = [0,191,255]
    elif color == 'purple':
        color_roi = [148,0,211]
    
    coronal_combined = sitk.LabelOverlay(image=img_coronal, 
                                         labelImage=roi_coronal,
                                         opacity=opacity, backgroundValue=0.0, colormap=color_roi)

    return coronal_combined



def labels_to_img_overlay(image, labelImage, z_slice, opacity=0.2):
    
    '''This function returns the image of a coronal section taken at the specified z_slice of an labelImage 
    (containing more than one label) alpha blended 
    at the given opacity value on top of the given image in grayscale. 
    The color of the labels is arbitrarily assigned. 
    The default value of opacity is 0.2.'''
    
#    stats = sitk.LabelShapeStatisticsImageFilter()
#    stats.Execute(roi_coronal)
#    s = stats.GetNumberOfLabels()    
    changelabel = sitk.ChangeLabelImageFilter()
    changelabel.SetChangeMap({0:0, 1:2, 2:4, 3:6, 4:8, 5:10})
    labelImage = changelabel.Execute(labelImage)    
    img_255 = window_img(image)
    roi_255 = window_img(labelImage)
    img_coronal = img_255[:,:,z_slice]
    roi_coronal = roi_255[:,:,z_slice]    
    coronal_combined = sitk.LabelOverlay(image=img_coronal, labelImage=roi_coronal, opacity=opacity, backgroundValue=0.0)

    return coronal_combined



def roibound_to_img_overlay(image, roi, z_slice, color, opacity=1):
    
    '''This function returns the image of a coronal section taken at the specified z_slice of the roi boundaries 
    on top of the given image in grayscale. 
    The color of the roi can be selected amongst the following options:
        red, green, blue, pink, magenta, cyan, orange, yellow, aqua marine, sky blue, purple.
        The default value of opacity is 0.2.'''
    
    img_255 = window_img(image)
    roi_255 = window_img(roi)
    img_coronal = img_255[:,:,z_slice]
    roi_coronal = roi_255[:,:,z_slice]
    
    if color == 'red':
        color_roi = [255, 0, 0]
    elif color == 'green':
        color_roi = [0, 255, 0]
    elif color == 'blue':
        color_roi = [0, 0, 255]
    elif color == 'pink':
        color_roi = [255,105,180]
    elif color == 'gold':
        color_roi = [255,215,0]
    elif color == 'magenta':
        color_roi = [255,0,255]
    elif color == 'cyan':
        color_roi = [0,255,255]
    elif color == 'orange':
        color_roi = [255,165,0]
    elif color == 'yellow':
        color_roi = [255,255,0]
    elif color == 'aqua marine':
        color_roi = [127,255,212]
    elif color == 'sky blue':
        color_roi = [0,191,255]
    elif color == 'purple':
        color_roi = [148,0,211]
    
    
    contour_overlaid_image = sitk.LabelMapContourOverlay(sitk.Cast(roi_coronal, sitk.sitkLabelUInt8), 
                                                         img_coronal, 
                                                         opacity = opacity, 
                                                         contourThickness=[2,2],
                                                         dilationRadius= [3,3],
                                                         colormap=color_roi)
    
    return contour_overlaid_image





def roidiff_to_img_overlay(image, roi1, roi2, z_slice, color, opacity=0.2):
    
    '''This function returns the image of a coronal section taken at the specified z_slice of the difference of two rois 
    alpha blended at the given opacity value on top of the given image in grayscale. 
    The color of the final roi can be selected amongst the following options:
        red, green, blue, pink, magenta, cyan, orange, yellow, aqua marine, sky blue, purple.
    The default value of opacity is 0.2.'''
    
    img_255 = window_img(image)
    roi1_255 = window_img(roi1)
    roi2_255 = window_img(roi2)
    img_coronal = img_255[:,:,z_slice]
    roi1_coronal = roi1_255[:,:,z_slice]
    roi2_coronal = roi2_255[:,:,z_slice]
    
    diff_roi = (roi1_coronal!=roi2_coronal)
    
    
    if color == 'red':
        color_roi = [255, 0, 0]
    elif color == 'green':
        color_roi = [0, 255, 0]
    elif color == 'blue':
        color_roi = [0, 0, 255]
    elif color == 'pink':
        color_roi = [255,105,180]
    elif color == 'gold':
        color_roi = [255,215,0]
    elif color == 'magenta':
        color_roi = [255,0,255]
    elif color == 'cyan':
        color_roi = [0,255,255]
    elif color == 'orange':
        color_roi = [255,165,0]
    elif color == 'yellow':
        color_roi = [255,255,0]
    elif color == 'aqua marine':
        color_roi = [127,255,212]
    elif color == 'sky blue':
        color_roi = [0,191,255]
    elif color == 'purple':
        color_roi = [148,0,211]
    
    coronal_combined = sitk.LabelOverlay(image=img_coronal, 
                                         labelImage=diff_roi,
                                         opacity=opacity, backgroundValue=0.0, colormap=color_roi)
    return coronal_combined