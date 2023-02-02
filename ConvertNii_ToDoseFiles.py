# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:58:57 2021

@author: cbri3325
"""



import math
import re
import struct
import numpy as np
import os
import time
import pydicom
import SimpleITK as sitk
from pathlib import Path



GTransferSyntaxUID = "1.2.840.10008.1.2"
GImplementationClassUID = "1.2.826.0.1.3680043.8.498.75006884747854523615841001"

RTDOSEModality = "RTDOSE"
RTPLANModality = "RTPLAN"
RTSTRUCTModality = "RTSTRUCT"

RTDoseSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.2"
RTStructSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.3"
RTPlanSOPClassUID = "1.2.840.10008.5.1.4.1.1.481.5"

Manufacturer = "RayStation"


def convert_nii_to_dicom_RTdosefile(nii_image, dcm_dir, output_directory=".", out_filename="dose.dcm"):
    """Converts a NII image to a Dicom RT dose file
    Args:
        nii_image: the nii image given as a simple ITK image
        dcm_dir: the directory containing the DICOM series
        output_directory (str, optional): The directory in which to place the generated Dicom
                                          files. Defaults to ".".
        out_filename: the name of the dicom RT dose file ending in .dcm
    """

    
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(dcm_dir)
    reader.SetFileNames(dicom_names)
    image = reader.Execute()
    
    dcm_dir_pat = Path(dcm_dir)
    
    first_slice_loc = 100000
    first_slice_obj = None
    for dcm_file in dcm_dir_pat.glob('*.dcm'):
            dicom_object = pydicom.read_file(dcm_file)
            if dicom_object.ImagePositionPatient[2] < first_slice_loc:
                first_slice_loc = dicom_object.ImagePositionPatient[2]
                first_slice_obj = dicom_object

    dicom_object = first_slice_obj
    x=dicom_object.RealWorldValueMappingSequence[0]
    
    # image= sitk.Cast(image, sitk.sitkFloat32)
    # image = image*x.RealWorldValueSlope + x.RealWorldValueIntercept
    # sitk.WriteImage(image, output_directory+'/image.nii')

    # Generate some new UIDs
    doseInstanceUID = pydicom.uid.generate_uid()
    
    # Populate required values for file meta information
    file_meta = pydicom.dataset.Dataset()
    file_meta.MediaStorageSOPClassUID = RTDoseSOPClassUID
    file_meta.TransferSyntaxUID = GTransferSyntaxUID
    file_meta.MediaStorageSOPInstanceUID = doseInstanceUID
    file_meta.ImplementationClassUID = GImplementationClassUID
    
    # Create the pydicom.dataset.FileDataset instance (initially no data elements, but file_meta supplied)
    RDfilename = f"RD.{file_meta.MediaStorageSOPInstanceUID}.dcm"
    ds = pydicom.dataset.FileDataset(
        RDfilename, {}, file_meta=file_meta, preamble=b"\x00" * 128
    )
    ds.InstanceCreationDate = time.strftime("%Y%m%d")
    ds.InstanceCreationTime = time.strftime("%H%M%S")

    ds.SOPClassUID = RTDoseSOPClassUID  # RT Dose Storage
    ds.SOPInstanceUID = doseInstanceUID
    ds.StudyDate = dicom_object.StudyDate
    ds.StudyTime = dicom_object.StudyTime
    ds.AccessionNumber = ""
    ds.Modality = "RTDOSE"
    ds.Manufacturer = dicom_object.Manufacturer
    ds.ReferringPhysicianName = dicom_object.ReferringPhysicianName
    ds.StationName = dicom_object.StationName
    ds.StudyDescription = dicom_object.StudyDescription
    ds.SeriesDescription = dicom_object.SeriesDescription
    ds.ManufacturerModelName = dicom_object.ManufacturerModelName
    ds.PatientName = dicom_object.PatientName
    ds.PatientID = dicom_object.PatientID
    ds.PatientBirthDate = dicom_object.PatientBirthDate
    ds.PatientSex = dicom_object.PatientSex
    ds.SliceThickness = dicom_object.SliceThickness
    # ds.DeviceSerialNumber = dicom_object.DeviceSerialNumber
    # ds.SoftwareVersion = dicom_object.SoftwareVersion
    ds.StudyInstanceUID = dicom_object.StudyInstanceUID
    ds.SeriesInstanceUID = dicom_object.SeriesInstanceUID
    ds.StudyID = dicom_object.StudyID
    ds.SeriesNumber = dicom_object.SeriesNumber
    ds.InstanceNumber = dicom_object.InstanceNumber
    ds.ImagePositionPatient = dicom_object.ImagePositionPatient
    ds.ImageOrientationPatient = dicom_object.ImageOrientationPatient
    ds.FrameOfReferenceUID = dicom_object.FrameOfReferenceUID
    ds.PositionReferenceIndicator = dicom_object.PositionReferenceIndicator
    ds.SamplesPerPixel = dicom_object.SamplesPerPixel
    ds.PhotometricInterpretation = dicom_object.PhotometricInterpretation
    ds.NumberOfFrames = int(nii_image.GetSize()[2])
    ds.FrameIncrementPointer = pydicom.dataelem.Tag("GridFrameOffsetVector")
    ds.Rows = dicom_object.Rows
    ds.Columns = dicom_object.Columns
    ds.PixelSpacing = dicom_object.PixelSpacing
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0
    ds.DoseUnits = 'GY'
    ds.DoseType = 'PHYSICAL'
    ds.DoseSummationType = 'PLAN'
    slice_thickness = image.GetSpacing()[2]
    ds.GridFrameOffsetVector = [ds.ImagePositionPatient[2] + x*slice_thickness for x in range(ds.NumberOfFrames)]

    np_dose = sitk.GetArrayFromImage(nii_image)
    max_dose_val = np_dose.max()

    ds.DoseGridScaling = max_dose_val/65535 #Divide maximum dose value of the image by 2^16
    # print(ds.DoseGridScaling)
    np_dose_scaled = np_dose/ds.DoseGridScaling
    np_dose_scaled = np_dose_scaled.astype(np.uint16)
    ds.TissueHeterogeneityCorrection = "IMAGE"
    # ds.ReferencedRTPlanSequence = dicom_object.ReferencedRTPlanSequence #need to link to RayStation frame of reference
    # ds.ReferencedSOPClassUID = dicom_object.ReferencedSOPClassUID
    
    #Need to copy intensity voxel values to pixels into dose image
    
    ds.PixelData = np_dose_scaled.tobytes()
 
    # Save the RTDose Dicom File
    output_file = os.path.join(output_directory, out_filename)
    ds.save_as(output_file)

if __name__ == "__main__":
    convert_nii_to_dicom_RTdosefile("./Data")