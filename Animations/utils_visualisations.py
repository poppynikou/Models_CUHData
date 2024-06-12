import nibabel as nib
import numpy as np 
import os 
from datetime import datetime, date


def get_pixel_spacing(CT_img):
    
    return CT_img.header['pixdim'][1:3]

def get_slice_thickness(CT_img):

    return CT_img.header['pixdim'][3]



def clip_intensities(niftiImgIn, minInt=-1000, maxInt=1000):
    """
    Clip intensities of nifti image to set range
    """
    img = niftiImgIn.get_fdata()
    affine = niftiImgIn.affine
    hdr = niftiImgIn.header

    img[np.isnan(img)] = -1000

    img[img <= minInt] = minInt
    img[img >= maxInt] = maxInt
    niftiImgOut = nib.Nifti1Image(img, affine, hdr)
    return niftiImgOut


def rescale_intensities(niftiImgIn, newMinInt=0, newMaxInt=1):
    """
    Rescale intensities of nifti image between set boundaries
    """
    img = niftiImgIn.get_fdata()
    affine = niftiImgIn.affine
    hdr = niftiImgIn.header

    minIn = np.amin(img)
    maxIn = np.amax(img)

    rescaledImg = (((img - minIn)/(maxIn - minIn)) * (newMaxInt - newMinInt)) + newMinInt
    rescaledImg = rescaledImg * 255
    niftiImgOut = nib.Nifti1Image(rescaledImg, affine, hdr)
    return niftiImgOut

def wrapper_function(niftiImgIn, minInt=-1000, maxInt=1000, newMinInt=0, newMaxInt=1):
    """
    Wraps around the two image scaling functions
    """
    niftiImgOut = clip_intensities(niftiImgIn, minInt=-1000, maxInt=1000)

    RescaledNiftiImgOut = rescale_intensities(niftiImgOut, newMinInt=0, newMaxInt=1)

    return RescaledNiftiImgOut


def get_dates(patient_no):
    '''
    This function assumes that each image is stored in the folder w/ naming convention:

    CBCT_yyyymmdd - for CBCT images
    CT_yyyymmdd - for CT images 

    This function only lists the dates of the CBCT images.
    returns: yyyymmdd
    
    '''
    # function to find folder names of the CBCTs in the patient directory
    directory = 'D://UCLH_HN//HN_' + str(patient_no) + '//'
    dates = os.listdir(directory)
    CBCT_dates = [date[-8:]  for date in list(dates) if date[0:4] == 'CBCT' and date[-3:] != 'dcm']
    return CBCT_dates


def get_time_points(CBCT_dates):
    '''
    This assumes that the first time points is the first CBCT in the series
    '''

    # create an array to store the time points
    time_points = []

    # find the date for the first CBCT date 
    day_zero = CBCT_dates[0]
    date_zero = date(int(day_zero[0:4]), int(day_zero[4:6]),  int(day_zero[6:8]))

    # loop through dates and get the time points
    for CBCT_date in CBCT_dates:
        
        CBCT_date = CBCT_date
        date_i = date(int(CBCT_date[0:4]), int(CBCT_date[4:6]),  int(CBCT_date[6:8]))
        time_points.append((date_i-date_zero).days) 

    return time_points

def get_pCT_date(patient_no):
    # function to find the pCT date for a specific patient 
    directory = 'D://UCLH_HN//HN_' + str(patient_no) + '//'
    # list all contents of the directory 
    dates = os.listdir(directory)
    # find all folders which start with 'CT_'
    CT_dates = [date[-8:]  for date in list(dates) if date[0:2] == 'CT' and date[-3:] != 'dcm']
    CT_dates = np.asarray(CT_dates, dtype=object)
    # get the dates times, so that you can get the one which is earliest in time 
    CT_datetimes = [datetime(int(CT[0:4]), int(CT[4:6]), int(CT[6:8])) for CT in CT_dates]
    pCT_datetime = [min(CT_datetimes)]

    Boolean = np.in1d(CT_datetimes, pCT_datetime)

    pCT_date = CT_dates[Boolean][0]
    
    return pCT_date

def velocity_to_deformationfield(ref_img, input_transformation, output_transformation):
    reg_transform = 'C:/Users/poppy/Documents/Nifty/niftyreg_install/bin/reg_transform.exe'
    command = reg_transform + ' -ref ' + ref_img + ' -def ' + input_transformation + ' ' + output_transformation
    os.system(command)

    return 

def compose_transform(ref_img, composed_transformations, transformation2, transformation1):
    '''
    Compose two transformations of any recognised type and returns a deformation field.
    composed_transformations(x) = transformation2(transformation1(x)).
    '''
    reg_transform = 'C:/Users/poppy/Documents/Nifty/niftyreg_install/bin/reg_transform.exe'
    command = reg_transform + ' -ref ' + ref_img + ' -comp ' + transformation1 + ' ' + transformation2 + ' ' + composed_transformations

    os.system(command)

def Transport_Model_to_Patient(input_path, output_path, patient_to_modelspace_data, CBCT_img_modelspace, CBCT_img_patientspace):
    '''
    input_path: path to velocity field in model space
    output_path: path to deformation field in patient space
    patient_to_modelspace_data: folder in which transformations for deforming patient to model space data is stored (i.e. Transported patients)
    '''

    # convert velocity field to deformation field
    output_transformation = 'deformation_field.nii.gz'
    velocity_to_deformationfield(CBCT_img_modelspace, input_path, output_transformation)

    T2 = patient_to_modelspace_data + '/T2.nii.gz'
    T2_inv = patient_to_modelspace_data + '/inv_T2.nii.gz'
    
    intermediate_composed_transformations = 'intermediate_composed_transformations.nii.gz'
    compose_transform(CBCT_img_modelspace, intermediate_composed_transformations, output_transformation, T2_inv)

    compose_transform(CBCT_img_patientspace, output_path, T2, intermediate_composed_transformations)

    os.remove(intermediate_composed_transformations)
    os.remove(output_transformation)

def deform_CT(ref_img, float_img, transformation_path, resampled_img):
     # generic wrapper

    reg_resample = 'C:/Users/poppy/Documents/Nifty/niftyreg_install/bin/reg_resample.exe'

    command = reg_resample + ' -ref ' + ref_img + ' -flo ' + float_img + ' -trans ' + transformation_path + ' -res ' + resampled_img + ' -inter 3 -omp 12 -pad nan'
    os.system(command)
