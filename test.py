import os
import numpy as np 
from datetime import datetime, date
import shutil 


def get_dates(base_path, PatientID):
    '''
    This function assumes that each image is stored in the folder w/ naming convention:
    CBCT_yyyymmdd - for CBCT images
    CT_yyyymmdd - for CT images 
    This function only lists the dates of the CBCT images.
    returns: yyyymmdd
    '''
    # function to find folder names of the CBCTs in the patient directory
    directory = base_path + '/' + PatientID + '/'
    ignore = ['CBCT_GROUPWISE', 'pCT', 'CBCT_pCT', 'T_model.nii.gz', 'MASKED_pCT.nii.gz']
    dates = [CBCT.replace('CBCT_', '') for CBCT in os.listdir(directory) if CBCT not in ignore]
    CBCT_dates = [date for date in list(dates) if date[0:2] != 'CT']

    return CBCT_dates


def get_CBCT_relative_timepoints(CBCT_dates):
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

path = 'D:/CUH_HN/NIFTI/'
ignore = ['CUH-UCL-02-001', 'CUH-UCL-02-002', 'CUH-UCL-02-004', 'CUH-UCL-02-005']
patients = [folder for folder in os.listdir(path) if folder.startswith('CUH-UCL')]

CBCT_directories_path = open('C:/Users/poppy/Documents/CUH_Data/Longitudinal/CBCT_directories.txt', 'a')
cpp_CBCT_directories_path = open('C:/Users/poppy/Documents/CUH_Data/Longitudinal/cpp_CBCT_directories.txt', 'a')
CT_directories_path = open('C:/Users/poppy/Documents/CUH_Data/Longitudinal/CT_directories.txt', 'a')
DEF_CBCT_directories_path = open('C:/Users/poppy/Documents/CUH_Data/Longitudinal/DEF_CBCT_directories.txt', 'a')

for patient in patients:

    CBCT_dates = get_dates(path, patient)
    CBCT_timepoints = get_CBCT_relative_timepoints(CBCT_dates)

    for CBCT_timepoint in CBCT_timepoints:
         
        CBCT_directories_path.write('/cluster/project7/HN_RT/CUH_HN/UCLHMODELSPACE_REGS/'+str(patient)+'/CBCT_'+str(CBCT_timepoint)+'/MASKED_CBCT.nii.gz \n')
        cpp_CBCT_directories_path.write('/cluster/project7/HN_RT/CUH_HN/UCLHMODELSPACE_REGS/'+str(patient)+'/CBCT_'+str(CBCT_timepoint)+'/cpp_CBCT.nii.gz \n')
        CT_directories_path.write('/cluster/project7/HN_RT/CUH_HN/UCLHMODELSPACE_REGS/'+str(patient)+'/CBCT_'+str(CBCT_timepoint)+'/MASKED_CT.nii.gz \n')
        DEF_CBCT_directories_path.write('/cluster/project7/HN_RT/CUH_HN/UCLHMODELSPACE_REGS/'+str(patient)+'/CBCT_'+str(CBCT_timepoint)+'/DEF_CBCT.nii.gz \n')


