import os
import numpy as np 
# NIFTY REG WRAPPED FUNCTIONS 

def UpdSform(reg_transform_path, img_to_be_updated_path, affine_matrix_path, updated_img_path):

    command = reg_transform_path +' -updSform ' + img_to_be_updated_path + ' ' + affine_matrix_path + ' ' + updated_img_path
    os.system(command)


def rigidReg(reg_aladin_path, ref_img, float_img, affine_matrix, resampled_img, RigOnly):

    basic_command = reg_aladin_path + ' -ref ' + ref_img + ' -flo ' + float_img + ' -aff ' + affine_matrix + ' -res ' + resampled_img + ' -omp 12 '
    if RigOnly:
        command = basic_command + '-rigOnly'
    else:
        command = basic_command
    os.system(command)
    
def rigidReg_SINGLETREAD(reg_aladin_path, ref_img, float_img, affine_matrix, resampled_img, RigOnly):
    
    basic_command = reg_aladin_path + ' -ref ' + ref_img + ' -flo ' + float_img + ' -aff ' + affine_matrix + ' -res ' + resampled_img + ' -omp 1 '
    if RigOnly:
        command = basic_command + '-rigOnly'
    else:
        command = basic_command
    os.system(command)
    
    

def avgAff(reg_average_path, average_affine_path, affine_matrixes):

    '''
    average_affine_path: where to store the average affine matrix
    affine_matrixes: a list of affine matrixes to average 
    '''

    command = reg_average_path + ' ' + average_affine_path + ' -avg ' 
    for i in np.arange(0,len(affine_matrixes)):

        command = command + affine_matrixes[i] + ' '

    os.system(command)
 
def invAff(reg_transform_path, ref_img, affine, inv_affine):

    command = reg_transform_path + ' -ref ' + ref_img + ' -invAff ' + affine + ' ' +  inv_affine
    os.system(command)


def compAff(reg_transform_path, ref_img, ref_affine, affine_matrix, comp_matrix):

    command = reg_transform_path +  ' -ref ' + ref_img + ' -comp ' + ref_affine + ' ' + affine_matrix + ' ' + comp_matrix
    os.system(command)

def resampleImg(reg_resample_path, ref_img, float_img, transformation, resampled_img):

    command = reg_resample_path + ' -ref ' + ref_img +  ' -flo ' + float_img + ' -trans ' + transformation + ' -res ' + resampled_img + ' -inter 3 -pad nan'
    os.system(command)
    
def resampleBINImg(reg_resample_path, ref_img, float_img, transformation, resampled_img):
    
    command = reg_resample_path + ' -ref ' + ref_img +  ' -flo ' + float_img + ' -trans ' + transformation + ' -res ' + resampled_img + ' -inter 1 -pad 0'
    os.system(command)

def deformableReg(reg_f3d_path, ref_img, float_img, resampled_img, transformation):

    command_basic = reg_f3d_path + ' -ref ' + ref_img + ' -flo ' + float_img + ' -res ' + resampled_img + ' -cpp ' + transformation + ' -omp 12 '
    command_params = ' -sx -10 -sy -10 -sz -10 -be 0 --lncc -5 -ln 5 -vel -le 0.1 -pad nan'
    
    command = command_basic + command_params
    os.system(command)


def RigidToDeformation(reg_transform_path, ref_img, input_transformation, output_transformation):

    command = reg_transform_path + ' -ref ' + ref_img + ' -def ' + input_transformation + ' ' + output_transformation 
    os.system(command)


def ComposeTransformations(reg_transform_path, ref_img, transformation1, transformation2, output_transformation):
    
    command = reg_transform_path + ' -ref ' + ref_img + ' -comp ' + transformation1 + ' ' + transformation2 + ' '  + output_transformation
    os.system(command)
    
    
