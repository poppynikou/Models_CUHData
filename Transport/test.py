import os 
import numpy as np 


atlas = 'D:/CUH_HN/NIFTI/average_pCT.nii.gz'

def velocity_to_deformationfield(ref_img, input_transformation, output_transformation):

    reg_transform = 'C:/Users/poppy/Documents/Nifty/niftyreg_install/bin/reg_transform.exe'
    #print(ref_img)
    #print(input_transformation)
    #print(output_transformation)
    command = reg_transform + ' -ref ' + ref_img + ' -def ' + input_transformation + ' ' + output_transformation
    os.system(command)


patients = ['CUH-UCL-02-001','CUH-UCL-02-002','CUH-UCL-02-004','CUH-UCL-02-005']
timepoints = np.arange(0,39)

for patient in patients:
    for timepoint in timepoints:

        input_path = 'D:/CUH_HN/NIFTI/PSM_noLOO_T_All_CPS_4_ModelSpace/' + str(patient) + '/cpp_' + str(timepoint) + '.nii.gz'
        output_transformation = 'D:/CUH_HN/NIFTI/PSM_noLOO_T_All_CPS_4_ModelSpace/' + str(patient) + '/dff_' + str(timepoint) + '.nii.gz'
            

        velocity_to_deformationfield(atlas, input_path, output_transformation)