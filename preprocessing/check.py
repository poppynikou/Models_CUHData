import os
from natsort import natsorted
import nibabel as nib 
import numpy as np 

patients = natsorted(set([file[0:14] for file in os.listdir('D:/CUH_HN/NIFTI')]))

for p_index, patient in enumerate(patients):
    # preprocess the images 

    patient_path = 'D:/CUH_HN/NIFTI/' + str(patient)

    for path, subdirs, files in os.walk(patient_path):
    
        for index, name in enumerate(files):
            
            file_path = path + '/' + name 
            
            if name == 'MASKED_CBCT.nii.gz' or name == 'MASKED_CT.nii.gz':
                
                img = nib.load(file_path)
                img_obj = img.get_fdata()
                img_obj = np.array(img_obj)

                if img_obj.dtype != 'float64':
                    print(file_path)
                    print('not float')
                
               


