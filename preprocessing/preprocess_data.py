import os
from classes import *
from natsort import natsorted

base_path = 'D:/CUH_HN/NIFTI/'

ignore_patients = np.arange(0,3)
ignore = ['CUH-UCL-02-0' + str("%.2d" % i) for i in ignore_patients]

patients = [ 'CUH-UCL-02-018', 'CUH-UCL-02-029']

#patient_list = [15]
#patients = ['CUH-UCL-02-0' + str("%.2d" % i) for i in patient_list]
#CBCT_list = [[file[15:28] for file in os.listdir('D:/CUH_HN/2D_slices') if file[0:14] == patient] for patient in patients]
#print(CBCT_list)

# patients = ['CUH-UCL-02-0' + str("%.2d" % i) for i in patient_list]
#[file for file in os.listdir('D:/CUH_HN/NIFTI') if file.__contains__('CUH-UCL')]

flip_record = base_path + 'flip_record.txt'
flip_record = open(flip_record, 'a')
flip_record.write('Patient No_Left_CTV_Voxels No_Right_CTV_Voxels Difference Flip \n')

for p_index, patient in enumerate(patients):
    

    # check if the preprocessing for that patient has already been done   
    PatientID = patient
    
    PatientObj = PatientData(PatientID, base_path)
    
    CBCT_dates = PatientObj.get_dates()
    pCT_date = PatientObj.get_pCT_date()
    
    # create image class 
    ImgObj = Image(PatientID, base_path, pCT_date)

    # find patient folder to search through
    patient_path = PatientObj.get_patient_folder()

    flip_img_Bool = ImgObj.flip_img_Bool(flip_record)
    
    
    # preprocess the images 
    for path, subdirs, files in os.walk(patient_path):
        
        for index, name in enumerate(files):
                
                file_path = path + '/' + name 
                
                
                
                if name == 'CBCT.nii.gz':
                    if flip_img_Bool:
                        print(patient, flip_img_Bool)
                        ImgObj.flip_img(file_path, new_img_path=True)
                        ImgObj.rename_parotid(path, name, flip = True)
                    else:
                        ImgObj.rename_parotid(path, name, flip = False)
                '''
                
                if name == 'CT.nii.gz':
                    print(patient, flip_img_Bool)
                    ImgObj.convert_to_float(file_path)
                    ImgObj.clip_HU(file_path)
                    ImgObj.mask_CT(file_path)
                    ImgObj.mask_anatomy()
                
                '''


                if name == 'CBCT.nii.gz':
                    ImgObj.convert_to_float(file_path)
                    ImgObj.mask_CBCT(file_path, 4, mask_couch = True)
            
                        
                    
