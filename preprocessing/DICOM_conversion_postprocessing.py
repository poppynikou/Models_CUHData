'''
This is to combine ctv high, medium and low into three singular structures.
'''
import os
import nibabel as nib 
import numpy as np 


Base_path = 'D:/CUH_HN/DICOM/'

ignore_patients = []#]
ignore = ['CUH-UCL-02-0' + str("%.2d" % i) for i in ignore_patients]
patients = [str(folder) for folder in os.listdir(Base_path) if (folder.startswith('CUH-UCL')) and (folder not in ignore)]


for patient in patients:

    print(patient)

    patient_path = 'D:/CUH_HN/NIFTI/'+str(patient)

    CT_folders = [folder for folder in os.listdir(patient_path) if folder.startswith('CT_')]


    for CT_folder in CT_folders:

        CT_folder_path = patient_path + '/' + str(CT_folder) + '/CT.nii.gz'
        CT_obj = nib.load(CT_folder_path)
        CT_img = CT_obj.get_fdata()
        CT_shape = np.shape(CT_img)

        structures_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/'

        CTV_structures = [structure for structure in os.listdir(structures_path) if structure.__contains__('ctv') or structure.__contains__('CTV')]

        high_dose_CTV = np.zeros(shape = CT_shape, dtype=np.bool8)
        medium_dose_CTV = np.zeros(shape = CT_shape, dtype=np.bool8)
        low_dose_CTV = np.zeros(shape = CT_shape, dtype=np.bool8)
        
        for CTV_structure in CTV_structures:
            
            if CTV_structure.__contains__('70') or CTV_structure.__contains__('7000') or CTV_structure.__contains__('65') or CTV_structure.__contains__('6500'):
                
                CTV_structure_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/' + str(CTV_structure)

                structure_obj = nib.load(CTV_structure_path)
                structure_hdr = structure_obj.header
                structure_affine = structure_obj.affine
                structure_img = np.array(structure_obj.get_fdata(), dtype = np.bool8)
                
                high_dose_CTV = np.logical_or(structure_img, high_dose_CTV)
                
            if CTV_structure.__contains__('60') or CTV_structure.__contains__('6000'):
    
                CTV_structure_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/' + str(CTV_structure)

                structure_obj = nib.load(CTV_structure_path)
                structure_hdr = structure_obj.header
                structure_affine = structure_obj.affine
                structure_img = np.array(structure_obj.get_fdata(), dtype = np.bool8)

                medium_dose_CTV = np.logical_or(structure_img, medium_dose_CTV)

            if CTV_structure.__contains__('56') or CTV_structure.__contains__('5600') or CTV_structure.__contains__('54') or CTV_structure.__contains__('5400'):
    
                CTV_structure_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/' + str(CTV_structure)

                structure_obj = nib.load(CTV_structure_path)
                structure_hdr = structure_obj.header
                structure_affine = structure_obj.affine
                structure_img = np.array(structure_obj.get_fdata(), dtype = np.bool8)

                low_dose_CTV = np.logical_or(structure_img, low_dose_CTV)

        if np.sum(high_dose_CTV) != 0:
                
            high_dose_CTV_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/BIN_CTV_HIGH.nii.gz'
            NewNiftiObj = nib.Nifti1Image(high_dose_CTV, structure_affine, structure_hdr)
            nib.save(NewNiftiObj, high_dose_CTV_path)

        if np.sum(medium_dose_CTV) != 0:
                
            medium_dose_CTV_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/BIN_CTV_MEDIUM.nii.gz'
            NewNiftiObj = nib.Nifti1Image(medium_dose_CTV, structure_affine, structure_hdr)
            nib.save(NewNiftiObj, medium_dose_CTV_path)

        if np.sum(low_dose_CTV) != 0:
            
            low_dose_CTV_path = patient_path + '/' + str(CT_folder) + '/STRUCTURES/BIN_CTV_LOW.nii.gz'
            NewNiftiObj = nib.Nifti1Image(low_dose_CTV, structure_affine, structure_hdr)
            nib.save(NewNiftiObj, low_dose_CTV_path)

        
        
        structures = os.listdir(patient_path + '/' + str(CT_folder) + '/STRUCTURES/')

        for structure in structures:

            if structure not in ['BIN_BRAINSTEM.nii.gz', 'BIN_CORD.nii.gz', 'BIN_CTV_HIGH.nii.gz', 'BIN_CTV_LOW.nii.gz', 'BIN_CTV_MEDIUM.nii.gz', 'BIN_PAROTIDL.nii.gz', 'BIN_PAROTIDR.nii.gz', 'BIN_BODY.nii.gz' \
                                 'BIN_BRAINSTEM_OG.nii.gz', 'BIN_CORD_OG.nii.gz', 'BIN_CTV_HIGH_OG.nii.gz', 'BIN_CTV_LOW_OG.nii.gz', 'BIN_CTV_MEDIUM_OG.nii.gz', 'BIN_PAROTIDL_OG.nii.gz', 'BIN_PAROTIDR_OG.nii.gz', 'BIN_BODY_OG.nii.gz']:

                os.remove(patient_path + '/' + str(CT_folder) + '/STRUCTURES/'+ str(structure))
        
        