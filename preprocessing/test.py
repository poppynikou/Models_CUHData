import os 
import shutil 


path = 'D:/CUH_HN/NIFTI/'

patients = [folder for folder in os.listdir(path) if folder.startswith('CUH-UCL')]

for patient in patients:

    CT_folders = [folder for folder in os.listdir(path + str(patient)) if folder.startswith('CT_')]

    for CT_folder in CT_folders:

        structures_path = 'D:/CUH_HN/NIFTI/'+str(patient)+'/'+str(CT_folder) + '/STRUCTURES'

        shutil.rmtree(structures_path)