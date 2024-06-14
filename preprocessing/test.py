import os 
import shutil 
import nibabel as nib 
import numpy as np 


def flip_img(img_path, new_img_path = False):
    
    img_data, img_affine, img_header = get_img_objects(img_path)
        
    img_flipped = np.flip(img_data, axis = 0)
    NewNiftiObj = nib.Nifti1Image(img_flipped, img_affine, img_header)
    if new_img_path != False:
        os.rename(img_path, img_path[:-7] + '_OG.nii.gz')
        nib.save(NewNiftiObj, img_path)
    else:
        nib.save(NewNiftiObj, img_path)
    
    
def rename_parotid(file_path, name, flip = False):
    new_name = name
    
    if name.__contains__('PAROTIDL'):
        if flip:
            new_name=new_name.replace('PAROTIDL','RPAROTID')
        else:
            new_name=new_name.replace('PAROTIDL','LPAROTID')
        source = file_path + '/' + name
        destination = file_path + '/' + new_name
        os.rename(source, destination)   
    elif name.__contains__('PAROTIDR'):
        if flip:
            new_name=new_name.replace('PAROTIDR','LPAROTID')
        else:
            new_name=new_name.replace('PAROTIDR','RPAROTID')
        source = file_path + '/' + name
        destination = file_path + '/' + new_name
        os.rename(source, destination)

def get_img_objects(img_path, data_type = np.float32):
    
    img_obj = nib.load(img_path)
    img_data = np.array(img_obj.get_fdata(), dtype = data_type)
    img_affine = img_obj.affine
    img_header = img_obj.header

    return img_data, img_affine, img_header


path = 'D:/CUH_HN/NIFTI/'

patients = [folder for folder in os.listdir(path) if folder.startswith('CUH-UCL')]

for patient in patients:

    CT_folders = [folder for folder in os.listdir(path + str(patient)) if folder.startswith('CT_')]
    
    for CT_folder in CT_folders:
        
        structures = [structure for structure in os.listdir('D:/CUH_HN/NIFTI/'+str(patient)+'/'+str(CT_folder) + '/STRUCTURES/')]
        
        for structure in structures:

            CT_path = 'D:/CUH_HN/NIFTI/'+str(patient)+'/'+str(CT_folder) + '/CT_OG.nii.gz'

            if os.path.exists(CT_path):
                

                structure_path = 'D:/CUH_HN/NIFTI/'+str(patient)+'/'+str(CT_folder) + '/STRUCTURES/' + str(structure)
                flip_img(structure_path, new_img_path=True)
                rename_parotid('D:/CUH_HN/NIFTI/'+str(patient)+'/'+str(CT_folder) + '/STRUCTURES/', structure, flip = True)

            else:
                rename_parotid('D:/CUH_HN/NIFTI/'+str(patient)+'/'+str(CT_folder) + '/STRUCTURES/', structure, flip = False)
            