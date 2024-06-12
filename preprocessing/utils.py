from DicomRTTool.ReaderWriter import DicomReaderWriter, ROIAssociationClass
import SimpleITK as sitk 
import os
import pandas as pd
import numpy as np 
import nibabel as nib 
import shutil 

def convert_OARs(dicom_Reader, contour_names, associations, structures_path):
    

    for i, contour_name in enumerate(contour_names):

        #print(contour_name)
        #print(associations[i])
    
    
        dicom_Reader.set_contour_names_and_associations(contour_names=[contour_name], associations=[ROIAssociationClass(contour_name, associations[i])])
        #path = dicom_Reader.where_is_ROI()
        #print('-----')
        #print(path)
        
        dicom_Reader.get_mask()

        #load mask
        mask_sitk_handle = dicom_Reader.annotation_handle

        
        print(contour_name)
        # saves as boolean mask 
        roi_name =  structures_path + 'BIN_'+ contour_name + '.nii.gz'
        sitk.WriteImage(mask_sitk_handle, roi_name)
        

    # copy file name of RTSTRUCT into the DICOM folder 
        test_binary_obj = nib.load(roi_name)
        if np.sum(test_binary_obj.get_fdata()) ==0:
            os.remove(roi_name)
    

def convert_CTVs(dicom_Reader, structures_path):

    ROIS = dicom_Reader.return_rois()
    print(ROIS)


    CTV_ROIS = [ROI for ROI in ROIS if ROI.__contains__('CTV') or ROI.__contains__('ctv')]

    for contour_name in CTV_ROIS:

        if contour_name.__contains__('/'):

            contour_name = contour_name.replace('/', '_')

        print('Contour name: ..... ')
        print(contour_name)
        #print(associations[i])
    
    
        dicom_Reader.set_contour_names_and_associations(contour_names=[contour_name], associations=[ROIAssociationClass(contour_name, contour_name)])
        #path = dicom_Reader.where_is_ROI()
        #print('-----')
        #print(path)
        
        dicom_Reader.get_mask()

        #load mask
        mask_sitk_handle = dicom_Reader.annotation_handle

        
        print(contour_name)
        # saves as boolean mask 
        roi_name =  structures_path + 'BIN_'+ contour_name + '.nii.gz'
        sitk.WriteImage(mask_sitk_handle, roi_name)
        

    # copy file name of RTSTRUCT into the DICOM folder 
        test_binary_obj = nib.load(roi_name)
        if np.sum(test_binary_obj.get_fdata()) ==0:
            os.remove(roi_name)
    


def DICOM_CONVERT(Base_path, contour_names, associations, results_path, patient):

    '''
    function which converts both images and structures at the same time 

    Base_path: string. folder in which the HN_x folders are stored.
    contour_names: list of strings. List of names of contours which you want to save
    associations: list of lists of strings. List of possible associated contour names as contour_names.
    ''' 

    #The number of nested lists in associations should be equal to the length of contour_names
    if len(associations) != len(contour_names):
        raise Exception('Each contour should have a list of associations')


    # for catching patients which didnt work
    RTSTRUCT_conversion_log = open("RTSTRUCT_conversion_log.txt", mode="a")
  
    # indefies the number of unique series UID images within the patient specific folder
    dicom_Reader = DicomReaderWriter()
    dicom_Reader.walk_through_folders(Base_path)
    indexes = dicom_Reader.images_dictionary

    # loops through the number of unique series UIDs
    for img_index in range(0, len(indexes)):

        # gets image object 
        dicom_Reader.set_index(img_index)

        # creates folders within the patient folder for saving images
        # see above for exact folder architecture which is stored 
        nifti_folder = results_path + '/' + str(patient) + '/'
        if not os.path.exists(nifti_folder):
            os.mkdir(nifti_folder)

        # this returns the series date of the scan 
        series_date = str(dicom_Reader.return_key_info('0008|0021'))
        machine_imaging_type = str(dicom_Reader.return_key_info('0008|0070'))
        
        # determines, based on the scanner, whether image is CT or CBCT
        # TOSHIBA CT scanner 
        # 'Varian Medical Systems' is the CBCT 
        if machine_imaging_type.__contains__('TOSHIBA'):
            filename_prefix = 'CT'
        else:
            filename_prefix = 'CBCT'
        
        image_folder = nifti_folder +  '/' + str(filename_prefix) + '_' + str(series_date) + '/'
        if not os.path.exists(image_folder):
            os.mkdir(image_folder)

        
        # gives name of nifti image
        # if ile already exists in that folder, it just creates a second version
        # to catch maybe double imaging on the same day 
        filename =  image_folder + str(filename_prefix) + '.nii.gz'
        if os.path.exists(filename):
            filename = image_folder + str(filename_prefix) + '_2.nii.gz'

        '''    
        if filename_prefix == 'CBCT':
            dicom_Reader.get_images()
            dicom_sitk_handle = dicom_Reader.dicom_handle
            sitk.WriteImage(dicom_sitk_handle, filename)  
        
        '''
        
        if filename_prefix == 'CT':

            # create a directory for the structures
            structures_path = image_folder + '/STRUCTURES/'
            if not os.path.exists(structures_path):
                os.mkdir(structures_path)

            
            #convert_OARs(dicom_Reader, contour_names, associations, structures_path)
            convert_CTVs(dicom_Reader, structures_path)


        
    # if no conversions failed you note this down too 
    if os.path.getsize("RTSTRUCT_conversion_log.txt") == 0:
        RTSTRUCT_conversion_log.write('No failed conversions.')
            
            




def get_contour_names_and_associations(path):
    '''
    path: string. path to excel file which contains the contour names as a header
    and all the listed associations in columns 

    returns: contour names and associations in list formats
    '''
   
    contour_names_excel = pd.read_excel(path, sheet_name = 'H&N', header =0 , dtype = str)
    # these are the contour names you want
    contour_names = list(contour_names_excel.columns.values)
    # possible list of associated contour names
    associations = []
    filtered_contour_names = []
    for roi in contour_names:

        if roi in ['BRAINSTEM', 'CORD', 'PAROTIDL', 'PAROTIDR', 'BODY']:
            #print(r)
            associations.append(list(contour_names_excel.loc[:,roi].dropna()))
            filtered_contour_names.append(roi)
        


    return filtered_contour_names, associations

            


def get_image_objects(path):

    nifti_obj = nib.load(path)
    nifti_img = nifti_obj.get_fdata()
    nifti_affine = nifti_obj.affine
    nifti_header = nifti_obj.header

    del nifti_obj

    return nifti_img, nifti_affine, nifti_header


def check_DICOM_conversions(results_path):

        # for catching patients which didnt work
        RTSTRUCT_conversion_log = open("RTSTRUCT_conversion_log.txt", mode="a")

        CBCT_folders = [file for file in os.listdir(results_path) if file.startswith('CBCT')]

        for CBCT_folder in CBCT_folders:

            CBCT_path = os.path.join(results_path, CBCT_folder, 'CBCT.nii.gz')

            nifti_img, nifti_affine, nifti_header = get_image_objects(CBCT_path)

            RTSTRUCT_conversion_log.write(str(CBCT_folder)+ str(nifti_header['dim'][1:4]) + '\n')


def crop_Img(img_path, z_slice, new_img_path):
    
    reg_transform = 'C:/Users/poppy/Documents/Nifty/niftyreg_install/bin/reg_transform.exe'
    
    affine_matrix_path = 'Sform_matrix_path.txt'

    calc_Sform_matrix(img_path, affine_matrix_path, z_slice)

    img_data, img_affine, img_header = get_img_objects(img_path)
    

    img_data_copy = img_data[:,:,z_slice:]

    # save the cropped image 
    NewNiftiObj = nib.Nifti1Image(img_data_copy, img_affine, img_header)
    nib.save(NewNiftiObj, new_img_path) 

    UpdSform(reg_transform, new_img_path, affine_matrix_path, new_img_path)

    os.remove(affine_matrix_path)
    
def UpdSform(reg_transform_path, img_to_be_updated_path, affine_matrix_path, updated_img_path):
    
    command = reg_transform_path +' -updSform ' + img_to_be_updated_path + ' ' + affine_matrix_path + ' ' + updated_img_path
    os.system(command)

def calc_Sform_matrix(img_path, affine_matrix_path, z_slice):
        
        # calculates the z shift to move the cropped image to
        _, _, img_header = get_img_objects(img_path)
        del _
        slice_width = img_header['pixdim'][3]
        z_shift = z_slice * slice_width

        # creates the transformation matrix to use in updating the Sform 
        identity_matrix = np.identity(4)
        identity_matrix[2][3] = -z_shift 
        #print(identity_matrix)

        # saves the transformation matrix to use later on with updating the Sform 
        identity_matrix = pd.DataFrame(identity_matrix)
        np.savetxt(affine_matrix_path, identity_matrix, fmt='%d')


def get_img_objects(img_path):
    
    img_obj = nib.load(img_path)
    img_data = np.array(img_obj.get_fdata())
    img_affine = img_obj.affine
    img_header = img_obj.header

    return img_data, img_affine, img_header
