import os 
import matplotlib.pyplot as plt 
import nibabel as nib 
from utils_imaging import * 
from math import isnan
from natsort import natsorted

def get_img_objects(img_path, data_type = np.float32):
    
        img_obj = nib.load(img_path)
        img_data = np.array(img_obj.get_fdata(), dtype = data_type)
        img_affine = img_obj.affine
        img_header = img_obj.header

        return img_data, img_affine, img_header

path = 'D:/CUH_HN/NIFTI/'

patients = [folder for folder in os.listdir('D:/CUH_HN/NIFTI') if folder.__contains__('CUH')]
CBCTS = [[file for file in os.listdir(path + str(patient) + '/CBCT_GROUPWISE/postprocessing')] for patient in patients]


for index, patient in enumerate(patients):
    print(patient)
    Patient_No = str(patient)

    patient_path = os.path.join(path, Patient_No)
    print(patient_path)

    for image_ in CBCTS[index]:
       
        CT_path = os.path.join(patient_path, 'CBCT_GROUPWISE/postprocessing/', image_)
        print(CT_path)
    
        fig = plt.figure(figsize = (10,8))
        
        # (nrows, ncolumns, index)
        ax1 = plt.subplot2grid(shape = (2,2), loc = (0,0), rowspan = 2, colspan =1)
        ax2 = plt.subplot2grid(shape = (2,2), loc = (0,1))
        ax3 = plt.subplot2grid(shape = (2,2), loc = (1,1))
        
        CT_img = nib.load(CT_path)

        image = wrapper_function(CT_img, minInt=-1000, maxInt=1000, newMinInt=0, newMaxInt=1)

        image = np.transpose(np.asarray(image.get_fdata(),  dtype=np.float32))
        img_shape = np.shape(image)

        Axial_slice = int(img_shape[0]/2)
        Coronal_slice = int(img_shape[1]/2)
        Sagital_slice = int(img_shape[2]/2)

        Pixel_spacing = get_pixel_spacing(CT_img)
        Slice_thickness = get_slice_thickness(CT_img)
        ax_aspect = Pixel_spacing[1]/Pixel_spacing[0]
        sag_aspect = Slice_thickness/Pixel_spacing[1]
        cor_aspect = Slice_thickness/Pixel_spacing[0]


        ax1.imshow(np.flipud(image[Axial_slice, :,:].astype('float32')), cmap = 'gray', origin = 'lower')
        ax2.imshow(image[:,Coronal_slice,:].astype('float32'), cmap='gray', origin = 'lower')
        ax3.imshow(image[:,:,Sagital_slice].astype('float32'), cmap='gray', origin = 'lower') 

        
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)
        ax2.get_xaxis().set_visible(False)
        ax2.get_yaxis().set_visible(False)
        ax3.get_xaxis().set_visible(False)
        ax3.get_yaxis().set_visible(False)


        ax1.set_aspect(ax_aspect)
        ax2.set_aspect(cor_aspect)
        ax3.set_aspect(sag_aspect)


        file_name = 'D:/CUH_HN/2D_slices/' + str(patient) + '_' + str(image_)  + '.png'

        #fig.suptitle(file_name, fontsize=8)
        #plt.subplot_tool()
        plt.tight_layout()
        plt.savefig(file_name)
        plt.close()

        
