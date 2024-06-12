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
#ignore_patients = [1,2,4,5,6,7,8,9,10,11,12,13,13,25,26,29,30,35,36,38,40,43,46,45,44,42,39,33,32,31,28,27,24,23,22]
#ignore = ['CUH-UCL-02-0' + str("%.2d" % i) for i in ignore_patients]

#patients = [folder for folder in os.listdir('D:/CUH_HN/NIFTI') if folder.__contains__('CUH')]

#patients = natsorted(set([file[0:14] for file in os.listdir('D:/CUH_HN/2D_slices')]))
#BCT_list = [[file[15:28] for file in os.listdir('D:/CUH_HN/2D_slices') if file[0:14] == patient] for patient in patients]

patients = ['CUH-UCL-02-011', 'CUH-UCL-02-018', 'CUH-UCL-02-029']

for index, patient in enumerate(patients):
    
    Patient_No = str(patient)

    patient_path = os.path.join(path, Patient_No)

    images = [file for file in os.listdir(patient_path) if file not in ['CBCT_GROUPWISE', 'CBCT_pCT']]

    for image_ in images:

            img_path = os.path.join(patient_path,image_)

            niftis = os.listdir(img_path)
            
            for nifti in niftis:
                
                if  nifti == 'MASKED_CBCT.nii.gz':# or nifti == 'MASKED_CBCT.nii.gz': 
                    try:
                        CT_path = os.path.join(img_path, nifti)
                    
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


                        file_name = 'D:/CUH_HN/2D_slices/' + str(patient) + '_' + str(image_) + '_' + str(nifti) + '.png'

                        #fig.suptitle(file_name, fontsize=8)
                        #plt.subplot_tool()
                        plt.tight_layout()
                        plt.savefig(file_name)
                        plt.close()
                    except:
                        print(CT_path)
                        
                        
                    
