import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np 
import os
import nibabel as nib
from utils_visualisations import * 


data_path = os.getcwd() + '\Animations'


#not_files = ['CHRISTIE_Average', 'CHRISTIE_PCA1_PCA2_plus1std']
#patients = ['CHRISTIE_PCA1_PCA2_plus1std']#[file for file in os.listdir(data_path) if file not in not_files]

Axial_slice = 65
Coronal_slice = 234
Sagital_slice = 264
lower_y = 100
upper_y = 330
lower_x = 150
upper_x = 375
lower_z = 20

Imgs = np.arange(0,39)
#[file for file in os.listdir(data_path + str(Patient_no)) if file.startswith('Img')]



for index, Img in enumerate(Imgs):    

    fig = plt.figure(figsize = (5,4))
    
    # (nrows, ncolumns, index)
    ax1 = plt.subplot2grid(shape = (2,2), loc = (0,0), rowspan = 2, colspan =1)
    ax2 = plt.subplot2grid(shape = (2,2), loc = (0,1))
    ax3 = plt.subplot2grid(shape = (2,2), loc = (1,1))
    

    CT_path = 'D:/CUH_HN/NIFTI/PSM_noLOO_T_All_CPS_4_PatientSpace/CUH-UCL-02-005/Img_'+str(Img)+'.nii.gz'
    print(CT_path)

    CT_img = nib.load(CT_path)
    #CT_img_day_0 = nib.load(CT_day_0_path)
    minInt = np.amin(CT_img.get_fdata())
    maxInt = np.amax(CT_img.get_fdata())

    image = wrapper_function(CT_img, minInt=minInt, maxInt=maxInt, newMinInt=0, newMaxInt=255)

    #heatmap = np.asarray(image_0.get_fdata(),  dtype=np.float32) - np.asarray(image.get_fdata(),  dtype=np.float32)

    image = np.transpose(np.asarray(image.get_fdata(),  dtype=np.float32))
    #heatmap = np.transpose(np.asarray(heatmap,  dtype=np.float32))

    #image_0 = wrapper_function(CT_img_day_0, minInt=-1000, maxInt=1000, newMinInt=0, newMaxInt=1)
    #image_0 = np.transpose(np.asarray(image_0.get_fdata(),  dtype=np.float32))

    Pixel_spacing = get_pixel_spacing(CT_img)
    Slice_thickness = get_slice_thickness(CT_img)
    ax_aspect = Pixel_spacing[1]/Pixel_spacing[0]
    sag_aspect = Slice_thickness/Pixel_spacing[1]
    cor_aspect = Slice_thickness/Pixel_spacing[0]

    #difference_img = image - image_0

    #print(np.amax(difference_img), np.amin(difference_img))

    ax1.imshow(np.flipud(image[Axial_slice, lower_y:upper_y, lower_x:upper_x].astype('uint8')), cmap = 'gray', origin = 'lower')
    ax2.imshow(image[lower_z:,Coronal_slice,10:502].astype('uint8'), cmap='gray', origin = 'lower')
    ax3.imshow(image[lower_z:,10:502,Sagital_slice].astype('uint8'), cmap='gray', origin = 'lower') 

    #ax1.imshow(np.flipud(heatmap[Axial_slice, lower_y:upper_y, lower_x:upper_x].astype('uint8')), cmap = 'Wistia', origin = 'lower', alpha = 0.5)
    #ax2.imshow(heatmap[lower_z:,Coronal_slice,10:502].astype('uint8'), cmap='Wistia', origin = 'lower', alpha = 0.5)
    #ax3.imshow(heatmap[lower_z:,10:502,Sagital_slice].astype('uint8'), cmap='Wistia', origin = 'lower', alpha = 0.5) 

    #ax1.imshow(np.flipud(difference_img[Axial_slice, lower_y:upper_y, lower_x:upper_x].astype('uint8')), origin = 'lower', norm = colors.Normalize(vmin=-300, vmax = 300), cmap = 'PiYG')
    #ax2.imshow(difference_img[lower_z:,Coronal_slice,10:502].astype('uint8'),  origin = 'lower', norm = colors.Normalize(vmin=-300, vmax = 300), cmap = 'PiYG')
    #ax3.imshow(difference_img[lower_z:,10:502,Sagital_slice].astype('uint8'),  origin = 'lower', norm = colors.Normalize(vmin=-300, vmax = 300), cmap = 'PiYG') 


    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)


    ax1.set_aspect(ax_aspect)
    ax2.set_aspect(cor_aspect)
    ax3.set_aspect(sag_aspect)

    fig.suptitle('Day ' + str(index+1), fontsize=10)
    #plt.subplot_tool()
    plt.tight_layout()
    plt.savefig(data_path + str(index) + ".jpg")
    plt.close()
    

            
#pillow to save all frames as an animation in a gif file
from PIL import Image

images = [Image.open(data_path + str(timepoint) + ".jpg") for timepoint in range(0, 35,1)]

figname = data_path + '/CUH-UCL-02-005.gif'
#data_path + '/' + str(Patient_no) +'.gif'
images[0].save(figname, save_all=True, append_images=images[1:], duration=200, loop=0)


test = os.listdir(data_path)

for item in test:
    if item.endswith(".jpg"):
            os.remove(os.path.join(data_path, item))


    