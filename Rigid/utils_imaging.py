import numpy as np 
import nibabel as nib 

def wrapper_function(niftiImgIn, minInt=-1000, maxInt=1000, newMinInt=0, newMaxInt=1):
    """
    Wraps around the two image scaling functions
    """
    niftiImgOut = clip_intensities(niftiImgIn, minInt=-1000, maxInt=1000)

    RescaledNiftiImgOut = rescale_intensities(niftiImgOut, newMinInt=0, newMaxInt=1)

    return RescaledNiftiImgOut

def clip_intensities(niftiImgIn, minInt=-1000, maxInt=1000):
    """
    Clip intensities of nifti image to set range
    """
    img = niftiImgIn.get_fdata()
    affine = niftiImgIn.affine
    hdr = niftiImgIn.header

    img[np.isnan(img)] = 700

    img[img <= minInt] = minInt
    img[img >= maxInt] = maxInt
    niftiImgOut = nib.Nifti1Image(img, affine, hdr)
    return niftiImgOut


def rescale_intensities(niftiImgIn, newMinInt=0, newMaxInt=1):
    """
    Rescale intensities of nifti image between set boundaries
    """
    img = niftiImgIn.get_fdata()
    affine = niftiImgIn.affine
    hdr = niftiImgIn.header

    minIn = np.amin(img)
    maxIn = np.amax(img)

    rescaledImg = (((img - minIn)/(maxIn - minIn)) * (newMaxInt - newMinInt)) + newMinInt
    rescaledImg = rescaledImg * 255
    niftiImgOut = nib.Nifti1Image(rescaledImg, affine, hdr)
    return niftiImgOut

def get_pixel_spacing(CT_img):
    
    return CT_img.header['pixdim'][1:3]

def get_slice_thickness(CT_img):

    return CT_img.header['pixdim'][3]

