import os 
from datetime import date
import numpy 
from scipy.ndimage import _ni_support
from scipy.ndimage.morphology import distance_transform_edt, binary_erosion, generate_binary_structure
import nibabel as nib 
import numpy as np 


def get_time_points(CBCT_dates):
    '''
    This assumes that the first time points is the first CBCT in the series
    '''
    # find the date for the first CBCT date 
    CBCT_dates = [date(int(CBCT_date[0:4]), int(CBCT_date[4:6]),  int(CBCT_date[6:8])) for CBCT_date in CBCT_dates]
    day_zero = min(CBCT_dates)

    time_points = [(timepoint - day_zero).days for timepoint in CBCT_dates]

    return time_points

def inv_Dff(ref_img, input_transformation, float_img, inverse_transform):
    
    command = 'reg_transform -ref ' + str(ref_img) + ' -invNrr ' + str(input_transformation) + ' ' + str(float_img) + ' ' + str(inverse_transform)
    os.system(command)

def check_list_contained(A,B):
        # check if list A is contained in list B 
        A_arr = np.array(A)
        B_arr = np.array(B)

        for i in range(len(B_arr)):

            if np.array_equal(A_arr, B_arr[i:i+len(A_arr)]):

                return True

        return False

def save_as_int(img_path):
    
    img_obj, img_affine, img_header = get_image_objects(img_path)
    img_obj_copy = np.array(img_obj, dtype = np.bool8)

    NewNiftiObj = nib.Nifti1Image(img_obj_copy, img_affine, img_header)
    NewNiftiObj.set_data_dtype(np.int8)
    nib.save(NewNiftiObj, img_path)

def remove_double_time_points(time_points):

    # If in a set of time points there are daily imaging one week
    # remove the tuesday and thursday 
    for i in np.arange(0,time_points[-1]-5):

        count_list = np.array([0,1,2,3,4]) + int(i)

        if check_list_contained(count_list, time_points):

            time_points.remove(count_list[1])
            time_points.remove(count_list[3])

    return time_points

def velocity_to_deformationfield(ref_img, input_transformation, output_transformation):
    
    print(ref_img)
    print(input_transformation)
    print(output_transformation)
    command = 'reg_transform -ref ' + ref_img + ' -def ' + input_transformation + ' ' + output_transformation
    os.system(command)

    return 

def compose_transform(ref_img, composed_transformations, transformation2, transformation1):


    '''
    Compose two transformations of any recognised type and returns a deformation field.
    composed_transformations(x) = transformation2(transformation1(x)).
    '''
    
    command = 'reg_transform -ref ' + ref_img + ' -comp ' + transformation1 + ' ' + transformation2 + ' ' + composed_transformations

    os.system(command)


def Transport_Model_to_Patient(model, input_path, output_path, output_transformation, patient_to_modelspace_data, CBCT_img_modelspace, CBCT_img_patientspace):
    '''
    input_path: path to velocity field in model space
    output_path: path to deformation field in patient space
    patient_to_modelspace_data: folder in which transformations for deforming patient to model space data is stored (i.e. Transported patients)
    '''

    # convert velocity field to deformation field
    #velocity_to_deformationfield(CBCT_img_modelspace, input_path, output_transformation)

    T2 = patient_to_modelspace_data + '/T_model.nii.gz'
    T2_inv = patient_to_modelspace_data + '/inv_T_model.nii.gz'
    
    intermediate_composed_transformations = patient_to_modelspace_data + '/' + str(model) + '_intermediate_composed_transformations.nii.gz'
    compose_transform(CBCT_img_modelspace, intermediate_composed_transformations, output_transformation, T2_inv)

    compose_transform(CBCT_img_patientspace, output_path, T2, intermediate_composed_transformations)

    #os.remove(intermediate_composed_transformations)
    #os.remove(output_transformation)

def resampleBINImg(ref_img, float_img, transformation, resampled_img):
    
    command = 'reg_resample -ref ' + ref_img +  ' -flo ' + float_img + ' -trans ' + transformation + ' -res ' + resampled_img + ' -inter 0 -omp 12 -pad 0'
    os.system(command)

def resampleImg(ref_img, float_img, transformation, resampled_img):
    
    command = 'reg_resample -ref ' + ref_img +  ' -flo ' + float_img + ' -trans ' + transformation + ' -res ' + resampled_img + ' -inter 3 -omp 12 -pad -1000'
    os.system(command)

#medpy functions 

def __surface_distances(result, reference, voxelspacing=None, connectivity=1):
    """
    The distances between the surface voxel of binary objects in result and their
    nearest partner surface voxel of a binary object in reference.
    """
    result = numpy.atleast_1d(result.astype(numpy.bool))
    reference = numpy.atleast_1d(reference.astype(numpy.bool))
    if voxelspacing is not None:
        voxelspacing = _ni_support._normalize_sequence(voxelspacing, result.ndim)
        voxelspacing = numpy.asarray(voxelspacing, dtype=numpy.float64)
        if not voxelspacing.flags.contiguous:
            voxelspacing = voxelspacing.copy()
            
    # binary structure
    footprint = generate_binary_structure(result.ndim, connectivity)
    
    # test for emptiness
    if 0 == numpy.count_nonzero(result): 
        raise RuntimeError('The first supplied array does not contain any binary object.')
    if 0 == numpy.count_nonzero(reference): 
        raise RuntimeError('The second supplied array does not contain any binary object.')    
            
    # extract only 1-pixel border line of objects
    result_border = result ^ binary_erosion(result, structure=footprint, iterations=1)
    reference_border = reference ^ binary_erosion(reference, structure=footprint, iterations=1)
    
    # compute average surface distance        
    # Note: scipys distance transform is calculated only inside the borders of the
    #       foreground objects, therefore the input has to be reversed
    dt = distance_transform_edt(~reference_border, sampling=voxelspacing)
    sds = dt[result_border]
    
    return sds

def hd(result, reference, voxelspacing=None, connectivity=1):
    """
    Hausdorff Distance.
    
    Computes the (symmetric) Hausdorff Distance (HD) between the binary objects in two
    images. It is defined as the maximum surface distance between the objects.
    
    Parameters
    ----------
    result : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    reference : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    voxelspacing : float or sequence of floats, optional
        The voxelspacing in a distance unit i.e. spacing of elements
        along each dimension. If a sequence, must be of length equal to
        the input rank; if a single number, this is used for all axes. If
        not specified, a grid spacing of unity is implied.
    connectivity : int
        The neighbourhood/connectivity considered when determining the surface
        of the binary objects. This value is passed to
        `scipy.ndimage.morphology.generate_binary_structure` and should usually be :math:`> 1`.
        Note that the connectivity influences the result in the case of the Hausdorff distance.
        
    Returns
    -------
    hd : float
        The symmetric Hausdorff Distance between the object(s) in ```result``` and the
        object(s) in ```reference```. The distance unit is the same as for the spacing of 
        elements along each dimension, which is usually given in mm.
        
    See also
    --------
    :func:`assd`
    :func:`asd`
    
    Notes
    -----
    This is a real metric. The binary images can therefore be supplied in any order.
    """
    hd1 = __surface_distances(result, reference, voxelspacing, connectivity).max()
    hd2 = __surface_distances(reference, result, voxelspacing, connectivity).max()
    hd = max(hd1, hd2)
    return hd

def dc(result, reference):
    r"""
    Dice coefficient
    
    Computes the Dice coefficient (also known as Sorensen index) between the binary
    objects in two images.
    
    The metric is defined as
    
    .. math::
        
        DC=/frac{2|A/cap B|}{|A|+|B|}
        
    , where :math:`A` is the first and :math:`B` the second set of samples (here: binary objects).
    
    Parameters
    ----------
    result : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    reference : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    
    Returns
    -------
    dc : float
        The Dice coefficient between the object(s) in ```result``` and the
        object(s) in ```reference```. It ranges from 0 (no overlap) to 1 (perfect overlap).
        
    Notes
    -----
    This is a real metric. The binary images can therefore be supplied in any order.
    """
    result = numpy.atleast_1d(result.astype(numpy.bool))
    reference = numpy.atleast_1d(reference.astype(numpy.bool))
    
    intersection = numpy.count_nonzero(result & reference)
    
    size_i1 = numpy.count_nonzero(result)
    size_i2 = numpy.count_nonzero(reference)
    
    try:
        dc = 2. * intersection / float(size_i1 + size_i2)
    except ZeroDivisionError:
        dc = 0.0
    
    return dc

def asd(result, reference, voxelspacing=None, connectivity=1):
    """
    Average surface distance metric.
    
    Computes the average surface distance (ASD) between the binary objects in two images.
    
    Parameters
    ----------
    result : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    reference : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    voxelspacing : float or sequence of floats, optional
        The voxelspacing in a distance unit i.e. spacing of elements
        along each dimension. If a sequence, must be of length equal to
        the input rank; if a single number, this is used for all axes. If
        not specified, a grid spacing of unity is implied.
    connectivity : int
        The neighbourhood/connectivity considered when determining the surface
        of the binary objects. This value is passed to
        `scipy.ndimage.morphology.generate_binary_structure` and should usually be :math:`> 1`.
        The decision on the connectivity is important, as it can influence the results
        strongly. If in doubt, leave it as it is.
    
    Returns
    -------
    asd : float
        The average surface distance between the object(s) in ``result`` and the
        object(s) in ``reference``. The distance unit is the same as for the spacing
        of elements along each dimension, which is usually given in mm.
        
    See also
    --------
    :func:`assd`
    :func:`hd`
    
    
    Notes
    -----
    This is not a real metric, as it is directed. See `assd` for a real metric of this.
    
    The method is implemented making use of distance images and simple binary morphology
    to achieve high computational speed.
    
    Examples
    --------
    The `connectivity` determines what pixels/voxels are considered the surface of a
    binary object. Take the following binary image showing a cross
    
    >>> from scipy.ndimage.morphology import generate_binary_structure
    >>> cross = generate_binary_structure(2, 1)
    array([[0, 1, 0],
           [1, 1, 1],
           [0, 1, 0]])
           
    With `connectivity` set to `1` a 4-neighbourhood is considered when determining the
    object surface, resulting in the surface
    
    .. code-block:: python
    
        array([[0, 1, 0],
               [1, 0, 1],
               [0, 1, 0]])
           
    Changing `connectivity` to `2`, a 8-neighbourhood is considered and we get:
    
    .. code-block:: python
    
        array([[0, 1, 0],
               [1, 1, 1],
               [0, 1, 0]])
           
    , as a diagonal connection does no longer qualifies as valid object surface.
    
    This influences the  results `asd` returns. Imagine we want to compute the surface
    distance of our cross to a cube-like object:
    
    >>> cube = generate_binary_structure(2, 1)
    array([[1, 1, 1],
           [1, 1, 1],
           [1, 1, 1]])
           
    , which surface is, independent of the `connectivity` value set, always
    
    .. code-block:: python
    
        array([[1, 1, 1],
               [1, 0, 1],
               [1, 1, 1]])
           
    Using a `connectivity` of `1` we get
    
    >>> asd(cross, cube, connectivity=1)
    0.0
    
    while a value of `2` returns us
    
    >>> asd(cross, cube, connectivity=2)
    0.20000000000000001
    
    due to the center of the cross being considered surface as well.
    
    """
    sds = __surface_distances(result, reference, voxelspacing, connectivity)
    asd = sds.mean()
    return asd

def hd95(result, reference, voxelspacing=None, connectivity=1):
    """
    95th percentile of the Hausdorff Distance.

    Computes the 95th percentile of the (symmetric) Hausdorff Distance (HD) between the binary objects in two
    images. Compared to the Hausdorff Distance, this metric is slightly more stable to small outliers and is
    commonly used in Biomedical Segmentation challenges.

    Parameters
    ----------
    result : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    reference : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    voxelspacing : float or sequence of floats, optional
        The voxelspacing in a distance unit i.e. spacing of elements
        along each dimension. If a sequence, must be of length equal to
        the input rank; if a single number, this is used for all axes. If
        not specified, a grid spacing of unity is implied.
    connectivity : int
        The neighbourhood/connectivity considered when determining the surface
        of the binary objects. This value is passed to
        `scipy.ndimage.generate_binary_structure` and should usually be :math:`> 1`.
        Note that the connectivity influences the result in the case of the Hausdorff distance.

    Returns
    -------
    hd : float
        The symmetric Hausdorff Distance between the object(s) in ```result``` and the
        object(s) in ```reference```. The distance unit is the same as for the spacing of 
        elements along each dimension, which is usually given in mm.

    See also
    --------
    :func:`hd`

    Notes
    -----
    This is a real metric. The binary images can therefore be supplied in any order.
    """
    hd1 = __surface_distances(result, reference, voxelspacing, connectivity)
    hd2 = __surface_distances(reference, result, voxelspacing, connectivity)
    hd95 = numpy.percentile(numpy.hstack((hd1, hd2)), 95)
    return hd95


def get_image_objects(Img_path):
    
    img = nib.load(Img_path)
    img_obj = img.get_fdata()
    img_affine = img.affine
    img_header = img.header

    return img_obj, img_affine, img_header

def is_empty(img):
    
    if np.sum(np.sum(np.sum(img))) ==0:

        return True 