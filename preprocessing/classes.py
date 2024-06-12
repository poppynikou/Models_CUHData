import os
import shutil 
import pandas as pd
import nibabel as nib 
import numpy as np 
from functions import *
import scipy.ndimage as ndimage
from datetime import datetime, date

class GSTTData():
    
    def __init__(self, base_path, atlas_path):

        self.base_path = base_path
        self.reg_transform = 'reg_transform'
        self.reg_average = 'reg_average'
        self.reg_aladin = 'reg_aladin'
        self.reg_resample = 'reg_resample'
        self.reg_f3d = 'reg_f3d'
        self.atlas_path = atlas_path
        
    def test__filenotexsits(self, file_path):
        if not os.path.exists(file_path):
            return True 
    
    def get_img_objects(self, img_path, data_type = np.float32):
    
        img_obj = nib.load(img_path)
        img_data = np.array(img_obj.get_fdata(), dtype = data_type)
        img_affine = img_obj.affine
        img_header = img_obj.header

        return img_data, img_affine, img_header

class PatientData(GSTTData):

    def __init__(self, PatientID, base_path):
        self.PatientID = PatientID
        GSTTData.__init__(self, base_path, '')

    def get_PatientNo(self):
        
        return self.PatientID
    
    def get_dates(self):
        '''
        This function assumes that each image is stored in the folder w/ naming convention:

        CBCT_yyyymmdd - for CBCT images
        CT_yyyymmdd - for CT images 

        This function only lists the dates of the CBCT images.
        returns: yyyymmdd
        '''
        # function to find folder names of the CBCTs in the patient directory
        directory = self.base_path + '/' + self.PatientID + '/'
        ignore = ['CBCT_GROUPWISE', 'pCT', 'CBCT_pCT', 'T_model.nii.gz', 'MASKED_pCT.nii.gz']
        dates = [CBCT.replace('CBCT_', '') for CBCT in os.listdir(directory) if CBCT not in ignore]
        self.CBCT_dates = [date for date in list(dates) if date[0:2] != 'CT']

        return self.CBCT_dates


    def get_CBCT_relative_timepoints(self):
        '''
        This assumes that the first time points is the first CBCT in the series
        '''

        # create an array to store the time points
        time_points = []
        
        # find the date for the first CBCT date 
        day_zero = self.CBCT_dates[0]
        date_zero = date(int(day_zero[0:4]), int(day_zero[4:6]),  int(day_zero[6:8]))

        # loop through dates and get the time points
        for CBCT_date in self.CBCT_dates:
            
            CBCT_date = CBCT_date
            date_i = date(int(CBCT_date[0:4]), int(CBCT_date[4:6]),  int(CBCT_date[6:8]))
            time_points.append((date_i-date_zero).days) 

        return time_points
        

        #directory = self.base_path + '/' + self.PatientID + '/CBCT_GROUPWISE/affine_0/'
        #time_points = [int(CBCT.replace('CBCT_', '').replace('.nii.gz', '')) for CBCT in os.listdir(directory)]
        #return time_points


        
    

    def get_pCT_date(self):

        # function to find the pCT date for a specific patient 
        directory = self.base_path + '/' + self.PatientID + '/'
        # list all contents of the directory 
        dates = os.listdir(directory)
        # find all folders which start with 'CT_'
        CT_dates = [date[-8:]  for date in list(dates) if date[0:2] == 'CT']
        CT_dates = np.asarray(CT_dates, dtype=object)
        # get the dates times, so that you can get the one which is earliest in time 
        CT_datetimes = [datetime(int(CT[0:4]), int(CT[4:6]), int(CT[6:8])) for CT in CT_dates]
        pCT_datetime = [min(CT_datetimes)]

        Boolean = np.in1d(CT_datetimes, pCT_datetime)

        pCT_date = CT_dates[Boolean][0]
        
        return pCT_date
    
    def get_patient_folder(self):

        return self.base_path + '/' + self.PatientID + '/'

    

class Image(GSTTData):
    
    def __init__(self, PatientNo, base_path, pCT_date):

        self.PatientNo = PatientNo
        self.pCT_date = pCT_date
        GSTTData.__init__(self, base_path, '')
        return
    
    def read_meta_info(self, csv_file):

        meta_info =  pd.read_csv(csv_file, header=0)
        self.patient_masking_info = meta_info.loc[(meta_info['PatientNo']==self.PatientNo)]

    
    def flip_img_Bool(self, flip_record):

        path_to_CT = self.base_path + str(self.PatientNo) + '/CT_' + str(self.pCT_date) + '/CT.nii.gz'
        if os.path.exists(self.base_path + str(self.PatientNo) + '/CT_' + str(self.pCT_date) + '/STRUCTURES/BIN_CTVHIGH_OG.nii.gz'):
            path_to_CTV = self.base_path + str(self.PatientNo) + '/CT_' + str(self.pCT_date) + '/STRUCTURES/BIN_CTVHIGH_OG.nii.gz'
        else:
            path_to_CTV = self.base_path + str(self.PatientNo) + '/CT_' + str(self.pCT_date) + '/STRUCTURES/BIN_CTVHIGH.nii.gz'
        path_to_CORD = self.base_path + str(self.PatientNo) + '/CT_' + str(self.pCT_date) + '/STRUCTURES/BIN_CORD.nii.gz'

        # import imgs which you need
        CT_img, _, _ = self.get_img_objects(path_to_CT)
        CTV_img, _, _ = self.get_img_objects(path_to_CTV)
        
        CORD_img, _, _ = self.get_img_objects(path_to_CORD)
        del _ 

        # shape of the images 
        x_shape, _, _ = np.shape(CT_img)
        del _

        # sum the voxel values in each slice
        sum_slices = np.sum(np.sum(CORD_img, axis =0), axis=0)
        # first the first non zero value, this corresponds to the first axial slice in which non zero values are stored
        z_coord = np.where(np.array(sum_slices) != 0)[0][-1]

        # the number of non zero voxels in the indexed slice
        no_nonzero_in_indexslice = len(np.where(CORD_img[:,:,z_coord]!=0)[0])
        # the middle non zero voxel in the indexed slice
        middle_nonzero_index = int(no_nonzero_in_indexslice/2)
        # coordinated of the middle of the spinal cord 
        x_coord = np.where(CORD_img[:,:,z_coord]!=0)[0][middle_nonzero_index]

        # find number of voxels in left and right hand side
        # to determine whether you need to flip images 
        CTV_left = int(np.sum(np.sum(np.sum(CTV_img[0:x_coord,:,:], axis = 2), axis=1)))
        CTV_right = int(np.sum(np.sum(np.sum(CTV_img[x_coord:x_shape,:,:], axis = 2), axis=1)))
        Difference_in_Voxels = np.abs(CTV_left - CTV_right)
        
        # flip so that all patients have majority of CTV of the right hand side 
        if CTV_left > CTV_right:
            flip_record.write(str(self.PatientNo) + ' ' + str(CTV_left) + ' ' + str(CTV_right) + ' ' + str(Difference_in_Voxels) + ' True \n')
            return True
        
        else:
            flip_record.write(str(self.PatientNo) + ' ' + str(CTV_left) + ' ' + str(CTV_right) + ' ' + str(Difference_in_Voxels) + ' False \n')
            return False
        
    def flip_img(self, img_path, new_img_path = False):
    
        img_data, img_affine, img_header = self.get_img_objects(img_path)
            
        img_flipped = np.flip(img_data, axis = 0)
        NewNiftiObj = nib.Nifti1Image(img_flipped, img_affine, img_header)
        if new_img_path != False:
            os.rename(img_path, img_path[:-7] + '_OG.nii.gz')
            nib.save(NewNiftiObj, img_path)
        else:
            nib.save(NewNiftiObj, img_path)
        
        
    def rename_parotid(self, file_path, name, flip = False):
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
        elif name.__contains__('PAROTIDIPSI'):
            if flip:
                new_name=new_name.replace('PAROTIDIPSI','CONTRAPAROTID')
            else:
                new_name=new_name.replace('PAROTIDIPSI','IPSIPAROTID')
            source = file_path + '/' + name
            destination = file_path + '/' + new_name
            os.rename(source, destination)
        elif name.__contains__('PAROTIDCONTRA'):
            if flip:
                new_name=new_name.replace('PAROTIDCONTRA','IPSIPAROTID')
            else:
                new_name=new_name.replace('PAROTIDCONTRA','CONTRAPAROTID')
            source = file_path + '/' + name
            destination = file_path + '/' + new_name
            os.rename(source, destination)


    def convert_to_float(self, img_path, new_img_path = False):

        img_data, img_affine, img_header = self.get_img_objects(img_path)
        img_data_copy = np.array(img_data, dtype = np.float32)

        NewNiftiObj = nib.Nifti1Image(img_data_copy, img_affine, img_header)
        NewNiftiObj.set_data_dtype('float32')
        if new_img_path != False:
            nib.save(NewNiftiObj, new_img_path)#
        else:
            nib.save(NewNiftiObj, img_path)#

    def clip_HU(self, img_path, new_img_path = False):
    
        img_data, img_affine, img_header = self.get_img_objects(img_path)
        img_data_copy = img_data.copy()
        img_data_copy[img_data_copy>1500] = 1500
        
        NewNiftiObj = nib.Nifti1Image(img_data_copy, img_affine, img_header)
       
        if new_img_path != False:
            nib.save(NewNiftiObj, new_img_path)#
        else:
            nib.save(NewNiftiObj, img_path)#

   

    def get_couch_slice(self, img_path):
    
        img_data, _, _ = self.get_img_objects(img_path, data_type= np.int8)
        del _
        print(np.sum(img_data))
        indexes = np.where(img_data==1)
        print(indexes)
        couch_slice = np.amin(indexes[1]) - 20  

        return couch_slice
    
    def mask_CT(self, img_path):
    
        masked_CT_path = img_path[:-9] + 'MASKED_CT.nii.gz'
        #print(masked_CT_path)
        # read in data 
        img_data, img_affine, img_header = self.get_img_objects(img_path)
        
        # try this 
        # 360 was base
        couch_slice = 360
        img_data[:,couch_slice:,:] = np.NaN

        NewNiftiObj = nib.Nifti1Image(img_data, img_affine, img_header)
        nib.save(NewNiftiObj, masked_CT_path)

    def get_vertabrae_segmentations(self, input_img):
        
        # create output folder 
        output_folder = input_img[:-9] + '/TOTAL_SEGMENTATOR_STRUCTURES/'
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
             
        if not os.path.exists(output_folder + '/vertebrae_T4.nii.gz'):
            os.system('TotalSegmentator -i ' + str(input_img) + ' -o ' + str(output_folder) + ' --fast')

    def get_cropping_slice(self, vertabrae_segmentation_path):
        
        img_data, _, _ = self.get_img_objects(vertabrae_segmentation_path)
        del _ 

        binary_img = np.sum(np.sum(img_data, axis=0), axis =0)
        index = np.where(binary_img!=0)[0][0]
        cropping_slice = int(index)

        return cropping_slice
    
    def mask_anatomy(self):

        
        img_path = self.base_path + '/' + self.PatientNo + '/CT_' + str(self.pCT_date) + '/CT.nii.gz'

        masked_img_path = img_path[:-9] + 'MASKED_CT.nii.gz'
        vertabrae_segmentation = self.base_path + '/' + self.PatientNo + '/CT_' + str(self.pCT_date) + '/TOTAL_SEGMENTATOR_STRUCTURES/vertebrae_T4.nii.gz'
        masked_img = img_path[:-9] + 'atlas_MASKED_CT.nii.gz' 

        # read in data 
        img_data, img_affine, img_header = self.get_img_objects(masked_img_path)
        img_data_copy = np.array(img_data.copy(), dtype = np.float32)
        (_,y,_) = np.shape(img_data_copy)
        del(_)

        self.get_vertabrae_segmentations(img_path)
        min_z = self.get_cropping_slice(vertabrae_segmentation)

        img_data_copy[:,:,0:min_z] = np.NaN  

        # overide the CT 
        newNiftiObj = nib.Nifti1Image(img_data_copy, img_affine, img_header)
        newNiftiObj.set_data_dtype('float32')
        nib.save(newNiftiObj, masked_img)
                
        

    def create_circular_mask(self, h, w, center=None, radius=None):
        
        if center is None: # use the middle of the image
            center = (int(w/2), int(h/2))
        if radius is None: # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w-center[0], h-center[1]) 

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

        mask = dist_from_center <= radius
        return mask    
    


    def get_body_segmentation(self, input_img, output_segmentation):


        os.system('TotalSegmentator -i ' + str(input_img) + ' -o ' + str(output_segmentation) + ' -ta body ')

        os.remove(output_segmentation + '/body_extremities.nii.gz')
        os.remove(output_segmentation + '/body_trunc.nii.gz')
        os.remove(output_segmentation + '/skin.nii.gz')

    def mask_CBCT(self, input_CBCT_path, slices_above, mask_couch = True):
        
        # reads in CBCT 
        img_data, img_affine, img_header = self.get_img_objects(input_CBCT_path)
        img_data = np.array(img_data, dtype= np.float32)
        l, w, h = np.shape(img_data)

        # mask for the circular field of view 
        mask = np.zeros(shape = np.shape(img_data))
        
        # creates a circular mask with radius = h/2 
        # assumes the slices are square
        for slice in np.arange(0, np.shape(img_data)[2]):
            mask[:,:, slice] = self.create_circular_mask(l,w)

        img_data[mask == 0] = np.NaN
        img_data[:,:,0:slices_above] = np.NaN
        img_data[:,:,h-slices_above:h] = np.NaN

        segmentation_path = input_CBCT_path[:-11] 
        bodysegmentation_path = segmentation_path + 'body.nii.gz'
        if not os.path.exists(bodysegmentation_path):
            self.get_body_segmentation(input_CBCT_path, segmentation_path)
        # reads in body segmentation 
        body_img_data, img_affine, img_header = self.get_img_objects(bodysegmentation_path, data_type= np.int8)
        #body_img_data = ndimage.binary_dilation(body_img_data, iterations=15)
        
        if mask_couch:
            # this is in a try loop since some CBCTs dont have a couch
            # mask for the couch 
            couch_mask = np.zeros(shape = (l,w,h))
            boolean = (body_img_data == 0)  & (img_data >= -1000) & (img_data != np.NaN)
            couch_mask[boolean] = 1

            # erode the couch mask a few times
            # usually 9 
            # but keep decreasing a little bit incase there are patients where the couch is small 
        
            for i in np.arange(9,15)[::-1]:
                couch_mask_new = ndimage.binary_erosion(couch_mask, iterations=i)
                if np.sum(couch_mask_new) != 0:
                    
                    # find starting slice of the couch 
                    index = np.where(couch_mask_new == 1)
                    yslice = np.amin(index[1])   
                    # mask out the couch 
                    img_data[:,yslice:,:] = np.NaN
                    break

            
        else:
            pass

        masked_CBCT_path = input_CBCT_path[:-11] + '/MASKED_' + input_CBCT_path[-11:]
        NewNiftiObj = nib.Nifti1Image(img_data, img_affine, img_header)
        NewNiftiObj.set_data_dtype(np.float32)
        nib.save(NewNiftiObj, masked_CBCT_path)


class GroupwiseRegs(GSTTData):

    def __init__(self, PatientNo, pCT_date, CBCT_dates, CBCT_relative_timepoints, no_itterations, base_path, log_path):
        self.PatientNo = PatientNo 
        self.no_itterations = no_itterations
        GSTTData.__init__(self, base_path, '')
        self.pCT_date = pCT_date
        self.CBCT_dates = CBCT_dates
        self.CBCT_relative_timepoints = CBCT_relative_timepoints
        return

    def refactor(self):

        # function that makes sure the folders are correctly set-up
        self.results_folder = str(self.base_path) + '/' + str(self.PatientNo) + '/CBCT_GROUPWISE/'
        print(self.results_folder)
        if not os.path.exists(self.results_folder):
            os.mkdir(self.results_folder)
        
        for no_itteration in np.arange(0, self.no_itterations+1):
            affine_folder = self.results_folder + 'affine_' + str(no_itteration)
            if not os.path.exists(affine_folder):
                os.mkdir(affine_folder)

        self.postprocessing_path = self.results_folder + 'postprocessing/'
        if not os.path.exists(self.postprocessing_path):
            os.mkdir(self.postprocessing_path)

        self.CBCT_pCT_path = str(self.base_path) + '/' + str(self.PatientNo) + '/CBCT_pCT/'
        if not os.path.exists(self.CBCT_pCT_path):
            os.mkdir(self.CBCT_pCT_path)


        for index, CBCT_timepoint in enumerate(self.CBCT_relative_timepoints):
                # cropped imgs 
            source = str(self.base_path) + '/' + str(self.PatientNo) + '/CBCT_' + str(self.CBCT_dates[index]) + '/MASKED_CBCT.nii.gz'
            destination = self.results_folder + 'affine_0/CBCT_' + str(CBCT_timepoint) + '.nii.gz'
            shutil.copy(source, destination)

        for index, CBCT_timepoint in enumerate(self.CBCT_relative_timepoints):
            # full imgs 
            source = str(self.base_path) + '/' + str(self.PatientNo) + '/CBCT_' + str(self.CBCT_dates[index]) + '/MASKED_CBCT.nii.gz'
            destination = self.postprocessing_path + '/CBCT_' + str(CBCT_timepoint) + '.nii.gz'
            shutil.copy(source, destination)
        

    
    def set_itteration(self, itteration):

        self.itteration = itteration
        
    def set_ref_img(self):

        # set the target/ reference img for that iteration 
        if self.itteration == 0:
            self.ref_img = str(self.base_path) + '/' + str(self.PatientNo) + '/CBCT_GROUPWISE/affine_0/CBCT_0.nii.gz'
        else:
            self.ref_img = self.results_folder + '/affine_' + str(self.itteration) + '/average_CBCT.nii.gz'
    
    def set_float_imgs(self):
        self.float_imgs = []
        for CBCT_timepoint in self.CBCT_relative_timepoints:
            self.float_imgs.append(self.results_folder + '/affine_' + str(self.itteration) + '/CBCT_'+str(CBCT_timepoint)+'.nii.gz')

    def set_rigidReg_affinematrixes(self):
    
        self.affine_matrixes = []
        for CBCT_timepoint in self.CBCT_relative_timepoints:
            self.affine_matrixes.append(self.results_folder + '/affine_' + str(self.itteration+1) + '/affine_'+str(CBCT_timepoint)+'.txt')

    def set_rigidReg_resampledimgs(self):

        self.rigidReg_resampled_imgs = []
        for CBCT_timepoint in self.CBCT_relative_timepoints:
            self.rigidReg_resampled_imgs.append(self.results_folder + '/affine_' + str(self.itteration +1) + '/res_CBCT_'+str(CBCT_timepoint)+'.nii.gz')

    def set_composed_affinematrixes(self):

        self.compaffine_matrixes = []
        for CBCT_timepoint in self.CBCT_relative_timepoints:
            self.compaffine_matrixes.append(self.results_folder + '/affine_' + str(self.itteration+1) + '/comp_affine_'+str(CBCT_timepoint)+'.txt')

    def set_resampledimgs(self):

        self.resampled_imgs = []
        for CBCT_timepoint in self.CBCT_relative_timepoints:
            self.resampled_imgs.append(self.results_folder + '/affine_' + str(self.itteration+1) + '/CBCT_'+str(CBCT_timepoint)+'.nii.gz')

    def set__itteration(self, itteration):

        self.set_itteration(itteration)
        self.set_ref_img()
        self.set_float_imgs()
        self.set_rigidReg_affinematrixes()
        self.set_rigidReg_resampledimgs()
        self.set_composed_affinematrixes()
        self.set_resampledimgs()
        


    def rigidGroupReg(self):
    
        for index in np.arange(0, len(self.CBCT_relative_timepoints)):

            float_img = self.float_imgs[index]
            affine_matrix = self.affine_matrixes[index]
            resampled_img = self.rigidReg_resampled_imgs[index]
                
            rigidReg(self.reg_aladin, self.ref_img, float_img, affine_matrix, resampled_img, RigOnly = True)

            
    
    def avgAffine(self):

        self.average_affine = self.results_folder + '/affine_' + str(self.itteration + 1) + '/average_affine.txt'

        avgAff(self.reg_average, self.average_affine, self.affine_matrixes)

    def invAffine(self):
        
        self.inv_affine = self.results_folder + '/affine_' + str(self.itteration + 1) + '/inv_average_affine.txt'

        invAff(self.reg_transform, self.ref_img, self.average_affine, self.inv_affine)

    def compAffine(self):

        for index in np.arange(0, len(self.CBCT_relative_timepoints)):
            
            affine_matrix = self.affine_matrixes[index]
            comp_matrix = self.compaffine_matrixes[index]

            compAff(self.reg_transform, self.ref_img, self.inv_affine, affine_matrix, comp_matrix)

    def resampleImages(self):

        for index in np.arange(0, len(self.CBCT_relative_timepoints)):
            
            float_img = self.float_imgs[index]
            resampled_img = self.resampled_imgs[index]
            transformation = self.compaffine_matrixes[index]
            
            resampleImg(self.reg_resample, self.ref_img, float_img, transformation, resampled_img)
    
    def avgImage(self):

        avg_img_path = self.results_folder + '/affine_' + str(self.itteration +1) + '/average_CBCT.nii.gz'

        img_data, _, _ = self.get_img_objects(self.resampled_imgs[0])
        img_shape = np.shape(np.array(img_data))
        del(img_data)
        del(_)

        # create an array of nans for storage
        Nan_map = np.empty(shape=(img_shape[0],img_shape[1],img_shape[2],len(self.resampled_imgs)))
        #del(img_shape)

        for index in np.arange(0,len(self.resampled_imgs)):

            # imports the data 
            img_data, img_affine, img_header = self.get_img_objects(self.resampled_imgs[index])

            #store
            Nan_map[:,:,:,index] = img_data

        del(img_data)

        # calculate the average image, ignoring the nans
        Masked_Average = np.empty(shape = (img_shape[0],img_shape[1],img_shape[2]))
        for slice in np.arange(0, img_shape[2]):
            Masked_Average[:,:,slice] = np.nanmean(Nan_map[:,:,slice,:], axis = 2)
        del(Nan_map)
        Average = Masked_Average.copy()
        del(Masked_Average)  

        # save 
        Avg_Niftiobj = nib.Nifti1Image(Average, img_affine, img_header)
        nib.save(Avg_Niftiobj, avg_img_path)


    def UpdateGroupSform(self):

        for CBCT_timepoint in self.CBCT_relative_timepoints:
            for itteration in np.arange(0, self.no_itterations):
                
                img = self.postprocessing_path + '/CBCT_' + str(CBCT_timepoint) + '.nii.gz'
                affine_matrix_path = self.results_folder + '/affine_' + str(itteration+1) + '/comp_affine_'+str(CBCT_timepoint)+'.txt'

                UpdSform(self.reg_transform, img, affine_matrix_path, img)

    def rigidpCTReg(self):

        ref_img = str(self.base_path) + '/' + str(self.PatientNo) + '/CT_'+str(self.pCT_date)+'/MASKED_CT.nii.gz'
        float_img = self.postprocessing_path + '/CBCT_0.nii.gz'
        transformation = self.CBCT_pCT_path + '/affine.txt'
        resampled_img = self.CBCT_pCT_path + '/CBCT_0.nii.gz'

        rigidReg(self.reg_aladin, ref_img, float_img, transformation, resampled_img, RigOnly = True)


    def UpdateSform(self):

        for CBCT_timepoint in self.CBCT_relative_timepoints:

            img_to_be_updated = self.postprocessing_path + '/CBCT_' + str(CBCT_timepoint) + '.nii.gz'
            affine_matrix_path = self.CBCT_pCT_path + '/affine.txt'
            updated_img = self.CBCT_pCT_path + '/CBCT_' + str(CBCT_timepoint) + '.nii.gz'

            UpdSform(self.reg_transform, img_to_be_updated, affine_matrix_path, updated_img)
    
class AtlasRegs(GSTTData):

    def __init__(self, PatientNo, pCT_date, base_path, atlas_path):
        GSTTData.__init__(self, base_path, atlas_path)
        self.PatientNo = PatientNo
        self.PatientPath = self.base_path + '/' + str(self.PatientNo) +'/'
        self.PatientCTPath = self.PatientPath + '/CT_'+str(pCT_date) + '/'
        self.atlas_path = atlas_path

        

    def refactor(self):
        self.ModelSpacePath = self.PatientCTPath + '/model_space/'
        if not os.path.exists(self.ModelSpacePath):
            os.mkdir(self.ModelSpacePath)

        UCLHRegsPath = self.base_path + '/UCLHMODELSPACE_REGS/'
        if not os.path.exists(UCLHRegsPath):
            os.mkdir(UCLHRegsPath)
        
        self.PatientUCLHRegsPath = UCLHRegsPath + str(self.PatientNo) +'/'
        if not os.path.exists(self.PatientUCLHRegsPath):
            os.mkdir(self.PatientUCLHRegsPath)
    

    def InitAlignment(self):

        float_img = self.PatientCTPath + 'atlas_MASKED_CT.nii.gz'
        affine_matrix = self.ModelSpacePath + 'InitAlignment_atlas.txt'
        resampled_img = self.ModelSpacePath + 'resampled_InitAlignment_pCT.nii.gz'

        rigidReg(self.reg_aladin, self.atlas_path, float_img, affine_matrix, resampled_img, RigOnly= True)
        os.remove(resampled_img)
        


        img_to_be_updated = self.PatientCTPath + 'atlas_MASKED_CT.nii.gz' 
        updated_img = self.ModelSpacePath + 'InitAlignment_pCT.nii.gz'
        UpdSform(self.reg_transform, img_to_be_updated, affine_matrix, updated_img)
        return

    def RigidReg(self):

        float_img = self.ModelSpacePath + 'InitAlignment_pCT.nii.gz'
        affine_matrix = self.ModelSpacePath + 'Rigid_atlas.txt'
        resampled_img = self.ModelSpacePath + 'resampled_pCT_atlas_rigid.nii.gz'

        rigidReg(self.reg_aladin, self.atlas_path, float_img, affine_matrix, resampled_img, RigOnly= True)
        os.remove(resampled_img)

        img_to_be_updated = self.ModelSpacePath + 'InitAlignment_pCT.nii.gz'
        updated_img = self.ModelSpacePath + 'Rigid_pCT.nii.gz'
        UpdSform(self.reg_transform, img_to_be_updated, affine_matrix, updated_img)

        return
    
    def AffineReg(self):

        float_img = self.ModelSpacePath + 'Rigid_pCT.nii.gz'
        affine_matrix = self.ModelSpacePath + 'Affine_atlas.txt'
        resampled_img =self.ModelSpacePath + 'resampled_pCT_atlas_affine.nii.gz'

        rigidReg(self.reg_aladin, self.atlas_path, float_img, affine_matrix, resampled_img, RigOnly= False)
        os.remove(resampled_img)
        
        
        img_to_be_updated = self.ModelSpacePath + 'Rigid_pCT.nii.gz'
        updated_img = self.ModelSpacePath + 'Affine_pCT.nii.gz'
        UpdSform(self.reg_transform, img_to_be_updated, affine_matrix, updated_img)
        
        return
    
    def DefReg(self):

        float_img = self.ModelSpacePath + 'Affine_pCT.nii.gz'
        resampled_img = self.ModelSpacePath + 'DEF_pCT.nii.gz'
        transformation = self.ModelSpacePath + 'cpp_pCT.nii.gz'
        if not os.path.exists(transformation):
            deformableReg(self.reg_f3d, self.atlas_path, float_img, resampled_img, transformation)

             
        return
        
    def Calc_Tatlas(self):
        
        # calculate deformation fields of affine matrixes
        input_transformation = self.ModelSpacePath + 'InitAlignment_atlas.txt'
        output_transformation = self.ModelSpacePath + 'InitAlignment_atlas.nii.gz'
        RigidToDeformation(self.reg_transform, self.atlas_path, input_transformation, output_transformation)
        input_transformation = self.ModelSpacePath + 'Rigid_atlas.txt'
        output_transformation = self.ModelSpacePath + 'Rigid_atlas.nii.gz'
        RigidToDeformation(self.reg_transform, self.atlas_path, input_transformation, output_transformation)
        input_transformation = self.ModelSpacePath + 'Affine_atlas.txt'
        output_transformation = self.ModelSpacePath + 'Affine_atlas.nii.gz'
        RigidToDeformation(self.reg_transform, self.atlas_path, input_transformation, output_transformation)

        transformation1 = self.ModelSpacePath + 'Rigid_atlas.nii.gz'
        transformation2 = self.ModelSpacePath + 'InitAlignment_atlas.nii.gz'
        output_transformation1 = self.ModelSpacePath + 'comp1.nii.gz'
        ComposeTransformations(self.reg_transform, self.atlas_path, transformation1, transformation2, output_transformation1)
        os.remove(transformation1)
        os.remove(transformation2)

        transformation1 = self.ModelSpacePath + 'Affine_atlas.nii.gz'
        transformation2 = self.ModelSpacePath + 'comp1.nii.gz'
        output_transformation2 = self.ModelSpacePath + 'comp2.nii.gz'
        ComposeTransformations(self.reg_transform, self.atlas_path, transformation1, transformation2, output_transformation2)
        os.remove(transformation1)
        os.remove(transformation2)

        transformation1 = self.ModelSpacePath + 'cpp_pCT.nii.gz'
        transformation2 = self.ModelSpacePath + 'comp2.nii.gz'
        output_transformation3 = self.PatientPath + 'T_model.nii.gz'
        ComposeTransformations(self.reg_transform, self.atlas_path, transformation1, transformation2, output_transformation3)
        os.remove(transformation2)
        

        return


    def ResampleImgs(self, CBCT_timepoints):

        # function which uses T_model to resample all patient images into the model space
        T_model = self.PatientPath + 'T_model.nii.gz'
        
        float_img = self.PatientCTPath + 'MASKED_pCT.nii.gz'
        resampled_img = self.PatientUCLHRegsPath + '/MASKED_pCT.nii.gz'
        resampleImg(self.reg_resample, self.atlas_path, float_img, T_model, resampled_img)
        
        for CBCT_timepoint in CBCT_timepoints:

            CBCT_RegsPath = self.PatientUCLHRegsPath + '/CBCT_' + str(CBCT_timepoint)
            if not os.path.exists(CBCT_RegsPath):
                os.mkdir(CBCT_RegsPath)
 
            float_img =  str(self.base_path) + '/' + str(self.PatientNo) + '/CBCT_pCT/CBCT_' + str(CBCT_timepoint) + '.nii.gz'
            resampled_img = CBCT_RegsPath + '/MASKED_CBCT.nii.gz'
            resampleImg(self.reg_resample, self.atlas_path, float_img, T_model, resampled_img)

        
class DefromableRegs(GSTTData):

    def __init__(self, PatientNo, base_path):
        GSTTData.__init__(self, base_path, '')
        self.PatientNo = PatientNo
        self.PatientUCLHRegsPath = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.PatientNo)

    def set__CBCTtimepoint(self, CBCT_timepoint):

        self.CBCT_timepoint = CBCT_timepoint

    def DefReg(self):
        
        float_img = self.PatientUCLHRegsPath + '/MASKED_CT.nii.gz'
        ref_img = self.PatientUCLHRegsPath + '/CBCT_' + str(self.CBCT_timepoint) + '/MASKED_CBCT.nii.gz'
        resampled_img = self.PatientUCLHRegsPath + '/CBCT_' + str(self.CBCT_timepoint) + '/DEF_CBCT.nii.gz'
        cpp = self.PatientUCLHRegsPath + '/CBCT_' + str(self.CBCT_timepoint) +'/cpp_CBCT.nii.gz'

        deformableReg(self.reg_f3d, ref_img, float_img, resampled_img, cpp)
