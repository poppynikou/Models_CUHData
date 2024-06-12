import os
import shutil 
import pandas as pd
import nibabel as nib 
import numpy as np 
import scipy.linalg

class PSM():

    def __init__(self, base_path, Patient_No, no_cps, anonymisation_key_path, results_path):
        
        '''self.Patient_No == HN_0 for example'''
        '''self.base_path == /home/pnikou/Documents/Manc_Data/'''

        self.base_path = base_path
        self.Patient_No = Patient_No
        self.no_cps = no_cps
        self.anonymisation_key_path = anonymisation_key_path
        self.structure_list = ['BRAINSTEM', 'CORD', 'CTVHIGH', 'CTVLOW', 'CTVMEDIUM', 'LPAROTID', 'RPAROTID']
        self.reg_resample = 'reg_resample'
        self.atlas_path = '/cluster/project7/HN_RT/CHRISTIE_HN/average_pCT.nii.gz'
        self.results_path = results_path
        return 
    
    def get_CBCT_relative_timepoints(self):
        
        # since we did not run all the registrations because of time constraints
        # i will just search for the CBCTs which exist in the folder 

        patient_path = self.base_path + 'UCLHMODELSPACE_REGS/' + str(self.Patient_No)
        CBCT_relative_timepoints = []
        for folder in os.listdir(patient_path):
            if folder.startswith('CBCT'):
                path = patient_path + '/' + folder + '/cpp_CBCT.nii.gz'
                if os.path.exists(path):
                    CBCT_relative_timepoints.append(int(folder.replace('CBCT_', '')))

        CBCT_relative_timepoints = list(np.sort(CBCT_relative_timepoints))
        return CBCT_relative_timepoints
    

    def get_CBCT_time_points(self):
        training_time_points = self.get_CBCT_relative_timepoints()
        return training_time_points

    def set_training_time_points(self, training_time_points = []):
        if len(training_time_points) == 0:
            self.training_time_points = self.get_CBCT_relative_timepoints()
        else:
            self.training_time_points = training_time_points
        
        print(self.training_time_points)


    def set_testing_time_points(self, testing_time_points = []):
        if len(testing_time_points) == 0:
            self.testing_time_points = np.arange(0, self.training_time_points[-1]+1)
            self.save_test_points = np.arange(0, self.training_time_points[-1]+1)
        else:
            self.testing_time_points = np.arange(0, self.training_time_points[-1]+1)
            self.save_test_points = testing_time_points
        

    def get_img_objects(self, img_path):
        
        img_obj = nib.load(img_path)
        img_data = np.array(img_obj.get_fdata())
        img_affine = img_obj.affine
        img_header = img_obj.header

        return img_data, img_affine, img_header
    
    def set_reference_data(self):

        first_training_time_point = self.training_time_points[0]

        ref_cpp_path = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.Patient_No)+ '/CBCT_' + str(first_training_time_point) + '/cpp_CBCT.nii.gz'
        img_data, self.ref_affine, self.ref_header = self.get_img_objects(ref_cpp_path)
        img_shape = np.shape(img_data)
        self.cpp_shape = [img_shape[0], img_shape[1], img_shape[2], img_shape[4]]
        self.no_voxels = self.cpp_shape[0]* self.cpp_shape[1] * self.cpp_shape[2] * self.cpp_shape[3]

    def set_derivative_reference_data(self):
        
        ref_cpp_path = self.results_path + '/' + str(self.Patient_No) + '/cpp_0.nii.gz'
        img_data, self.ref_affine, self.ref_header = self.get_img_objects(ref_cpp_path)
        img_shape = np.shape(img_data)
        self.cpp_shape = [img_shape[0], img_shape[1], img_shape[2], img_shape[4]]
        self.no_voxels = self.cpp_shape[0]* self.cpp_shape[1] * self.cpp_shape[2] * self.cpp_shape[3]


    def get_training_data(self):

        # creates an array to store cpp data
        training_cpp = np.empty(shape=((len(self.training_time_points), ) + tuple(self.cpp_shape)), dtype=np.float32)
        # print(coordinates.shape)
        # gets all the data from all the CBCTs in that patient series

        for i, training_time_point in enumerate(self.training_time_points):
            
            cpp_directory = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.Patient_No)+ '/CBCT_' + str(training_time_point) + '/cpp_CBCT.nii.gz'
            cpp = nib.load(cpp_directory)
            cpp = np.squeeze(cpp.get_fdata())
            training_cpp[i] = cpp

        return training_cpp
    
    def get_derivative_training_data(self):

        # creates an array to store cpp data
        training_cpp = np.empty(shape=((len(self.training_time_points), ) + tuple(self.cpp_shape)), dtype=np.float32)
        #print(coordinates.shape)
        # gets all the data from all the CBCTs in that patient series
        
        for i, training_time_point in enumerate(self.training_time_points):
            
            # find a path 
            cpp_directory = self.base_path + '/PSM_noLOO_T_All_CPS_4_ModelSpace/' + str(self.Patient_No) + '/cpp_' + str(training_time_point) + '.nii.gz'
            cpp = nib.load(cpp_directory)
            cpp = np.squeeze(cpp.get_fdata())
            training_cpp[i] = cpp


        return training_cpp

    
    def Calc_BSplineCoeffs(self, timepoints, numcp):

        '''
        function fits normal BSplines 
        param = parameters
        numcp = number of control points
        ''' 

        param = timepoints
        minval = timepoints[0]
        maxval = timepoints[-1]
        
        # creates a storage array
        bspline_coeffs = np.zeros(shape = (len(param), numcp+1))
        
        # raises error if the number of control points are less than four
        # you need four in order to fit cubic BSpline
        if numcp <4:
            ValueError('must have at least four control points')

        ## add scenarios where you dont have to calculate anything 

        # define interval in control point space, and then their positions
        cpspace = (maxval-minval)/(numcp-3)
        
        # you need the plus one to create 4 cpp
        cplocs = np.arange(minval-cpspace, maxval+cpspace+1, cpspace)
        
        #print(cplocs)
        s = [(param[i]-cplocs[cplocs <= param[i]][-1])/cpspace for i in range(len(param))]
        #print(s)
        index = [np.sum([cplocs <= param[i]]) for i in range(len(param))]
        #print(index)

        # calculate the matrix of BSpline coefficients 
        for i in range(len(param)):

            j = index[i]

            x = np.array([np.power(s[i],3), np.power(s[i],2), s[i], 1]).reshape((1,4))

            constants = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,0,3,0],[1,4,1,0]])

            bspline_coeffs[i, j-2:j+2] = 1/6 * np.matmul(x, constants)
        
        bspline_coeffs = bspline_coeffs[:, 0:numcp]     

        BSplineCoeffs = np.array(bspline_coeffs, dtype = 'float32')

        return BSplineCoeffs
    
    def Calc_ControlPointGrid(self, gt_velocity_fields, BSplineCoeffs):

        '''
        This function fits the BSpline Functions 
        and outputs the predictions of the velocity fields at the specified test_time_points

        gt_velocity_fields: ground truth velocity fields obtained from the registrations. Numpy array of size (no_of_time_points, x_dim, y_dim, z_dim, 3)
        gt_BS_cps: ground truth BSpline control point positions. Numpy array of size (no_of_time_points, no_control_points_in_time)
        test_time_points: Array of time points from 0 to the end of treatment. Day 0 is the day of the first CBCT
        cpp_shape: shape of the velocity field obtained from deformable registraion. Tuple
        numcp: Number of control points to fit the BSpline functions to. Float


        save_BSpline_coeffs: Boolean. Whether to save the BSpline coefficients for modelling later on. 
        Patient: Boolean or Integer. If save_BSpline_coeffs is True, then Patient must be an interger.

        returns:
        Predicted velocity field at each test time point. Numpy array of shape (no_of_test_time_points, x_dim, y_dim, z_dim, 3)

        '''
        
        # no of training time points
        no_time_points = np.shape(gt_velocity_fields)[0]

        # first, we reshape the cpp_files
        # new shape (no_time_points, no_voxels)
        gt_velocity_fields = np.reshape(gt_velocity_fields, newshape=(no_time_points, self.no_voxels))

        #### ----- Step 1 ------ ###
        # we use Ax = B, to find x = A\B. 
        # A = gt_velocity_fields
        # B = gt_BS_coeff

        # calculate the control point positions
        # this has shape (numcp, no_voxels)
        control_point_grid = scipy.linalg.lstsq(BSplineCoeffs, gt_velocity_fields)[0]

        return control_point_grid
    
    def Calc_NewTransformation(self, BSpline_coeffs, Control_Point_Grid):
    
        # we use Ax = B, to find B
        # A = test_BSpline_coefs
        # x = model_velocity_field

        BSpline_Coefficients_shape = np.shape(BSpline_coeffs)

        # final predicted velocity field 
        # this has shape (no_test_time_points, no_voxels)
        final_cpp = np.matmul(BSpline_coeffs, Control_Point_Grid)

        # reshape to be the shape of a velocity field again
        test_cpp = np.reshape(final_cpp, newshape=(BSpline_Coefficients_shape[0], self.cpp_shape[0], self.cpp_shape[1], self.cpp_shape[2], self.cpp_shape[3]))
        del final_cpp

        return test_cpp
   
    def fit_SM(self):

        BSplineCoeffs = self.Calc_BSplineCoeffs(self.training_time_points, self.no_cps)
        training_cpp = self.get_training_data()
        self.controlpointgrid = self.Calc_ControlPointGrid(training_cpp, BSplineCoeffs)
        return self.controlpointgrid
    
    def fit_derivative_SM(self):
        BSplineCoeffs = self.Calc_BSplineCoeffs(self.training_time_points, self.no_cps)
        training_cpp = self.get_derivative_training_data()
        self.controlpointgrid = self.Calc_ControlPointGrid(training_cpp, BSplineCoeffs)
        return self.controlpointgrid
    
    def save_temporal_CPG(self):

        CPG = np.reshape(self.controlpointgrid, newshape=(4, self.cpp_shape[0], self.cpp_shape[1], self.cpp_shape[2], self.cpp_shape[3]))
        NewNiftiObj = nib.Nifti1Image(CPG, self.ref_affine, None)
     
        if not os.path.exists(self.results_path + '/' + str(self.Patient_No)):
            os.mkdir(self.results_path + '/' + str(self.Patient_No))
        directory = self.results_path + '/' + str(self.Patient_No) + '/TEMP_CPPGRID_T35.nii.gz'
        nib.save(NewNiftiObj, directory)

        
    def test_SM(self):

        BSplineCoeffs = self.Calc_BSplineCoeffs(self.testing_time_points, self.no_cps)
        self.testing_cpp = self.Calc_NewTransformation(BSplineCoeffs, self.controlpointgrid)
        return self.testing_cpp

    def save_SM(self):
        
        for save_time_point in self.save_test_points:

            cpp_new = np.squeeze(self.testing_cpp[save_time_point, :, :, :, :])

            # need to save the cpp file
            #cpp_new = np.reshape(cpp_new, newshape = (cpp_shape[0], cpp_shape[1], cpp_shape[2], 1, cpp_shape[3]))
            cpp_new = np.expand_dims(cpp_new, axis = 3)
            # create new nifti object
            new_niftiobj = nib.Nifti1Image(cpp_new,  self.ref_affine,  self.ref_header)

            # save 
            patient_folder_path = self.results_path + '/' + str(self.Patient_No)
            if not os.path.exists(patient_folder_path):
                os.mkdir(patient_folder_path)
            file_name = patient_folder_path + '/cpp_' + str(save_time_point) + '.nii.gz'
            nib.save(new_niftiobj, file_name)
    
    def resampleBINImg(self, ref_img, float_img, transformation, resampled_img):
    
        command = self.reg_resample + ' -ref ' + ref_img +  ' -flo ' + float_img + ' -trans ' + transformation + ' -res ' + resampled_img + ' -inter 1 -omp 12 -pad 0'
        os.system(command)

    
    def resample_GT_Model(self):

        ref_img = self.atlas_path
        #make this
        new_folder =self.base_path + '/Masks/' + str(self.Patient_No) + '/atlas/'
        if not os.path.exists(new_folder):
            os.mkdir(new_folder)
        
        for structure in self.structure_list:
            float_img = self.base_path + '/Masks/' + str(self.Patient_No) + '/pCT/BIN_' + str(structure) + '.nii.gz'
            
            if os.path.exists(float_img):
                transformed_img = new_folder + '/BIN_' + str(structure) + '.nii.gz'
                if not os.path.exists(transformed_img):

                    model_cpp = self.base_path + '/Masks/' + str(self.Patient_No) + '/T_model.nii.gz'
                    self.resampleBINImg(ref_img, float_img, model_cpp, transformed_img)


    def resample_RTSTRUCTs(self, testing_time_point = None):
        ref_img = self.atlas_path
        for structure in self.structure_list:

            float_img = self.base_path + '/Masks/' + str(self.Patient_No) + '/atlas/BIN_' + str(structure) + '.nii.gz'
            
            if os.path.exists(float_img):

                if testing_time_point == None:
                    
                    for testing_time_point_index in self.testing_time_points:

                        model_cpp = self.results_path + '/' + str(self.Patient_No) + '/cpp_' + str(testing_time_point_index) + '.nii.gz'
                        transformed_img = self.results_path + '/' + str(self.Patient_No) + '/BIN_' + str(structure) + '_T_' + str(testing_time_point_index) +'.nii.gz'
                        self.resampleBINImg(ref_img, float_img, model_cpp, transformed_img)

                else:
                    model_cpp = self.results_path + '/' + str(self.Patient_No) + '/cpp_' + str(testing_time_point) + '.nii.gz'
                    transformed_img = self.results_path + '/' + str(self.Patient_No) + '/BIN_' + str(structure) + '_T_' + str(testing_time_point) +'.nii.gz'

                    self.resampleBINImg(ref_img, float_img, model_cpp, transformed_img)


    def resample_GT_RTSTRUCTs(self, training_time_point = None):

        ref_img = self.atlas_path
        for structure in self.structure_list:

            float_img = self.base_path + '/Masks/' + str(self.Patient_No) + '/atlas/BIN_' + str(structure) + '.nii.gz'
            
            if os.path.exists(float_img):
                transformed_img = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.Patient_No)+ '/CBCT_' + str(training_time_point) + '/BIN_' + str(structure) + '.nii.gz'
                #if not os.path.exists(transformed_img):
                if training_time_point != None:
                    model_cpp = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.Patient_No)+ '/CBCT_' + str(training_time_point) + '/cpp_CBCT.nii.gz'
                    self.resampleBINImg(ref_img, float_img, model_cpp, transformed_img)

                if training_time_point == None:
                    for T_time_point in self.training_time_points:
                        transformed_img = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.Patient_No)+ '/CBCT_' + str(T_time_point) + '/BIN_' + str(structure) + '.nii.gz'
                        #if not os.path.exists(transformed_img):
                        model_cpp = self.base_path + '/UCLHMODELSPACE_REGS/' + str(self.Patient_No)+ '/CBCT_' + str(T_time_point) + '/cpp_CBCT.nii.gz'
                        self.resampleBINImg(ref_img, float_img, model_cpp, transformed_img)
                