from utils import * 
import numpy as np 
import sys 

patient = sys.argv[1]
model_folder = sys.argv[2]
atlas = '/cluster/project7/HN_RT/CHRISTIE_HN/average_pCT.nii.gz'

# get pCT date 
pCT_folder = [folder for folder in os.listdir('/cluster/project7/HN_RT/CUH_HN/'+str(patient)) if folder[0:3] =='CT_'][0]

CBCT_dates = [file.replace('CBCT_', '') for file in os.listdir('/cluster/project7/HN_RT/CUH_HN/'+str(patient)) if file[0:4] == 'CBCT' and file not in ['CBCT_pCT', 'CBCT_GROUPWISE']]
timepoints = get_time_points(CBCT_dates)
timepoints = np.arange(0, timepoints[-1],1)


for timepoint in timepoints:
        
        # make paths 
        if not os.path.exists('/cluster/project7/HN_RT/CUH_HN/' + str(model_folder) + '_PatientSpace'):
            os.mkdir('/cluster/project7/HN_RT/CUH_HN/' + str(model_folder)+ '_PatientSpace')
        if not os.path.exists('/cluster/project7/HN_RT/CUH_HN/' + str(model_folder)+ '_PatientSpace'+'/'+str(patient)):
            os.mkdir('/cluster/project7/HN_RT/CUH_HN/' + str(model_folder)+ '_PatientSpace'+'/'+str(patient))

        # create inv_T_model
        ref_img = atlas
        input_transformation = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) + '/T_model.nii.gz'
        float_img = atlas
        inverse_transform = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) + '/inv_T_model.nii.gz'
        if not os.path.exists(inverse_transform):
            inv_Dff(ref_img, input_transformation, float_img, inverse_transform)


        # transport model to patient space
        input_path = '/cluster/project7/HN_RT/CUH_HN/' + str(model_folder) + '_ModelSpace/' + str(patient) + '/cpp_' + str(timepoint) + '.nii.gz'
        output_transformation = '/cluster/project7/HN_RT/CUH_HN/' + str(model_folder) + '_ModelSpace/' + str(patient) + '/dff_' + str(timepoint) + '.nii.gz'
        
        transported_transformation = '/cluster/project7/HN_RT/CUH_HN/' + str(model_folder) + '_PatientSpace/'+str(patient) + '/dff_'+ str(timepoint) + '.nii.gz'
        patient_to_modelspace_data = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) 
        CBCT_img_modelspace = atlas
        CT_img_patientspace = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) +'/'+str(pCT_folder)+'/MASKED_CT.nii.gz'
        Transport_Model_to_Patient(model_folder, input_path, transported_transformation, output_transformation,  patient_to_modelspace_data, CBCT_img_modelspace, CT_img_patientspace) 


        transported_transformation = '/cluster/project7/HN_RT/CUH_HN/' + str(model_folder) + '_PatientSpace/'+str(patient) + '/dff_'+ str(timepoint) + '.nii.gz'
        CT_img_patientspace = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) +'/'+str(pCT_folder)+'/CT.nii.gz'    
        Resampled_CT_img_patientspace = '/cluster/project7/HN_RT/CUH_HN/' + str(model_folder) + '_PatientSpace/'+str(patient) + '/Img_'+ str(timepoint) + '.nii.gz'
        resampleImg(CT_img_patientspace, CT_img_patientspace, transported_transformation, Resampled_CT_img_patientspace)




#### === THE REST NEEDS EDITING FOR NEW FOLDER STRUCTURE === ###




#structures = ['BRAINSTEM', 'CORD', 'CTVHIGH', 'CTVLOW', 'CTVMEDIUM', 'LPAROTID', 'RPAROTID']




'''
# transport cpp to patient space 
input_path = '/cluster/project7/HN_RT/CUH_HN/UCLHMODELSPACE_REGS/'+str(patient)+'/CBCT_'+str(timepoint)+'/cpp_CBCT.nii.gz'
transported_transformation = '/cluster/project7/HN_RT/CUH_HN/CPPs_PatientSpace/' + '/'+str(patient) + '/CBCT_' + str(timepoint) + '/dff_'+ str(timepoint) + '.nii.gz'
patient_to_modelspace_data = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) 
CBCT_img_modelspace = atlas
CBCT_img_patientspace = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) +'/'+str(pCT_folder)+'/MASKED_CT.nii.gz'
if not os.path.exists(transported_transformation):
    Transport_Model_to_Patient(model, input_path, transported_transformation, patient_to_modelspace_data, CBCT_img_modelspace, CBCT_img_patientspace) 


for structure in structures:

    float_structure = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) +'/'+str(pCT_folder)+'/STRUCTURES/BIN_' + str(structure) + '.nii.gz'


    if os.path.exists(float_structure):
        
        
        ref_img = '/cluster/project7/HN_RT/CUH_HN/'+str(patient) +'/'+str(pCT_folder)+'/MASKED_CT.nii.gz'

        # resample using cpp 
        # gt 
        transformation = '/cluster/project7/HN_RT/CUH_HN/CPPs_PatientSpace/' + '/'+str(patient) + '/CBCT_' + str(timepoint) + '/dff_'+ str(timepoint) + '.nii.gz'
        resampled_img_cpp = '/cluster/project7/HN_RT/CUH_HN/CPPs_PatientSpace/' + '/'+str(patient) + '/CBCT_' + str(timepoint) + '/BIN_' + str(structure) + '.nii.gz'
        resampleBINImg(ref_img, float_structure, transformation, resampled_img_cpp)

        # resample using model 
        # test case
        transformation = '/cluster/project7/HN_RT/CUH_HN/' + str(model) + '/'+str(patient) + '/CBCT_' + str(timepoint) + '/dff_'+ str(timepoint) + '.nii.gz'
        resampled_img_model = '/cluster/project7/HN_RT/CUH_HN/' + str(model) + '/'+str(patient) + '/CBCT_' + str(timepoint) + '/BIN_' + str(structure) + '.nii.gz'
        resampleBINImg(ref_img, float_structure, transformation, resampled_img_model)
        img_obj, img_affine, img_header = get_image_objects(resampled_img_model)
        if img_obj.dtype != np.int8:
            save_as_int(resampled_img_model)
        del img_obj
        del img_affine
        del img_header


        gt_binary_obj, _, _ = get_image_objects(resampled_img_cpp)
        gt_binary = np.array(gt_binary_obj, dtype = np.bool8).copy()
        del gt_binary_obj
        del _
        print(np.shape(gt_binary))

        test_binary_obj, _, test_binary_hdr = get_image_objects(resampled_img_model)
        test_binary = np.array(test_binary_obj, dtype = np.bool8).copy()
        del test_binary_obj
        del _
        print(np.shape(test_binary))


        # define pixel dimensions for the distance calculations
        pix_dimensions = test_binary_hdr['pixdim']
        voxel_Spacing = [pix_dimensions[1], pix_dimensions[2], pix_dimensions[3]]


        # catch any images which now contain no image because of error 
        #Empty = no voxels 
        if (is_empty(test_binary)) or (is_empty(gt_binary)):
            # write to file 
            file.write(str(patient) + ' ' + str(model) + ' ' + str(timepoint) + ' ' + str(structure) + ' Empty Empty Empty \n')
        else:
            
            # calculate statistics 
            DICE_score = dc(test_binary, gt_binary)
            HausdorfDistance = hd95(test_binary, gt_binary, voxelspacing=voxel_Spacing, connectivity=1)
            AverageDistance = asd(test_binary, gt_binary, voxelspacing=voxel_Spacing, connectivity=1)

            # write to file 
            file.write(str(patient) + ' ' + str(model) + ' ' + str(timepoint) + ' ' + str(structure) + ' ' + str(DICE_score) + ' ' + str(AverageDistance) + ' ' + str(HausdorfDistance)  + '\n')

        


'''
        
        
'''
results = '/cluster/project7/HN_RT/CUH_HN/HN' + str(patient) +'_' + str(model) + '_geometric_results.txt'
if os.path.exists(results):
    os.remove(results)
file = open(results, 'a')
file.write('Patient Model Timepoint Structure DICE_Score AverageDistance HausdorfDistance \n')

'''


