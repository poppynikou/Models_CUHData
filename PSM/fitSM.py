import numpy as np 
from classes_CUH import *
import sys
import os 

'''
Add in the different options for the training time points
'random' 
'weekly' 
'specific time points' 
'''

# base path to where the batches of patient images are stored
base_path = sys.argv[1]

#patient 
patient = str(sys.argv[2])

print(patient)

# number of control points 
numcp = int(sys.argv[3])

results_path = "/cluster/project7/HN_RT/CUH_HN/PSM_noLOO_T_All_CPS_" + str(numcp) + '_ModelSpace'
if not os.path.exists(results_path):
        os.mkdir(results_path)

PSM_Model = PSM(base_path, patient, numcp, '', results_path)
PSM_Model.set_training_time_points()
PSM_Model.set_testing_time_points()
PSM_Model.set_reference_data()

PSM_Model.fit_SM()
PSM_Model.test_SM()
PSM_Model.save_SM() 

results_path = "/cluster/project7/HN_RT/CUH_HN/PSM_noLOO_T_35_CPS_" + str(numcp)+ '_ModelSpace'
if not os.path.exists(results_path):
        os.mkdir(results_path)

PSM_Model = PSM(base_path, patient, numcp, '', results_path)   
training_time_points = np.arange(0,35)            
PSM_Model.set_training_time_points(training_time_points = training_time_points)
PSM_Model.set_testing_time_points(testing_time_points=training_time_points)
PSM_Model.set_reference_data()

PSM_Model.fit_derivative_SM()
PSM_Model.save_temporal_CPG()
PSM_Model.test_SM()
PSM_Model.save_SM()

