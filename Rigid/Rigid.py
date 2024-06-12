import numpy as np 
from classes import *
import sys

# base path to where the batches of patient images are stored
base_path = sys.argv[1]
#patient 
patient =  str(sys.argv[2])
log_path = base_path + '/' + str(patient) + '.txt'

no_itterations = 4
atlas_path = '/home/pnikou/Documents/ATLAS_IMG/NIFTI/MASKED_average_pCT.nii.gz'

# if a certain path exists and is not empty 
# check before you do it again 
PatientData = PatientData(patient, base_path)
CBCT_dates = PatientData.get_dates()
CBCT_relative_timepoints = PatientData.get_CBCT_relative_timepoints()
pCT_date = PatientData.get_pCT_date()

'''
GroupwiseReg = GroupwiseRegs(patient, pCT_date, CBCT_dates, CBCT_relative_timepoints, no_itterations, base_path, log_path)

GroupwiseReg.refactor()

for itteration in np.arange(0, no_itterations):

    GroupwiseReg.set__itteration(itteration)
    GroupwiseReg.rigidGroupReg()
    GroupwiseReg.avgAffine()
    GroupwiseReg.invAffine()
    GroupwiseReg.compAffine()
    GroupwiseReg.resampleImages()
    GroupwiseReg.avgImage()


GroupwiseReg.UpdateGroupSform()

GroupwiseReg.rigidpCTReg()
GroupwiseReg.UpdateSform()

'''
AtlasAlignment = AtlasRegs(patient, pCT_date, base_path, atlas_path)

AtlasAlignment.refactor()

AtlasAlignment.InitAlignment()
AtlasAlignment.RigidReg()
AtlasAlignment.AffineReg()
AtlasAlignment.DefReg()

AtlasAlignment.Calc_Tatlas()

AtlasAlignment.ResampleImgs(CBCT_relative_timepoints)

