import os
import numpy as np 
from datetime import datetime, date

path1 = 'D:/CUH_HN/NIFTI/'

differences = []


patients = [patient for patient in os.listdir(path1) if patient.startswith('CUH-UCL')]

for patient in patients: 

    
    CBCTs = [CBCT.replace('CBCT_', '') for CBCT in os.listdir(path1 + '/' + str(patient)) if CBCT.startswith('CBCT') and CBCT not in ['CBCT_pCT', 'CBCT_GROUPWISE']]
    
    
    dates = [CT[-8:] for CT in os.listdir(path1+ '/' + str(patient)) if CT.startswith('CT_')]
    '''CT_datetimes = [datetime(int(CT[0:4]), int(CT[4:6]), int(CT[6:8])) for CT in dates]
    pCT_datetime = min(CT_datetimes)
    
    
    CBCT = CBCTs[0]
    CBCT_datetime = datetime(int(CBCT[0:4]), int(CBCT[4:6]), int(CBCT[6:8]))
    difference = CBCT_datetime - pCT_datetime
    print(difference)
    '''
    
    if len(dates) ==2:
        print(patient)


#print(min(differences))
#print(max(differences))