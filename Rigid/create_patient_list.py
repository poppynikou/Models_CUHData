import os 



path = 'D:/CUH_HN/NIFTI/'

patients = [folder for folder in os.listdir(path) if folder.__contains__('CUH-UCL')]



text_file = open('Patients.txt', 'a')

for patient in patients:

    text_file.write(str(patient) + '\n')