# Scheduler directives
#$ -S /bin/bash
#$ -l h_rt=12:00:00
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -j y
#$ -cwd
#$ -N CUH_Transport 
#$ -t 1-4

#export lib path
export LD_LIBRARY_PATH=/share/apps/gcc-8.3/lib64:$LD_LIBRARY_PATH

#path to niftireg executables 
export PATH=/SAN/medic/RTIC-MotionModel/software/niftyReg/install/bin:${PATH}
export LD_LIBRARY_PATH=/SAN/medic/RTIC-MotionModel/software/niftyReg/install/bin:${LD_LIBRARY_PATH}

# path to executable of python installed in my environment 
python_poppy="/home/pnikou/.conda/envs/hn_atlas/bin/python3.9"

patient=$(head -${SGE_TASK_ID} Patients.txt | tail -1)

model_folder='PSM_noLOO_T_All_CPS_4'

#path to save log file
path_to_log_file="/home/pnikou/Documents/CUH_Code/Results"

# command to run the code 
$python_poppy Transport.py  $patient $model_folder -u $path_to_log_file