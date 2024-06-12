# Scheduler directives
#$ -S /bin/bash
#$ -l h_rt=08:00:00
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -j y
#$ -cwd
#$ -N Atlas_CUH
#$ -R y 
#$ -t 1-45
#$ -pe smp 5

#export lib path
export LD_LIBRARY_PATH=/share/apps/gcc-8.3/lib64:$LD_LIBRARY_PATH

#path to niftireg executables 
export PATH=/SAN/medic/RTIC-MotionModel/software/niftyReg/install/bin:${PATH}
export LD_LIBRARY_PATH=/SAN/medic/RTIC-MotionModel/software/niftyReg/install/bin:${LD_LIBRARY_PATH}

# path to executable of python installed in my environment 
python_poppy="/home/pnikou/.conda/envs/hn_atlas/bin/python3.9"

#root data folder on the cluster
path_to_data="/cluster/project7/HN_RT/CUH_HN/"

patient=$(head -${SGE_TASK_ID} Patients.txt | tail -1)

#path to save log file
path_to_log_file="/home/pnikou/"

# command to run the code 
$python_poppy Rigid.py $path_to_data  $patient -u $path_to_log_file