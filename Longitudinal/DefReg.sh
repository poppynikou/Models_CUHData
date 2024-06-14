# Scheduler directives
#$ -S /bin/bash
#$ -l h_rt=08:00:00
#$ -l tmem=6G
#$ -l h_vmem=6G
#$ -j y
#$ -cwd
#$ -N DefRegs_GSTT
#$ -pe smp 5
#$ -R y 
#$ -t 1-1211

#export lib path
export LD_LIBRARY_PATH=/share/apps/gcc-8.3/lib64:$LD_LIBRARY_PATH

#source file location
source /SAN/medic/RTIC-MotionModel/software/niftyReg/niftyReg.source

#f3d command
F3D_CMD="/SAN/medic/RTIC-MotionModel/software/niftyReg/install/bin/reg_f3d"

#root data folder on the cluster
#path_to_data="/home/pnikou/Documents"

# each parallel job loops over one line in each of these files
CBCT=$(head -${SGE_TASK_ID} CBCT_directories.txt | tail -1)
CT=$(head -${SGE_TASK_ID} CT_directories.txt | tail -1)

# prints the paths
echo ${CBCT}
echo ${CT}

# creates a folder to store the outputs
# only if this path doesnt already exist
#path_to_out=$(head -${SGE_TASK_ID} resampled_imgs.txt | tail -1)
#mkdir -p ${path_to_out}

# defines the names of the output files
CBCT_DEF=$(head -${SGE_TASK_ID} DEF_CBCT_directories.txt | tail -1)
cpp_grid=$(head -${SGE_TASK_ID} cpp_CBCT_directories.txt | tail -1)

# parameters for registrations
parameters="-be 0 --lncc -5 -ln 5 -vel -le 0.01 -sx -10 -sy -10 -sz -10 -omp 10 -pad nan" 

if [ -f cpp_grid ]; then
    echo 'File exists.'
else
    #submit registration
    ${F3D_CMD} ${parameters} -flo ${CT} -ref ${CBCT} -res ${CBCT_DEF} -cpp ${cpp_grid}
fi



