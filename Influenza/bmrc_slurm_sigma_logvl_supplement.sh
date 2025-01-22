#!/bin/bash

#SBATCH -A moru-batty.prj
#SBATCH -D /well/moru-batty/users/lcd199/Viral_clearance/Influenza
#SBATCH -J res
#SBATCH -n 4
#SBATCH -o /well/moru-batty/users/lcd199/Viral_clearance/Influenza/o_and_e_files/output.o%A_%a.out
#SBATCH -e /well/moru-batty/users/lcd199/Viral_clearance/Influenza/o_and_e_files/output.e%A_%a.out
#SBATCH -p short
#SBATCH --array 3501-3600


echo started=`date`
module purge
module load R/4.3.2-gfbf-2023a



echo "job=$SLURM_JOB_ID"
echo "hostname="`hostname`
echo "OS="`uname -s`
echo "username="`whoami`
Rscript /well/moru-batty/users/lcd199/Viral_clearance/Influenza/run_sample_size_calculations_sigma_logvl_supplement.R ${SLURM_ARRAY_TASK_ID} --no-save --no-restore
echo "finished="`date`
exit 0