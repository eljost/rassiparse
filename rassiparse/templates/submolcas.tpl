#!/bin/bash

#SBATCH -J {{ job_name }}
#SBATCH -t 200:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=0
#SBATCH --hint=nomultithread
    
export CurrDir={{ curr_dir }}
export MOLCAS_PROJECT='{{ project }}'

module load molcas/8.2own

export MOLCAS_NPROCS=24
export MOLCAS_PRINT=2
export MOLCAS_MEM=4000
export MOLCAS_EXE=$MOLCAS/bin/molcas.exe
export MOLCAS_MOLDEN=ON
export MOLCAS_OUTPUT=$CurrDir
export MOLCAS_WORKDIR=/beegfs/si93yaq/mc_scratch/molcas${SLURM_JOB_ID}
 
$MOLCAS_EXE $CurrDir/{{ input_fn_no_ext }}.in &> $CurrDir/{{ input_fn_no_ext }}.out

rm -rf $MOLCAS_WORKDIR
