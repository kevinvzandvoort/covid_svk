#!/bin/bash
#SBATCH --job-name=covid_svk --array=1-200

/nfs/general/covid_svk/single_run.sh $SLURM_ARRAY_TASK_ID
