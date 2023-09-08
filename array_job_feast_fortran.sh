#!/bin/bash

#SBATCH -N1             # jeden
#SBATCH -c24            # 50 rdzeni per proces
#SBATCH --mem=192gb	   # 10 RAM
#SBATCH --time=24:00:00   # limit czasu ustawiony na 1h
#SBATCH --mail-user=jakub.jastrzebski99@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --array=10-24:2 # lista ID podzadań

#SBATCH --export=none
#SBATCH --job-name=${$2}


module load intel/2022a

./heisenberg_main $SLURM_ARRAY_TASK_ID 1.0 > log_feast_16.out



