#!/bin/bash

#SBATCH -N1             # jeden
#SBATCH -c24            # 50 rdzeni per proces
#SBATCH --mem=100gb	   # 10 RAM
#SBATCH --time=2:00:00   # limit czasu ustawiony na 1h
#SBATCH --mail-user=jakub.jastrzebski99@gmail.com
#SBATCH --mail-type=ALL

# Dodatkowo a
#SBATCH --export=none
#SBATCH --job-name=${$2}
# Nascja.

module load intel/2022a
module load OpenMPI/4.1.4-GCC-11.3.0

./heisenberg_main 16 1.0 > log_feast_16.out

#for i in {4..14..2}
#do
 #   ./a.out $i 1.0 > log_$i.out
#done



