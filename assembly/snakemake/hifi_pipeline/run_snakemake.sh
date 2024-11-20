#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=100MB
#SBATCH --qos=standard
#SBATCH --time=2-00:00:00
#SBATCH --job-name=beetle_assembly
#SBATCH -o /scratch/galeov97/begendiv/beetle.out
#SBATCH -e /scratch/galeov97/begendiv/beetle.err
#SBATCH --mail-user=galeov97@zedat.fu-berlin.de
#SBATCH --mail-type=FAIL

snakemake --profile submit_config/slurm/
