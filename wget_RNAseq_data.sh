#!/bin/bash

# options for sbatch
#SBATCH --job-name=Abbasi # Job name
#SBATCH --nodes=1 # should never be anything other than 1
#SBATCH --cpus-per-task=16 # number of cpus to use
#SBATCH --mem=2000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=biopatic # cluster partition
##SBATCH --account=workshop # cluster account to use for the job
##SBATCH --array=1-16 # Task array indexing, see https://slurm.schedmd.com/job_array.html, the double # means this line is commented out
#SBATCH --output=stdout.out # File to which STDOUT will be written
#SBATCH --error=stderr.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nasir.abbasi@univ-tours.fr



wget https://objectstorage.uk-london-1.oraclecloud.com/p/-vsBVa5It4JipgFbath3ERXpmTIJ_FIGIt9CUEpy2Tg5P8usGf8wo41Hu-7faQf9/n/cnyr09qj8zbo/b/england-data/o/out/CP2022053000078/X204SC22123284-Z01-F001/X204SC22123284-Z01-F001.tar
