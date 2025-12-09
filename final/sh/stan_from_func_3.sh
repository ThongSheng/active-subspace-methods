#!/bin/sh
 # computing info
#SBATCH --job-name=as
#SBATCH --mail-user=angt@purdue.edu
#SBATCH --mail-type=END
# fives cpus per task gives 5 cores; stan can run 1 chain per core in parallel
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000m
#SBATCH --time=16:00:00
#SBATCH --account=statdept
#SBATCH --output=/home/%u/logs/%x-%j.log
#SBATCH --array=000-999
#SBATCH --partition=cpu
#SBATCH --qos=standby

VALUES=({2001..3000}) # negishi has limit that job array can only have 1000 jobs at a time
    # if we want to run more than 1000, can replace this to 1001-2000, ..., then resubmit with --array=000-999
export THISJOBVALUE=${VALUES[$SLURM_ARRAY_TASK_ID]} # following https://github.com/HALFpipe/HALFpipe/discussions/438
export out_file=/home/angt/output/stan_$THISJOBVALUE.out # where the R output goes
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export BLAS_NUM_THREADS=2
# The application(s) to execute along with its input arguments and options:
cd /home/angt
module load r/4.4.1
R --no-save --no-restore CMD BATCH final/stan_from_func.R $out_file
