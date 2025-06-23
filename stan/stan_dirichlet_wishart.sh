#!/bin/sh
 # computing info
#SBATCH --job-name=as
#SBATCH --mail-user=anyarger@purdue.edu
#SBATCH --mail-type=END
# fives cpus per task gives 5 cores; stan can run 1 chain per core in parallel
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1200m
#SBATCH --time=03:59:00
#SBATCH --account=standby
#SBATCH --output=/home/%u/shit/%x-%j.log
#SBATCH --array=000-799

VALUES=({1..1000}) # negishi has limit that job array can only have 1000 jobs at a time
    # if we want to run more than 1000, can replace this to 1001-2000, ..., then resubmit with --array=000-999
export THISJOBVALUE=${VALUES[$SLURM_ARRAY_TASK_ID]} # following https://github.com/HALFpipe/HALFpipe/discussions/438
export out_file=/home/anyarger/shit/stan_$THISJOBVALUE.out # where the R output goes
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export BLAS_NUM_THREADS=2
# The application(s) to execute along with its input arguments and options:
cd /home/anyarger/active-subspace-methods/
module load r/4.4.1
R --no-save --no-restore CMD BATCH stan/stan_dirichlet_wishart.R $out_file
