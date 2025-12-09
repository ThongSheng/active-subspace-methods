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
#SBATCH --array=000-119
#SBATCH --partition=cpu
#SBATCH --qos=normal

VALUES=({2975 2976 2978 2979 2984 2985 2987 2988 2993 2994 2996 2997 3002 3003 3005 3006 3011 3012 3014 3015 \
         3020 3021 3023 3024 3029 3030 3032 3033 3038 3039 3041 3042 3047 3048 3050 3051 3056 3057 3059 3060 \
         3065 3066 3068 3069 3074 3075 3077 3078 3083 3084 3086 3087 3092 3093 3095 3096 3101 3102 3104 3105 \
         3110 3111 3113 3114 3119 3120 3122 3123 3128 3129 3131 3132 3137 3138 3140 3141 3146 3147 3149 3150 \
         3155 3156 3158 3159 3164 3165 3167 3168 3173 3174 3176 3177 3182 3183 3185 3186 3191 3192 3194 3195 \
         3200 3201 3203 3204 3209 3210 3212 3213 3218 3219 3221 3222 3227 3228 3230 3231 3236 3237 3239 3240})
export THISJOBVALUE=${VALUES[$SLURM_ARRAY_TASK_ID]} # following https://github.com/HALFpipe/HALFpipe/discussions/438
export out_file=/home/angt/output/stan_$THISJOBVALUE.out # where the R output goes
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export BLAS_NUM_THREADS=2
# The application(s) to execute along with its input arguments and options:
cd /home/angt
module load r/4.4.1
R --no-save --no-restore CMD BATCH final/stan_from_func.R $out_file
