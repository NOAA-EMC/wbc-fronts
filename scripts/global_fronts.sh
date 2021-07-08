#!/bin/bash -l

#-------------------------------------------------------
# WBC Frontal Analysis Validation System
# Version 1.0
# Todd Spindler
# 24 Aug 2020
#-------------------------------------------------------

# set up module environment
module use /scratch2/NCEPDEV/marine/Todd.Spindler/save/modulefiles
module load anaconda-work/1.0.0 mmab/1.0.0

THE_DATE=${1:-`date --date yesterday +%Y%m%d`}

NCORES=12

TASK_QUEUE='batch'
WALL='0:30:00'
PROJ='marine-cpu'
LOGPATH=/scratch2/NCEPDEV/stmp1/Todd.Spindler/logs/class-4/fronts
JOB="fronts"

mkdir -p $LOGPATH
rm -f $LOGPATH/*.log

export SRCDIR='/scratch2/NCEPDEV/marine/Todd.Spindler/save/VPPPG/Global_RTOFS/EMC_ocean-verification/fronts'

job1=$(sbatch --parsable -J ${JOB}_${THE_DATE} -o $LOGPATH/${JOB}_${THE_DATE}.log -q $TASK_QUEUE --account=$PROJ --time $WALL --ntasks=$(($NCORES + 1)) --nodes=1 --wrap "python $SRCDIR/ush/global_fronts.py $THE_DATE")

