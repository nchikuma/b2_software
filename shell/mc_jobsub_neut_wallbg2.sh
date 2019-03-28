#!/bin/sh


NUMJOB=10


if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "======================================="
  echo "Usage:                                 "
  echo "  ./mc_jobsub.sh <run_id> <srun_id/$NUMJOB (start)>"
  echo " NumJob = $NUMJOB"
  echo "======================================="
  exit
fi

RUNID=`expr $1 "*" $NUMJOB`
SRUNID=$2

for i in `seq 1 $NUMJOB`
do 
  eval "timeout 5m ./mc_neut_wallbg2.sh $RUNID $SRUNID"
  #echo "timeout 3m ./mc_neut_wallbg.sh $RUNID $SRUNID"
  SRUNID=`expr $SRUNID + 1`
done
