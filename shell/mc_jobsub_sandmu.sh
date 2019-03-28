#!/bin/sh


NUMJOB=100


if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "======================================="
  echo "Usage:                                 "
  echo "  ./mc_sandmu.sh <run_id> <srun_id/$NUMJOB (start)>"
  echo " NumJob = $NUMJOB"
  echo "======================================="
  exit
fi

RUNID=$1
SRUNID=`expr $2 "*" $NUMJOB`

for i in `seq 1 $NUMJOB`
do 
  eval "timeout 5m ./mc_sandmu.sh $RUNID $SRUNID"
  SRUNID=`expr $SRUNID + 1`
done
