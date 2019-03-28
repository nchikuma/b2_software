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

RUNID=$1
SRUNID=`expr $2 "*" $NUMJOB`

for i in `seq 1 $NUMJOB`
do 
  eval "timeout 5m ./mc_muontest.sh $RUNID $SRUNID"
  #echo "timeout 5m ./mc_muontest.sh $RUNID $SRUNID"
  SRUNID=`expr $SRUNID + 1`
done
