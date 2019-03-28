#!/bin/sh


NEUT_HEAD="/gpfs/fs03/t2k/beam/work/nchikuma/data/neutfile/nd2_7_8/nd7/nega250/13a_nd7_nega250"
for i in `seq 0 1200`;
do 
  NEUTFILE1=`printf "%s_numubar_h2o_%d.nt" $NEUT_HEAD $i`
  NEUTFILE2=`printf "%s_numu_h2o_%d.nt" $NEUT_HEAD $i`

  if [ -f $NEUTFILE1 -a -f $NEUTFILE2 ];
  then
    printf "0 %04d\n" $i
  fi
done

