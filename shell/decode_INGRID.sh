#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "======================================="
  echo "Usage:                                 "
  echo "  ./decode_INGRID.sh <run_id> <srun_id>"
  echo "======================================="
  exit
fi

RUNID=$1;
SRUNID=$2;

#if [ $RUNID -eq 31097 -a $SRUNID -eq 1 ]; then exit; fi
#if [ $RUNID -eq 31119 -a $SRUNID -eq 1 ]; then exit; fi
#if [ $RUNID -eq 31148 -a $SRUNID -eq 1 ]; then exit; fi
#if [ $RUNID -eq 31158 -a $SRUNID -eq 1 ]; then exit; fi
#if [ $RUNID -eq 31167 -a $SRUNID -eq 2 ]; then exit; fi
#if [ $RUNID -eq 31169 -a $SRUNID -eq 3 ]; then exit; fi
#if [ $RUNID -eq 31186 -a $SRUNID -eq 1 ]; then exit; fi

echo "yes"

expr "${RUNID}"  + 1 > /dev/null 2>&1;status1=$?
expr "${SRUNID}" + 1 > /dev/null 2>&1;status2=$?
if [ $status1 -ne 0 -o $status2 -ne 0 ] ; then
  echo "run_id must be an integer."
  exit
fi


WORK_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2"
INGRID_SOFT_DIR="$WORK_DIR/b2_software/INGRID/INGRID/v1r1"  
INGRID_FORM_DIR="$WORK_DIR/b2_software/INGRID/ingrid_format"

#Set up
cmd="source ${WORK_DIR}/b2_software/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


#midas_file_kyoto=`printf "/export/scraid5/data1/raw_tmp/00030000_00039999/ingrid_%08d_%04d.daq.mid.gz" ${RUNID} ${SRUNID}`
#midas_file_local=`printf "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_midas/ingrid_%08d_%04d.daq.mid.gz" ${RUNID} ${SRUNID}`
#cmd="scp kyoto:${midas_file_kyoto} ${midas_file_local}"
#echo -e "---------\n ${cmd}\n"; eval $cmd


##Calibration
cmd=`printf "%s/amd64_linux26/Calc_MPPC_new.exe -r %d -s %d -t 1" ${INGRID_SOFT_DIR} ${RUNID} ${SRUNID}`
echo -e "---------\n ${cmd}\n"; eval $cmd
#
###Beam trigger
#cmd=`printf "%s/amd64_linux26/DSTMaker.exe -r %d -s %d -t 1" ${INGRID_SOFT_DIR} ${RUNID} ${SRUNID}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#cmd=`printf "%s/app/IngCalib_new -r %d -s %d" ${INGRID_FORM_DIR} ${RUNID} ${SRUNID}`
#echo -e "---------\n ${cmd}\n"; eval $cmd


##Cosmic trigger
#cmd=`printf "%s/amd64_linux26/DSTMaker.exe -r %d -s %d -t 128" ${INGRID_SOFT_DIR} ${RUNID} ${SRUNID}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#cmd=`printf "%s/app/IngCalib_new -r %d -s %d -c" ${INGRID_FORM_DIR} ${RUNID} ${SRUNID}`
#echo -e "---------\n ${cmd}\n"; eval $cmd



#
#CALIB_DATA=`printf "${WORK_DIR}/data/data_calib/ingrid_%08d_%04d_Calib00.root" ${RUNID} ${SRUNID}`
#DST_DATA=`printf "${WORK_DIR}/data/data_dst/ingrid_%08d_%04d.root" ${RUNID} ${SRUNID}`
#DST_CALIB_DATA=`printf "${WORK_DIR}/data/data_dst/ingrid_%08d_%04d_calib.root" ${RUNID} ${SRUNID}`
#COSMIC_DATA=`printf "${WORK_DIR}/data/data_cosmic/ingrid_%08d_%04d.root" ${RUNID} ${SRUNID}`
#COSMIC_CALIB_DATA=`printf "${WORK_DIR}/data/data_cosmic/ingrid_%08d_%04d_calib.root" ${RUNID} ${SRUNID}`
#
#
#sleep 1
#if [ -f $DST_CALIB_DATA -a -f $DST_DATA ]; then
#  echo "The following files are removed."
#  echo $DST_DATA
#  rm -f $DST_DATA
#fi
#if [ -f $COSMIC_CALIB_DATA -a -f $COSMIC_DATA ]; then
#  echo "The following files are removed."
#  echo $COSMIC_DATA
#  rm -f $COSMIC_DATA
#fi

#rm $midas_file_local
