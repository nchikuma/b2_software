#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "==============================================="
  echo "Usage:                                         "
  echo "  ./mc_cosmic.sh <run_id> <srun_id>"
  echo "==============================================="
  exit
fi

RUNID=$1;
SRUNID=$2;

expr "${RUNID}"  + 1 > /dev/null 2>&1;status1=$?
expr "${SRUNID}" + 1 > /dev/null 2>&1;status2=$?
if [ $status1 -ne 0 -o $status2 -ne 0 ] ; then
  echo "run_id must be an integer."
  exit
fi


WORK_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software"
WAGASCI_SOFT_DIR="$WORK_DIR/wagasci_software/Analysis/bin"
INGRID_SOFT_DIR="$WORK_DIR/INGRID/INGRID/v1r1"  
INGRID_FORM_DIR="$WORK_DIR/INGRID/ingrid_format"
ANALYSIS_DIR="$WORK_DIR/analysis/app"
MC_DIR="$WORK_DIR/b2mc"
NEUT_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/data/neutfile/nd2_7_8/nd7/nega250"
OUTPUT_SANDMU_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_sandmu"
OUTPUT_COSMIC_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_cosmic"
OUTPUTFILE=`printf "%s/ingrid_%08d_%04d_calib.root" $OUTPUT_COSMIC_DIR $RUNID $SRUNID`
RECONFILE=`printf "%s/ingrid_%08d_%04d_recon.root" $OUTPUT_COSMIC_DIR $RUNID $SRUNID`

##Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${MC_DIR}/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


#MC Job
#cosmic muon
cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f12 -o %s -e 1000" $MC_DIR $OUTPUTFILE`
echo -e "---------\n ${cmd}\n"; eval $cmd

#Recon
cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE} ${RECONFILE}`
echo -e "---------\n ${cmd}\n"; eval $cmd

#Hit efficiency analysis (except for cosmic muon with B2 WAGASCI)
for i in `seq 0 7`; do
  thres_pe=$i
  cmd=`printf "%s/HitEfficiency -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
done

#PE
for i in `seq 4 4`; do
  thres_pe=$i
  cmd=`printf "%s/B2pe -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
done

#Hit efficiency in Scinti
for i in `seq 4 4`; do
  thres_pe=$i
  cmd=`printf "%s/HitEfficiency_2 -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
done

#PE in scinti
thres_pe=0
cmd=`printf "%s/B2pe_Scinti -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
echo -e "---------\n ${cmd}\n"; eval $cmd
