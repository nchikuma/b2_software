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
OUTPUT_MUONTEST_DIR="/home/t2k/nchikuma/b2_data2/mc_muon_test"
OUTPUT_PITEST_DIR="/home/t2k/nchikuma/b2_data2/mc_pion_test"
OUTPUT_PROTEST_DIR="/home/t2k/nchikuma/b2_data2/mc_proton_test"
OUTPUTFILE1=`printf "%s/ingrid_%08d_%04d_calib.root"  $OUTPUT_MUONTEST_DIR $RUNID $SRUNID`
OUTPUTFILE2=`printf "%s/ingrid_%08d_%04d_calib.root"  $OUTPUT_PITEST_DIR   $RUNID $SRUNID`
OUTPUTFILE3=`printf "%s/ingrid_%08d_%04d_calib.root"  $OUTPUT_PROTEST_DIR  $RUNID $SRUNID`
RECONFILE1=`printf "%s/ingrid_%08d_%04d_recon.root"   $OUTPUT_MUONTEST_DIR $RUNID $SRUNID`
RECONFILE2=`printf "%s/ingrid_%08d_%04d_recon.root"   $OUTPUT_PITEST_DIR   $RUNID $SRUNID`
RECONFILE3=`printf "%s/ingrid_%08d_%04d_recon.root"   $OUTPUT_PROTEST_DIR  $RUNID $SRUNID`
ANAFILE1=`printf "%s/ingrid_%08d_%04d_anas1.root"     $OUTPUT_MUONTEST_DIR $RUNID $SRUNID`
ANAFILE2=`printf "%s/ingrid_%08d_%04d_anas1.root"     $OUTPUT_PITEST_DIR   $RUNID $SRUNID`
ANAFILE3=`printf "%s/ingrid_%08d_%04d_anas1.root"     $OUTPUT_PROTEST_DIR  $RUNID $SRUNID`
DETEFFFILE1=`printf "%s/ingrid_%08d_%04d_deteff.root" $OUTPUT_MUONTEST_DIR $RUNID $SRUNID`
DETEFFFILE2=`printf "%s/ingrid_%08d_%04d_deteff.root" $OUTPUT_PITEST_DIR   $RUNID $SRUNID`
DETEFFFILE3=`printf "%s/ingrid_%08d_%04d_deteff.root" $OUTPUT_PROTEST_DIR  $RUNID $SRUNID`

##Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${MC_DIR}/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


###MC Job
##muon injection
#cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f13 -o %s -e 1000" $MC_DIR $OUTPUTFILE1`
#echo -e "---------\n ${cmd}\n"; eval $cmd
##pion injection
#cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f14 -o %s -e 1000" $MC_DIR $OUTPUTFILE2`
#echo -e "---------\n ${cmd}\n"; eval $cmd
##proton injection
#cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f15 -o %s -e 1000" $MC_DIR $OUTPUTFILE3`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#
##2DRecon
#cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${RECONFILE1}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE2} ${RECONFILE2}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE3} ${RECONFILE3}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#
##3DRecon
#cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -q" ${ANALYSIS_DIR} ${RECONFILE1} ${ANAFILE1}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -q" ${ANALYSIS_DIR} ${RECONFILE2} ${ANAFILE2}`
#echo -e "---------\n ${cmd}\n"; eval $cmd
#cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -q" ${ANALYSIS_DIR} ${RECONFILE3} ${ANAFILE3}`
#echo -e "---------\n ${cmd}\n"; eval $cmd

##Detection Efficiency
cmd=`printf "%s/DetEff -i %s -o %s" ${ANALYSIS_DIR} ${ANAFILE1} ${DETEFFFILE1} `
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd=`printf "%s/DetEff -i %s -o %s" ${ANALYSIS_DIR} ${ANAFILE2} ${DETEFFFILE2} `
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd=`printf "%s/DetEff -i %s -o %s" ${ANALYSIS_DIR} ${ANAFILE3} ${DETEFFFILE3} `
echo -e "---------\n ${cmd}\n"; eval $cmd



##Reconstruction efficiency
#thres_nhit=3
#cmd=`printf "%s/TrackEfficiency_sim -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_nhit}`
#echo -e "---------\n ${cmd}\n"; eval $cmd


##Hit efficiency analysis (except for cosmic muon with B2 WAGASCI)
#for i in `seq 0 7`; do
#  thres_pe=$i
#  cmd=`printf "%s/HitEfficiency -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#done

#for i in `seq 5 5`; do
#  thres_pe=$i
#  cmd=`printf "%s/B2pe -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#done
