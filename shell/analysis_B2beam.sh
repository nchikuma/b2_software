#!/bin/sh

SPILLMERGE=0
BSDMERGE=0
TWODRECON=0
THRDRECON=0
PLOT_ALL=0
DETRES=0
THRESHOLD=1
TRKEFF=0
PECHECK=0
HITEFF=0


ONLYNSEL=1

OPT_ONLYNSEL=""
if test $ONLYNSEL -eq 1
then
  OPT_ONLYNSEL="-p"
fi


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

expr "${RUNID}"  + 1 > /dev/null 2>&1;status1=$?
expr "${SRUNID}" + 1 > /dev/null 2>&1;status2=$?
if [ $status1 -ne 0 -o $status2 -ne 0 ] ; then
  echo "run_id must be an integer."
  exit
fi


WORK_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software"
DATA_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2/data"
WAGASCI_SOFT_DIR="$WORK_DIR/wagasci_software/Analysis/bin"
INGRID_SOFT_DIR="$WORK_DIR/INGRID/INGRID/v1r1"  
INGRID_FORM_DIR="$WORK_DIR/INGRID/ingrid_format"
ANALYSIS_DIR="$WORK_DIR/analysis/app"
CALIB_FILE=`printf      "/home/t2k/nchikuma/b2_data/data_dst/ingrid_%08d_%04d_bsd.root"     ${RUNID} ${SRUNID}`
TWODRECON_FILE=`printf  "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_recon.root"  ${RUNID} ${SRUNID}`
THRDRECON_FILE=`printf  "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_anas1.root"  ${RUNID} ${SRUNID}`
BSD_MERGED_FILE=`printf "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_bsd1.root"   ${RUNID} ${SRUNID}`
PLOT_FILE=`printf       "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_plot.root"   ${RUNID} ${SRUNID}`
TRKEFF_FILE=`printf     "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_trkeff.root" ${RUNID} ${SRUNID}`
PE_FILE=`printf         "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_pe.root"     ${RUNID} ${SRUNID}`
HITEFF_FILE=`printf     "/home/t2k/nchikuma/b2_data2/data_dst/ingrid_%08d_%04d_hiteff.root" ${RUNID} ${SRUNID}`
DETRES_DIR="/home/t2k/nchikuma/b2_data3/sysE_detRes/data"
THRES_DIR="/home/t2k/nchikuma/b2_data2/pe_thres/data"


##Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd

#Get corresponding WAGASCI rundid
#RUNLIST="${WORK_DIR}/shell/wagasci_runid.txt"
RUNLIST2="${WORK_DIR}/shell/wagasci_runid2.txt"
if [ ! -f $RUNLIST2 ];then
  echo "No such files : $RUNLIST2"
  exit
fi

cmd=`printf "cat %s | grep '%08d %04d' | awk '{print $%d}'" ${RUNLIST2} ${RUNID} ${SRUNID} 3`
wgrun_i=`eval $cmd`
cmd=`printf "cat %s | grep '%08d %04d' | awk '{print $%d}'" ${RUNLIST2} ${RUNID} ${SRUNID} 4`
wgrun_f=`eval $cmd`

#cmd=`printf "cat %s | grep '%08d %04d' | awk '{print $%d}'" ${RUNLIST} ${RUNID} ${SRUNID} 3`
#num_wgacq=`eval $cmd`
#cmd=`printf "cat %s | grep '%08d %04d' | awk '{print $%d}'" ${RUNLIST} ${RUNID} ${SRUNID} 5`
#wgacq_i=`eval $cmd`

echo "WAGASCI Run ID : $wgrun_i - $wgrun_f"


if test $SPILLMERGE -eq 1
then
  ##Merge WAGASCI data
  cmd=`printf "%s/wgMerge_spill -r %d -s %d -i %d -j %d" ${WAGASCI_SOFT_DIR} ${RUNID} ${SRUNID} ${wgrun_i} ${wgrun_f}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi

if test $BSDMERGE -eq 1
then
  cmd=`printf "%s/app/IngAddBSD -r %d -s %d -m 1 -v" ${INGRID_FORM_DIR} ${RUNID} ${SRUNID}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi

if test $TWODRECON -eq 1
then
  cmd=`printf "%s/TwoDimRecon -r %d -s %d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi

if test $THRDRECON -eq 1
then
  #Three Dim Recon
  #cmd=`printf "%s/ThreeDimRecon -r %d -s %d -m 0 -x" ${ANALYSIS_DIR} ${RUNID} ${SRUNID}`
  #echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -r %d -s %d -m 1 -x" ${ANALYSIS_DIR} ${RUNID} ${SRUNID}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

  #BSD Merge
  #cmd=`printf "%s/app/IngAddBSD -r %d -s %d -m 0 -v" ${INGRID_FORM_DIR} ${RUNID} ${SRUNID}`
  #echo -e "---------\n ${cmd}\n"; eval $cmd
  #cmd=`printf "%s/app/IngAddBSD -r %d -s %d -m 1 -v" ${INGRID_FORM_DIR} ${RUNID} ${SRUNID}`
  #echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $PLOT_ALL -eq 1
then
  ##Plot 
  #cmd=`printf "%s/B2Plot -i %s -o %s -m 0" ${ANALYSIS_DIR} ${BSD_MERGED_FILE} ${PLOT_FILE}`
  cmd=`printf "%s/B2Plot -i %s -o %s -m 0 %s" ${ANALYSIS_DIR} ${THRDRECON_FILE} ${PLOT_FILE} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $TRKEFF -eq 1
then
  ##Plot 
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 0" ${ANALYSIS_DIR} ${THRDRECON_FILE} ${TRKEFF_FILE}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi

if test $PECHECK -eq 1
then
  #P.E
  cmd=`printf "%s/B2pe -i %s -o %s -m 0" ${ANALYSIS_DIR} ${THRDRECON_FILE} ${PE_FILE}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $HITEFF -eq 1
then
  #P.E
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 0" ${ANALYSIS_DIR} ${THRDRECON_FILE} ${HITEFF_FILE}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $DETRES -eq 1
then
  echo "Plot for sys. err. from detectro response"

  for dial in `seq 0 39`
  do
    DET_ANAFILE=`printf  "%s/dial_%02d/data_dst/ingrid_%08d_%04d_anas1.root" $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE=`printf "%s/dial_%02d/data_dst/ingrid_%08d_%04d_plot.root"  $DETRES_DIR $dial $RUNID $SRUNID`
    if test $dial -lt 27
    then
      #Three Dim Recon
      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -x -n %d" ${ANALYSIS_DIR} ${TWODRECON_FILE} ${DET_ANAFILE} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      #Plot 
      cmd=`printf "%s/B2Plot -i %s -o %s -m 0 -p" ${ANALYSIS_DIR} ${DET_ANAFILE} ${DET_PLOTFILE}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
    else
      ##Plot 
      cmd=`printf "%s/B2Plot -i %s -o %s -m 0 -n %d -p" ${ANALYSIS_DIR} ${THRDRECON_FILE} ${DET_PLOTFILE} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
    fi
  done

fi


if test $THRESHOLD -eq 1
then
  echo "Plot for sys. err. from threshold pe"


  for dial in `seq 0 5`
  do

    THRES_RECONFILE1=`printf  "%s/dial_%d/data_dst/ingrid_%08d_%04d_recon.root" $THRES_DIR $dial $RUNID $SRUNID`
    THRES_ANAFILE1=`printf    "%s/dial_%d/data_dst/ingrid_%08d_%04d_anas1.root" $THRES_DIR $dial $RUNID $SRUNID`
    THRES_PLOTFILE1=`printf   "%s/dial_%d/data_dst/ingrid_%08d_%04d_plot.root"  $THRES_DIR $dial $RUNID $SRUNID`

    #Recon
    cmd=`printf "%s/TwoDimRecon -f %s -o %s -t %d" ${ANALYSIS_DIR} ${CALIB_FILE} ${THRES_RECONFILE1} ${dial}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -x %s" ${ANALYSIS_DIR} ${THRES_RECONFILE1} ${THRES_ANAFILE1} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_RECONFILE1}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 0 %s" ${ANALYSIS_DIR} ${THRES_ANAFILE1} ${THRES_PLOTFILE1} ${OPT_ONLYNSEL}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_ANAFILE1}`
    echo -e "---------\n ${cmd}\n"; eval $cmd

  done

fi
