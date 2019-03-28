#!/bin/sh

echo "Running... setting_wagasci.sh"

export WAGASCI_RAWDATADIR=/hsm/nu/wagasci/daqdata
export WAGASCI_DECODEDIR=/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/decode
#export WAGASCI_DECODEDIR=/home/t2k/nchikuma/b2_data2/decode
#export WAGASCI_LOGDIR=/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/log
export WAGASCI_LOGDIR=/home/t2k/nchikuma/b2_data2/log
export WAGASCI_CALIBDATADIR=/hsm/nu/wagasci/calibration
#export WAGASCI_RECONDIR=/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/recon
#export WAGASCI_SPILLDIR=/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/spill


###################

echo "WAGASCI_RAWDATADIR  =${WAGASCI_RAWDATADIR}  "
echo "WAGASCI_DECODEDIR   =${WAGASCI_DECODEDIR}   "
echo "WAGASCI_LOGDIR      =${WAGASCI_LOGDIR}      "
echo "WAGASCI_CALIBDATADIR=${WAGASCI_CALIBDATADIR}"
#echo "WAGASCI_RECONDIR    =${WAGASCI_RECONDIR}    "
#echo "WAGASCI_SPILLDIR    =${WAGASCI_SPILLDIR}    "

if [ ! -d ${WAGASCI_RAWDATADIR}   ]; then echo "No such a directory : ${WAGASCI_RAWDATADIR}    "; fi
if [ ! -d ${WAGASCI_DECODEDIR}    ]; then echo "No such a directory : ${WAGASCI_DECODEDIR}     "; fi
if [ ! -d ${WAGASCI_LOGDIR}       ]; then echo "No such a directory : ${WAGASCI_LOGDIR}        "; fi
if [ ! -d ${WAGASCI_CALIBDATADIR} ]; then echo "No such a directory : ${WAGASCI_CALIBDATADIR}  "; fi
#if [ ! -d ${WAGASCI_RECONDIR}     ]; then echo "No such a directory : ${WAGASCI_RECONDIR}      "; fi
#if [ ! -d ${WAGASCI_SPILLDIR}     ]; then echo "No such a directory : ${WAGASCI_SPILLDIR}      "; fi
