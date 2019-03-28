#!/bin/bash

WORK_DIR=/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software

#export CVSROOT=:ext:anoncvs@repo.nd280.org:/home/trt2kmgr/ND280Repository
#export CVSROOT=:ext:ingrid@scbn00.hepnet.scphys.kyoto-u.ac.jp:/export/scbn00/data1/INGRIDMC
export CVSROOT=:ext:ingrid@scbn00.hepnet.scphys.kyoto-u.ac.jp:/export/scbn00/data1/ingrid/mc
export CVS_RSH=ssh
unset CVS_SERVER
export CMTPATH=/gpfs/fs03/t2k/beam/work/nchikuma
source /gpfs/fs03/t2k/beam/work/nchikuma/CMT/v1r20p20081118/mgr/setup.sh

#CERNLIB
export CERN=/home/t2k/tatsuya1/cern/cernlib
export CERN_LEVEL=2006
export CERN_ROOT=${CERN}/${CERN_LEVEL}
export GCALOR=/home/t2k/nchikuma/gcalor/gcalor.o
export PATH=${CERN_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${CERN_ROOT}/lib:${LD_LIBRARY_PATH}

echo ""
echo "CERN LIB : ${CERN_ROOT}"
echo "GCALOR   : ${GCALOR}   "
if [ ! -f ${GCALOR}    ];then echo "ERROR: No such a file : ${GCALOR}";fi
if [ ! -d ${CERN_ROOT} ];then echo "ERROR: No such a directory :${CERN_ROOT}";fi

#Set Env various by cmt
BASESOFT_DIR1=/gpfs/fs03/t2k/beam/work/tatsuya1/INGRID/basesoft
BASESOFT_DIR2=/gpfs/fs03/t2k/beam/work/nchikuma
VERSION_ROOT=v5r24p00n02
VERSION_CLHEP=v2r0p3
VERSION_GEANT=v9r2p01n00
echo ""
export ND280ROOT=${BASESOFT_DIR1}/ROOT/${VERSION_ROOT}
export ND280CLHEP=${BASESOFT_DIR1}/CLHEP/${VERSION_CLHEP}
export ND280GEANT=${BASESOFT_DIR1}/GEANT/${VERSION_GEANT}
eval "source $ND280GEANT/cmt/setup.sh"
eval "source $ND280ROOT/cmt/setup.sh"
eval "source $ND280CLHEP/cmt/setup.sh"

echo "ROOT  Version : ${VERSION_ROOT} : PATH=${ND280ROOT}"
echo "CLHEP Version : ${VERSION_CLHEP}: PATH=${ND280CLHEP}"
echo "GEANT Version : ${VERSION_GEANT}: PATH=${ND280GEANT}"


echo ""
INGRID_SOFT=$WORK_DIR/INGRID/INGRID/v1r1
source $INGRID_SOFT/cmt/setup.sh
echo "INGRID software : PATH=${INGRID_SOFT}"

echo ""
if [ -f $WORK_DIR/setting_wagasci.sh ];then
  source $WORK_DIR/setting_wagasci.sh
else
  echo "ERROR: No setting_wagasci.sh"
fi
echo ""
