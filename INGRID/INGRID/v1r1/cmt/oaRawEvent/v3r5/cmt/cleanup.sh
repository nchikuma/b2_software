if test "${CMTROOT}" = ""; then
  CMTROOT=/gpfs/fs03/t2k/beam/work/nchikuma/CMT/v1r20p20081118; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=oaRawEvent -version=v3r5 -path=/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/INGRID/INGRID/v1r1/cmt $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

