if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /gpfs/fs03/t2k/beam/work/nchikuma/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh
set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=INGRID -version=v1r1 -path=/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/INGRID $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

