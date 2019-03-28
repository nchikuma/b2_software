# echo "Setting ROOT v5r24p00n02 in /gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/INGRID/INGRID/v1r1/cmt"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /gpfs/fs03/t2k/beam/work/nchikuma/CMT/v1r20p20081118
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=ROOT -version=v5r24p00n02 -path=/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/INGRID/INGRID/v1r1/cmt  -no_cleanup $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

