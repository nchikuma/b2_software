----------> uses
# use nd280Policy * 
# use ROOT * 
#   use EXTERN * 
#     use nd280Policy * 
#   use MYSQL * 
#     use EXTERN * 
# use oaEvent * 
#   use nd280Policy * 
#   use ROOT * 
#   use testBase * 
#     use nd280Policy * 
# use oaRawEvent * 
#   use nd280Policy * 
#   use ROOT * 
#   use oaEvent * 
#
# Selection :
use CMT v1r20p20081118 (/gpfs/fs03/t2k/beam/work/nchikuma)
use nd280Policy v2r21  (/gpfs/fs03/t2k/beam/work/nchikuma)
use testBase v1r5  (/gpfs/fs03/t2k/beam/work/nchikuma)
use EXTERN v3r1  (/gpfs/fs03/t2k/beam/work/nchikuma)
use MYSQL v5r051a  (/gpfs/fs03/t2k/beam/work/nchikuma)
use ROOT v5r24p00n02  (/gpfs/fs03/t2k/beam/work/nchikuma)
use oaEvent v7r3  (/gpfs/fs03/t2k/beam/work/nchikuma)
use oaRawEvent v3r5  (/gpfs/fs03/t2k/beam/work/nchikuma)
----------> tags
CMTv1 (from CMTVERSION)
CMTr20 (from CMTVERSION)
CMTp20081118 (from CMTVERSION)
Linux (from uname) package CMT implies [Unix]
amd64_linux26 (from CMTCONFIG)
work_config (from PROJECT) excludes [work_no_config]
work_root (from PROJECT) excludes [work_no_root]
work_cleanup (from PROJECT) excludes [work_no_cleanup]
work_prototypes (from PROJECT) excludes [work_no_prototypes]
work_without_installarea (from PROJECT) excludes [work_with_installarea]
work_with_version_directory (from PROJECT) excludes [work_without_version_directory]
work (from PROJECT)
x86_64 (from package CMT)
slc67 (from package CMT)
gcc447 (from package CMT)
Unix (from package CMT) excludes [WIN32 Win32]
----------> CMTPATH
# Add path /gpfs/fs03/t2k/beam/work/nchikuma from initialization
