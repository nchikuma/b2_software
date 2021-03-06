#-- start of make_header -----------------

#====================================
#  Document version
#
#   Generated Tue May 22 23:29:30 2018  by nchikuma
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_version_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_version_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_version

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_version = /tmp/CMT_$(INGRID_tag)_version.make$(cmt_lock_pid)
else
#cmt_local_tagfile_version = $(INGRID_tag)_version.make
cmt_local_tagfile_version = $(bin)$(INGRID_tag)_version.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_version = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_version = $(INGRID_tag).make
cmt_local_tagfile_version = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_version)

ifdef cmt_version_has_target_tag

ifdef READONLY
cmt_final_setup_version = /tmp/CMT_INGRID_versionsetup.make
cmt_local_version_makefile = /tmp/CMT_version$(cmt_lock_pid).make
else
cmt_final_setup_version = $(bin)INGRID_versionsetup.make
cmt_local_version_makefile = $(bin)version.make
endif

else

ifdef READONLY
cmt_final_setup_version = /tmp/CMT_INGRIDsetup.make
cmt_local_version_makefile = /tmp/CMT_version$(cmt_lock_pid).make
else
cmt_final_setup_version = $(bin)INGRIDsetup.make
cmt_local_version_makefile = $(bin)version.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

version ::


ifdef READONLY
version ::
	@echo tags=$(tags)
	@echo cmt_local_tagfile=$(cmt_local_tagfile)
endif


dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'version'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = version/
version::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

# -*- makefile -*-
# 
# A fragment used by CMT to build the library version file
#

version_output = $(src)

version :: 
	@echo "------> Updated $(src)$(package)_version.h"

ana_MPPC_cxx_dependencies = ../src/ana_MPPC.cxx
INGRID_Dimension_cxx_dependencies = ../src/INGRID_Dimension.cxx
WM_Ch_config_cxx_dependencies = ../src/WM_Ch_config.cxx
B2ING_Ch_config_cxx_dependencies = ../src/B2ING_Ch_config.cxx
PM_Ch_config_cxx_dependencies = ../src/PM_Ch_config.cxx
INGRID_Ch_config_cxx_dependencies = ../src/INGRID_Ch_config.cxx
LI_Ch_config_cxx_dependencies = ../src/LI_Ch_config.cxx
TINGRID_version_cxx_dependencies = ../src/TINGRID_version.cxx
INGRID_BadCh_mapping_cxx_dependencies = ../src/INGRID_BadCh_mapping.cxx
TINGRID_version_Dict_cxx_dependencies = ../dict/TINGRID_version_Dict.cxx
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)ana_MPPC.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)ana_MPPC.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)INGRID_Dimension.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)INGRID_Dimension.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)WM_Ch_config.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)WM_Ch_config.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)B2ING_Ch_config.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)B2ING_Ch_config.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)PM_Ch_config.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)PM_Ch_config.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)INGRID_Ch_config.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)INGRID_Ch_config.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)LI_Ch_config.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)LI_Ch_config.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)TINGRID_version.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)TINGRID_version.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: $(src)INGRID_BadCh_mapping.cxx
	@touch $(src)$(package)_version.h
	@echo $(src)INGRID_BadCh_mapping.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
#
# A fragment used by CMT to build the version include file.
#

version :: $(src)$(package)_version.h

$(src)$(package)_version.h :: ../dict/TINGRID_version_Dict.cxx
	@touch $(src)$(package)_version.h
	@echo ../dict/TINGRID_version_Dict.cxx
	@echo '/* Version information autogenerated by version fragment */' > $(src)$(package)_version.h
	@echo \#define $(package)_NAME \"$(package)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_VERSION \"$(version)\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DATE \"`date -u`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_HOST \"`uname -n`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_UNAME \"`uname -a`\" >> $(src)$(package)_version.h
	@echo \#define $(package)_COMPILE_DIR \"`pwd`\" >> $(src)$(package)_version.h
# -*- makefile -*-
# 
# A trailer fragment used by CMT to build the version include file.
#

clean :: versionclean

versionclean :: 
	/bin/rm $(src)$(package)_version.h || true

#-- start of cleanup_header --------------

clean :: versionclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(version.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(version.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_version)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(version.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(version.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(version.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

versionclean ::
#-- end of cleanup_header ---------------
