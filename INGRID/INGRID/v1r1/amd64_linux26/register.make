#-- start of make_header -----------------

#====================================
#  Document register
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

cmt_register_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_register_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_register

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_register = /tmp/CMT_$(INGRID_tag)_register.make$(cmt_lock_pid)
else
#cmt_local_tagfile_register = $(INGRID_tag)_register.make
cmt_local_tagfile_register = $(bin)$(INGRID_tag)_register.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_register = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_register = $(INGRID_tag).make
cmt_local_tagfile_register = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_register)

ifdef cmt_register_has_target_tag

ifdef READONLY
cmt_final_setup_register = /tmp/CMT_INGRID_registersetup.make
cmt_local_register_makefile = /tmp/CMT_register$(cmt_lock_pid).make
else
cmt_final_setup_register = $(bin)INGRID_registersetup.make
cmt_local_register_makefile = $(bin)register.make
endif

else

ifdef READONLY
cmt_final_setup_register = /tmp/CMT_INGRIDsetup.make
cmt_local_register_makefile = /tmp/CMT_register$(cmt_lock_pid).make
else
cmt_final_setup_register = $(bin)INGRIDsetup.make
cmt_local_register_makefile = $(bin)register.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

register ::


ifdef READONLY
register ::
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
	$(echo) 'register'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = register/
register::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

# -*- makefile -*-
# 
# A fragment used by CMT to register the library version with
# TOADatabase.
#

register_output = $(src)

register :: $(src)T$(package)_version.hxx $(src)T$(package)_version.cxx $(src)T$(package)_version_LinkDef.h
	@echo "------> Register $(package) version"

$(src)T$(package)_version.hxx :: $(OAEVENTROOT)/src/Tpackage_version.hxx.in
	sed s/%PACKAGE%/$(package)/g $^ > $@

$(src)T$(package)_version.cxx :: $(OAEVENTROOT)/src/Tpackage_version.cxx.in
	sed s/%PACKAGE%/$(package)/g $^ > $@

$(src)T$(package)_version_LinkDef.h :: $(OAEVENTROOT)/src/Tpackage_version_LinkDef.in
	sed s/%PACKAGE%/$(package)/g $^ > $@

# -*- makefile -*-
# 
# A fragment used by CMT to register the library with TOADatabase.
#

clean :: registerclean

registerclean :: 
	/bin/rm $(src)T$(package)_version.hxx || true
	/bin/rm $(src)T$(package)_version.cxx || true
	/bin/rm $(src)T$(package)_version_LinkDef.h || true

#-- start of cleanup_header --------------

clean :: registerclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(register.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(register.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_register)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(register.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(register.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(register.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

registerclean ::
#-- end of cleanup_header ---------------
