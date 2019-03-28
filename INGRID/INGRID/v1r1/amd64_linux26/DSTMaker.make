#-- start of make_header -----------------

#====================================
#  Application DSTMaker
#
#   Generated Tue May 22 23:29:31 2018  by nchikuma
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_DSTMaker_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_DSTMaker_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_DSTMaker

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_DSTMaker = /tmp/CMT_$(INGRID_tag)_DSTMaker.make$(cmt_lock_pid)
else
#cmt_local_tagfile_DSTMaker = $(INGRID_tag)_DSTMaker.make
cmt_local_tagfile_DSTMaker = $(bin)$(INGRID_tag)_DSTMaker.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_DSTMaker = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_DSTMaker = $(INGRID_tag).make
cmt_local_tagfile_DSTMaker = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_DSTMaker)

ifdef cmt_DSTMaker_has_target_tag

ifdef READONLY
cmt_final_setup_DSTMaker = /tmp/CMT_INGRID_DSTMakersetup.make
cmt_local_DSTMaker_makefile = /tmp/CMT_DSTMaker$(cmt_lock_pid).make
else
cmt_final_setup_DSTMaker = $(bin)INGRID_DSTMakersetup.make
cmt_local_DSTMaker_makefile = $(bin)DSTMaker.make
endif

else

ifdef READONLY
cmt_final_setup_DSTMaker = /tmp/CMT_INGRIDsetup.make
cmt_local_DSTMaker_makefile = /tmp/CMT_DSTMaker$(cmt_lock_pid).make
else
cmt_final_setup_DSTMaker = $(bin)INGRIDsetup.make
cmt_local_DSTMaker_makefile = $(bin)DSTMaker.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

DSTMaker ::


ifdef READONLY
DSTMaker ::
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
	$(echo) 'DSTMaker'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = DSTMaker/
DSTMaker::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of application_header

DSTMaker :: dirs  $(bin)DSTMaker${application_suffix}
	$(echo) "DSTMaker ok"

#-- end of application_header
#-- start of application

$(bin)DSTMaker${application_suffix} :: $(bin)DSTMaker.o $(use_stamps) $(DSTMakerstamps) requirements $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)DSTMaker.o $(cmt_installarea_linkopts) $(DSTMaker_use_linkopts) $(DSTMakerlinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
DSTMakerinstallname = DSTMaker${application_suffix}

DSTMaker :: DSTMakerinstall

install :: DSTMakerinstall

DSTMakerinstall :: $(install_dir)/$(DSTMakerinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(DSTMakerinstallname) :: $(bin)$(DSTMakerinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(DSTMakerinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##DSTMakerclean :: DSTMakeruninstall

uninstall :: DSTMakeruninstall

DSTMakeruninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(DSTMakerinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (DSTMaker.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),DSTMakerclean)

#$(bin)DSTMaker_dependencies.make :: dirs

ifndef QUICK
$(bin)DSTMaker_dependencies.make : ../app/DSTMaker.cxx $(use_requirements) $(cmt_final_setup_DSTMaker)
	$(echo) "(DSTMaker.make) Rebuilding $@"; \
	  $(build_dependencies) DSTMaker -all_sources -out=$@ ../app/DSTMaker.cxx
endif

#$(DSTMaker_dependencies)

-include $(bin)DSTMaker_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp ------

$(bin)DSTMaker_dependencies.make : $(DSTMaker_cxx_dependencies)

$(bin)$(binobj)DSTMaker.o : $(DSTMaker_cxx_dependencies)
	$(cpp_echo) ../app/DSTMaker.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(DSTMaker_pp_cppflags) $(app_DSTMaker_pp_cppflags) $(DSTMaker_pp_cppflags) $(use_cppflags) $(DSTMaker_cppflags) $(app_DSTMaker_cppflags) $(DSTMaker_cppflags) $(DSTMaker_cxx_cppflags) -I../app ../app/DSTMaker.cxx

#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: DSTMakerclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(DSTMaker.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(DSTMaker.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_DSTMaker)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(DSTMaker.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(DSTMaker.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(DSTMaker.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

DSTMakerclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) DSTMaker${application_suffix}
	-$(cleanup_silent) cd $(bin); /bin/rm -f DSTMaker${application_suffix}

#	@echo "------> (DSTMaker.make) Removing application files"
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects
	-$(cleanup_silent) /bin/rm -f $(bin)DSTMaker.o
	-$(cleanup_silent) cd $(bin); /bin/rm -rf DSTMaker_deps DSTMaker_dependencies.make
#-- end of cleanup_objects ------
