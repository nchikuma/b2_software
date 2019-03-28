#-- start of make_header -----------------

#====================================
#  Application Calc_MPPC_new
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

cmt_Calc_MPPC_new_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_Calc_MPPC_new_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_Calc_MPPC_new

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_Calc_MPPC_new = /tmp/CMT_$(INGRID_tag)_Calc_MPPC_new.make$(cmt_lock_pid)
else
#cmt_local_tagfile_Calc_MPPC_new = $(INGRID_tag)_Calc_MPPC_new.make
cmt_local_tagfile_Calc_MPPC_new = $(bin)$(INGRID_tag)_Calc_MPPC_new.make
endif

else

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile_Calc_MPPC_new = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile_Calc_MPPC_new = $(INGRID_tag).make
cmt_local_tagfile_Calc_MPPC_new = $(bin)$(INGRID_tag).make
endif

endif

-include $(cmt_local_tagfile_Calc_MPPC_new)

ifdef cmt_Calc_MPPC_new_has_target_tag

ifdef READONLY
cmt_final_setup_Calc_MPPC_new = /tmp/CMT_INGRID_Calc_MPPC_newsetup.make
cmt_local_Calc_MPPC_new_makefile = /tmp/CMT_Calc_MPPC_new$(cmt_lock_pid).make
else
cmt_final_setup_Calc_MPPC_new = $(bin)INGRID_Calc_MPPC_newsetup.make
cmt_local_Calc_MPPC_new_makefile = $(bin)Calc_MPPC_new.make
endif

else

ifdef READONLY
cmt_final_setup_Calc_MPPC_new = /tmp/CMT_INGRIDsetup.make
cmt_local_Calc_MPPC_new_makefile = /tmp/CMT_Calc_MPPC_new$(cmt_lock_pid).make
else
cmt_final_setup_Calc_MPPC_new = $(bin)INGRIDsetup.make
cmt_local_Calc_MPPC_new_makefile = $(bin)Calc_MPPC_new.make
endif

endif

ifdef READONLY
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
cmt_final_setup = $(bin)INGRIDsetup.make
endif

Calc_MPPC_new ::


ifdef READONLY
Calc_MPPC_new ::
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
	$(echo) 'Calc_MPPC_new'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = Calc_MPPC_new/
Calc_MPPC_new::
	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

#-- end of make_header ------------------

#-- start of application_header

Calc_MPPC_new :: dirs  $(bin)Calc_MPPC_new${application_suffix}
	$(echo) "Calc_MPPC_new ok"

#-- end of application_header
#-- start of application

$(bin)Calc_MPPC_new${application_suffix} :: $(bin)Calc_MPPC_new.o $(use_stamps) $(Calc_MPPC_newstamps) requirements $(use_requirements)
	$(link_echo) "application $@"
	$(link_silent) $(cpplink) -o $(@).new $(bin)Calc_MPPC_new.o $(cmt_installarea_linkopts) $(Calc_MPPC_new_use_linkopts) $(Calc_MPPC_newlinkopts) && mv -f $(@).new $(@)

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

install_dir = ${CMTINSTALLAREA}/$(tag)/bin
Calc_MPPC_newinstallname = Calc_MPPC_new${application_suffix}

Calc_MPPC_new :: Calc_MPPC_newinstall

install :: Calc_MPPC_newinstall

Calc_MPPC_newinstall :: $(install_dir)/$(Calc_MPPC_newinstallname)
ifdef CMTINSTALLAREA
	$(echo) "installation done"
endif

$(install_dir)/$(Calc_MPPC_newinstallname) :: $(bin)$(Calc_MPPC_newinstallname)
ifdef CMTINSTALLAREA
	$(install_silent) $(cmt_install_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(Calc_MPPC_newinstallname)" \
	    -out "$(install_dir)" \
	    -cmd "$(cmt_installarea_command)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

##Calc_MPPC_newclean :: Calc_MPPC_newuninstall

uninstall :: Calc_MPPC_newuninstall

Calc_MPPC_newuninstall ::
ifdef CMTINSTALLAREA
	$(cleanup_silent) $(cmt_uninstall_action) \
	    -source "`(cd $(bin); pwd)`" \
	    -name "$(Calc_MPPC_newinstallname)" \
	    -out "$(install_dir)" \
	    -cmtpath "$($(package)_cmtpath)"
endif

#	@echo "------> (Calc_MPPC_new.make) Removing installed files"
#-- end of application
#-- start of dependency ------------------
ifneq ($(MAKECMDGOALS),Calc_MPPC_newclean)

#$(bin)Calc_MPPC_new_dependencies.make :: dirs

ifndef QUICK
$(bin)Calc_MPPC_new_dependencies.make : ../app/Calc_MPPC_new.cxx $(use_requirements) $(cmt_final_setup_Calc_MPPC_new)
	$(echo) "(Calc_MPPC_new.make) Rebuilding $@"; \
	  $(build_dependencies) Calc_MPPC_new -all_sources -out=$@ ../app/Calc_MPPC_new.cxx
endif

#$(Calc_MPPC_new_dependencies)

-include $(bin)Calc_MPPC_new_dependencies.make

endif
#-- end of dependency -------------------
#-- start of cpp ------

$(bin)Calc_MPPC_new_dependencies.make : $(Calc_MPPC_new_cxx_dependencies)

$(bin)$(binobj)Calc_MPPC_new.o : $(Calc_MPPC_new_cxx_dependencies)
	$(cpp_echo) ../app/Calc_MPPC_new.cxx
	$(cpp_silent) $(cppcomp) -o $(@) $(use_pp_cppflags) $(Calc_MPPC_new_pp_cppflags) $(app_Calc_MPPC_new_pp_cppflags) $(Calc_MPPC_new_pp_cppflags) $(use_cppflags) $(Calc_MPPC_new_cppflags) $(app_Calc_MPPC_new_cppflags) $(Calc_MPPC_new_cppflags) $(Calc_MPPC_new_cxx_cppflags) -I../app ../app/Calc_MPPC_new.cxx

#-- end of cpp ------
#-- start of cleanup_header --------------

clean :: Calc_MPPC_newclean
	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(Calc_MPPC_new.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: No rule for such target" >&2; exit
else
.DEFAULT::
	$(echo) "(Calc_MPPC_new.make) PEDANTIC: $@: No rule for such target" >&2
	if test $@ = "$(cmt_final_setup)" -o\
	 $@ = "$(cmt_final_setup_Calc_MPPC_new)" ; then\
	 found=n; for s in 1 2 3 4 5; do\
	 sleep $$s; test ! -f $@ || { found=y; break; }\
	 done; if test $$found = n; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(Calc_MPPC_new.make) PEDANTIC: $@: Seems to be missing. Ignore it for now" >&2; exit 0 ; fi;\
	 elif test `expr index $@ '/'` -ne 0 ; then\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(Calc_MPPC_new.make) PEDANTIC: $@: Seems to be a missing file. Please check" >&2; exit 2 ; \
	 else\
	 test -z "$(cmtmsg)" ||\
	 echo "$(CMTMSGPREFIX)" "(Calc_MPPC_new.make) PEDANTIC: $@: Seems to be a fake target due to some pattern. Just ignore it" >&2 ; exit 0; fi
endif

Calc_MPPC_newclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_application ------
	$(cleanup_echo) Calc_MPPC_new${application_suffix}
	-$(cleanup_silent) cd $(bin); /bin/rm -f Calc_MPPC_new${application_suffix}

#	@echo "------> (Calc_MPPC_new.make) Removing application files"
#-- end of cleanup_application ------
#-- start of cleanup_objects ------
	$(cleanup_echo) objects
	-$(cleanup_silent) /bin/rm -f $(bin)Calc_MPPC_new.o
	-$(cleanup_silent) cd $(bin); /bin/rm -rf Calc_MPPC_new_deps Calc_MPPC_new_dependencies.make
#-- end of cleanup_objects ------
