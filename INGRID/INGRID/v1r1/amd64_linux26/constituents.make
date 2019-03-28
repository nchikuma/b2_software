
#-- start of constituents_header ------

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

tags      = $(tag),$(CMTEXTRATAGS)

INGRID_tag = $(tag)

ifdef READONLY
cmt_local_tagfile = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
else
#cmt_local_tagfile = $(INGRID_tag).make
cmt_local_tagfile = $(bin)$(INGRID_tag).make
endif

#-include $(cmt_local_tagfile)
include $(cmt_local_tagfile)

ifdef READONLY
cmt_local_setup = /tmp/CMT_INGRIDsetup$(cmt_lock_pid).make
cmt_final_setup = /tmp/CMT_INGRIDsetup.make
else
#cmt_local_setup = $(bin)INGRIDsetup$(cmt_lock_pid).make
cmt_local_setup = $(bin)$(package)setup$$$$.make
#cmt_final_setup = $(bin)INGRIDsetup.make
cmt_final_setup = $(bin)$(package)setup.make
endif

#--------------------------------------------------------

#cmt_lock_setup = /tmp/lock$(cmt_lock_pid).make
#cmt_temp_tag = /tmp/tag$(cmt_lock_pid).make

#first :: $(cmt_local_tagfile)
#	@echo $(cmt_local_tagfile) ok
ifndef QUICK
first :: $(cmt_final_setup) ;
else
first :: ;
endif

##	@bin=`$(cmtexe) show macro_value bin`

#$(cmt_local_tagfile) : $(cmt_lock_setup)
#	@echo "#CMT> Error: $@: No such file" >&2; exit 1
$(cmt_local_tagfile) :
	@echo "#CMT> Warning: $@: No such file" >&2; exit
#	@echo "#CMT> Info: $@: No need to rebuild file" >&2; exit

$(cmt_final_setup) : $(cmt_local_tagfile) 
	$(echo) "(constituents.make) Rebuilding $@"
	@if test ! -d $(@D); then $(mkdir) -p $(@D); fi; \
	  if test -f $(cmt_local_setup); then /bin/rm -f $(cmt_local_setup); fi; \
	  trap '/bin/rm -f $(cmt_local_setup)' 0 1 2 15; \
	  $(cmtexe) -tag=$(tags) show setup >>$(cmt_local_setup); \
	  if test ! -f $@; then \
	    mv $(cmt_local_setup) $@; \
	  else \
	    if /usr/bin/diff $(cmt_local_setup) $@ >/dev/null ; then \
	      : ; \
	    else \
	      mv $(cmt_local_setup) $@; \
	    fi; \
	  fi

#	@/bin/echo $@ ok   

config :: checkuses
	@exit 0
checkuses : ;

env.make ::
	printenv >env.make.tmp; $(cmtexe) check files env.make.tmp env.make

ifndef QUICK
all :: build_library_links
	$(echo) "(constituents.make) all done"
endif

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

build_library_links : dirs requirements
	$(echo) "(constituents.make) Rebuilding library links"; \
	if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi; \
	$(build_library_links)

.DEFAULT ::
	$(echo) "(constituents.make) $@: No rule for such target" >&2
#	@echo "#CMT> Warning: $@: Using default commands" >&2; exit

#	@if test "$@" = "$(cmt_lock_setup)"; then \
	#  /bin/rm -f $(cmt_lock_setup); \
	 # touch $(cmt_lock_setup); \
	#fi

#-- end of constituents_header ------
#-- start of group ------

all_groups :: all

all :: $(all_dependencies)  $(all_pre_constituents) $(all_constituents)  $(all_post_constituents)
	$(echo) "all ok."

#	@/bin/echo " all ok."

clean :: allclean

allclean ::  $(all_constituentsclean)
	$(echo) $(all_constituentsclean)
	$(echo) "allclean ok."

#	@echo $(all_constituentsclean)
#	@/bin/echo " allclean ok."

allclean :: makefilesclean

#-- end of group ------
#-- start of group ------

all_groups :: cmt_actions

cmt_actions :: $(cmt_actions_dependencies)  $(cmt_actions_pre_constituents) $(cmt_actions_constituents)  $(cmt_actions_post_constituents)
	$(echo) "cmt_actions ok."

#	@/bin/echo " cmt_actions ok."

clean :: allclean

cmt_actionsclean ::  $(cmt_actions_constituentsclean)
	$(echo) $(cmt_actions_constituentsclean)
	$(echo) "cmt_actionsclean ok."

#	@echo $(cmt_actions_constituentsclean)
#	@/bin/echo " cmt_actionsclean ok."

cmt_actionsclean :: makefilesclean

#-- end of group ------
#-- start of group ------

all_groups :: documentation

documentation :: $(documentation_dependencies)  $(documentation_pre_constituents) $(documentation_constituents)  $(documentation_post_constituents)
	$(echo) "documentation ok."

#	@/bin/echo " documentation ok."

clean :: allclean

documentationclean ::  $(documentation_constituentsclean)
	$(echo) $(documentation_constituentsclean)
	$(echo) "documentationclean ok."

#	@echo $(documentation_constituentsclean)
#	@/bin/echo " documentationclean ok."

documentationclean :: makefilesclean

#-- end of group ------
#-- start of constituent_lock ------

cmt_dictionary_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_dictionary_has_target_tag

ifdef READONLY
cmt_local_tagfile_dictionary = /tmp/CMT_$(INGRID_tag)_dictionary.make$(cmt_lock_pid)
cmt_final_setup_dictionary = /tmp/CMT_INGRID_dictionarysetup.make
cmt_local_dictionary_makefile = /tmp/CMT_dictionary$(cmt_lock_pid).make
else
#cmt_local_tagfile_dictionary = $(INGRID_tag)_dictionary.make
cmt_local_tagfile_dictionary = $(bin)$(INGRID_tag)_dictionary.make
cmt_final_setup_dictionary = $(bin)INGRID_dictionarysetup.make
cmt_local_dictionary_makefile = $(bin)dictionary.make
endif

dictionary_extratags = -tag_add=target_dictionary

#$(cmt_local_tagfile_dictionary) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_dictionary) ::
else
$(cmt_local_tagfile_dictionary) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_dictionary)"
	@if test -f $(cmt_local_tagfile_dictionary); then /bin/rm -f $(cmt_local_tagfile_dictionary); fi ; \
	  $(cmtexe) -tag=$(tags) $(dictionary_extratags) build tag_makefile >>$(cmt_local_tagfile_dictionary); \
	  if test -f $(cmt_final_setup_dictionary); then /bin/rm -f $(cmt_final_setup_dictionary); fi; \
	  $(cmtexe) -tag=$(tags) $(dictionary_extratags) show setup >>$(cmt_final_setup_dictionary)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_dictionary = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_dictionary = /tmp/CMT_INGRIDsetup.make
cmt_local_dictionary_makefile = /tmp/CMT_dictionary$(cmt_lock_pid).make
else
#cmt_local_tagfile_dictionary = $(INGRID_tag).make
cmt_local_tagfile_dictionary = $(bin)$(INGRID_tag).make
cmt_final_setup_dictionary = $(bin)INGRIDsetup.make
cmt_local_dictionary_makefile = $(bin)dictionary.make
endif

endif

ifndef QUICK
$(cmt_local_dictionary_makefile) :: $(dictionary_dependencies) $(cmt_local_tagfile_dictionary) build_library_links dirs
else
$(cmt_local_dictionary_makefile) :: $(cmt_local_tagfile_dictionary)
endif
	$(echo) "(constituents.make) Building dictionary.make"; \
	  $(cmtexe) -tag=$(tags) $(dictionary_extratags) build constituent_makefile -out=$(cmt_local_dictionary_makefile) dictionary

dictionary :: $(dictionary_dependencies) $(cmt_local_dictionary_makefile)
	$(echo) "(constituents.make) Creating dictionary${lock_suffix} and Starting dictionary"
	@${lock_command} dictionary${lock_suffix} || exit $$?; \
	  retval=$$?; \
	  trap '${unlock_command} dictionary${lock_suffix}; exit $${retval}' 1 2 15; \
	  $(MAKE) -f $(cmt_local_dictionary_makefile) cmt_lock_pid=$${cmt_lock_pid} dictionary; \
	  retval=$$?; ${unlock_command} dictionary${lock_suffix}; exit $${retval}
	$(echo) "(constituents.make) dictionary done"

clean :: dictionaryclean

dictionaryclean :: $(dictionaryclean_dependencies) ##$(cmt_local_dictionary_makefile)
	$(echo) "(constituents.make) Starting dictionaryclean"
	@-if test -f $(cmt_local_dictionary_makefile); then \
	  $(MAKE) -f $(cmt_local_dictionary_makefile) cmt_lock_pid=$${cmt_lock_pid} dictionaryclean; \
	fi

##	  /bin/rm -f $(cmt_local_dictionary_makefile) $(bin)dictionary_dependencies.make

install :: dictionaryinstall

dictionaryinstall :: $(dictionary_dependencies) $(cmt_local_dictionary_makefile)
	$(echo) "(constituents.make) Starting install dictionary"
	@-$(MAKE) -f $(cmt_local_dictionary_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install dictionary done"

uninstall :: dictionaryuninstall

dictionaryuninstall :: $(cmt_local_dictionary_makefile)
	$(echo) "(constituents.make) Starting uninstall dictionary"
	@-$(MAKE) -f $(cmt_local_dictionary_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall dictionary done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ dictionary"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ dictionary done"
endif


#-- end of constituent_lock ------
#-- start of constituent_lock ------

cmt_version_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_version_has_target_tag

ifdef READONLY
cmt_local_tagfile_version = /tmp/CMT_$(INGRID_tag)_version.make$(cmt_lock_pid)
cmt_final_setup_version = /tmp/CMT_INGRID_versionsetup.make
cmt_local_version_makefile = /tmp/CMT_version$(cmt_lock_pid).make
else
#cmt_local_tagfile_version = $(INGRID_tag)_version.make
cmt_local_tagfile_version = $(bin)$(INGRID_tag)_version.make
cmt_final_setup_version = $(bin)INGRID_versionsetup.make
cmt_local_version_makefile = $(bin)version.make
endif

version_extratags = -tag_add=target_version

#$(cmt_local_tagfile_version) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_version) ::
else
$(cmt_local_tagfile_version) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_version)"
	@if test -f $(cmt_local_tagfile_version); then /bin/rm -f $(cmt_local_tagfile_version); fi ; \
	  $(cmtexe) -tag=$(tags) $(version_extratags) build tag_makefile >>$(cmt_local_tagfile_version); \
	  if test -f $(cmt_final_setup_version); then /bin/rm -f $(cmt_final_setup_version); fi; \
	  $(cmtexe) -tag=$(tags) $(version_extratags) show setup >>$(cmt_final_setup_version)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_version = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_version = /tmp/CMT_INGRIDsetup.make
cmt_local_version_makefile = /tmp/CMT_version$(cmt_lock_pid).make
else
#cmt_local_tagfile_version = $(INGRID_tag).make
cmt_local_tagfile_version = $(bin)$(INGRID_tag).make
cmt_final_setup_version = $(bin)INGRIDsetup.make
cmt_local_version_makefile = $(bin)version.make
endif

endif

ifndef QUICK
$(cmt_local_version_makefile) :: $(version_dependencies) $(cmt_local_tagfile_version) build_library_links dirs
else
$(cmt_local_version_makefile) :: $(cmt_local_tagfile_version)
endif
	$(echo) "(constituents.make) Building version.make"; \
	  $(cmtexe) -tag=$(tags) $(version_extratags) build constituent_makefile -out=$(cmt_local_version_makefile) version

version :: $(version_dependencies) $(cmt_local_version_makefile)
	$(echo) "(constituents.make) Creating version${lock_suffix} and Starting version"
	@${lock_command} version${lock_suffix} || exit $$?; \
	  retval=$$?; \
	  trap '${unlock_command} version${lock_suffix}; exit $${retval}' 1 2 15; \
	  $(MAKE) -f $(cmt_local_version_makefile) cmt_lock_pid=$${cmt_lock_pid} version; \
	  retval=$$?; ${unlock_command} version${lock_suffix}; exit $${retval}
	$(echo) "(constituents.make) version done"

clean :: versionclean

versionclean :: $(versionclean_dependencies) ##$(cmt_local_version_makefile)
	$(echo) "(constituents.make) Starting versionclean"
	@-if test -f $(cmt_local_version_makefile); then \
	  $(MAKE) -f $(cmt_local_version_makefile) cmt_lock_pid=$${cmt_lock_pid} versionclean; \
	fi

##	  /bin/rm -f $(cmt_local_version_makefile) $(bin)version_dependencies.make

install :: versioninstall

versioninstall :: $(version_dependencies) $(cmt_local_version_makefile)
	$(echo) "(constituents.make) Starting install version"
	@-$(MAKE) -f $(cmt_local_version_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install version done"

uninstall :: versionuninstall

versionuninstall :: $(cmt_local_version_makefile)
	$(echo) "(constituents.make) Starting uninstall version"
	@-$(MAKE) -f $(cmt_local_version_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall version done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ version"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ version done"
endif


#-- end of constituent_lock ------
#-- start of constituent_lock ------

cmt_register_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_register_has_target_tag

ifdef READONLY
cmt_local_tagfile_register = /tmp/CMT_$(INGRID_tag)_register.make$(cmt_lock_pid)
cmt_final_setup_register = /tmp/CMT_INGRID_registersetup.make
cmt_local_register_makefile = /tmp/CMT_register$(cmt_lock_pid).make
else
#cmt_local_tagfile_register = $(INGRID_tag)_register.make
cmt_local_tagfile_register = $(bin)$(INGRID_tag)_register.make
cmt_final_setup_register = $(bin)INGRID_registersetup.make
cmt_local_register_makefile = $(bin)register.make
endif

register_extratags = -tag_add=target_register

#$(cmt_local_tagfile_register) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_register) ::
else
$(cmt_local_tagfile_register) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_register)"
	@if test -f $(cmt_local_tagfile_register); then /bin/rm -f $(cmt_local_tagfile_register); fi ; \
	  $(cmtexe) -tag=$(tags) $(register_extratags) build tag_makefile >>$(cmt_local_tagfile_register); \
	  if test -f $(cmt_final_setup_register); then /bin/rm -f $(cmt_final_setup_register); fi; \
	  $(cmtexe) -tag=$(tags) $(register_extratags) show setup >>$(cmt_final_setup_register)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_register = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_register = /tmp/CMT_INGRIDsetup.make
cmt_local_register_makefile = /tmp/CMT_register$(cmt_lock_pid).make
else
#cmt_local_tagfile_register = $(INGRID_tag).make
cmt_local_tagfile_register = $(bin)$(INGRID_tag).make
cmt_final_setup_register = $(bin)INGRIDsetup.make
cmt_local_register_makefile = $(bin)register.make
endif

endif

ifndef QUICK
$(cmt_local_register_makefile) :: $(register_dependencies) $(cmt_local_tagfile_register) build_library_links dirs
else
$(cmt_local_register_makefile) :: $(cmt_local_tagfile_register)
endif
	$(echo) "(constituents.make) Building register.make"; \
	  $(cmtexe) -tag=$(tags) $(register_extratags) build constituent_makefile -out=$(cmt_local_register_makefile) register

register :: $(register_dependencies) $(cmt_local_register_makefile)
	$(echo) "(constituents.make) Creating register${lock_suffix} and Starting register"
	@${lock_command} register${lock_suffix} || exit $$?; \
	  retval=$$?; \
	  trap '${unlock_command} register${lock_suffix}; exit $${retval}' 1 2 15; \
	  $(MAKE) -f $(cmt_local_register_makefile) cmt_lock_pid=$${cmt_lock_pid} register; \
	  retval=$$?; ${unlock_command} register${lock_suffix}; exit $${retval}
	$(echo) "(constituents.make) register done"

clean :: registerclean

registerclean :: $(registerclean_dependencies) ##$(cmt_local_register_makefile)
	$(echo) "(constituents.make) Starting registerclean"
	@-if test -f $(cmt_local_register_makefile); then \
	  $(MAKE) -f $(cmt_local_register_makefile) cmt_lock_pid=$${cmt_lock_pid} registerclean; \
	fi

##	  /bin/rm -f $(cmt_local_register_makefile) $(bin)register_dependencies.make

install :: registerinstall

registerinstall :: $(register_dependencies) $(cmt_local_register_makefile)
	$(echo) "(constituents.make) Starting install register"
	@-$(MAKE) -f $(cmt_local_register_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install register done"

uninstall :: registeruninstall

registeruninstall :: $(cmt_local_register_makefile)
	$(echo) "(constituents.make) Starting uninstall register"
	@-$(MAKE) -f $(cmt_local_register_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall register done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ register"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ register done"
endif


#-- end of constituent_lock ------
#-- start of constituent_lock ------

cmt_doxygen_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_doxygen_has_target_tag

ifdef READONLY
cmt_local_tagfile_doxygen = /tmp/CMT_$(INGRID_tag)_doxygen.make$(cmt_lock_pid)
cmt_final_setup_doxygen = /tmp/CMT_INGRID_doxygensetup.make
cmt_local_doxygen_makefile = /tmp/CMT_doxygen$(cmt_lock_pid).make
else
#cmt_local_tagfile_doxygen = $(INGRID_tag)_doxygen.make
cmt_local_tagfile_doxygen = $(bin)$(INGRID_tag)_doxygen.make
cmt_final_setup_doxygen = $(bin)INGRID_doxygensetup.make
cmt_local_doxygen_makefile = $(bin)doxygen.make
endif

doxygen_extratags = -tag_add=target_doxygen

#$(cmt_local_tagfile_doxygen) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_doxygen) ::
else
$(cmt_local_tagfile_doxygen) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_doxygen)"
	@if test -f $(cmt_local_tagfile_doxygen); then /bin/rm -f $(cmt_local_tagfile_doxygen); fi ; \
	  $(cmtexe) -tag=$(tags) $(doxygen_extratags) build tag_makefile >>$(cmt_local_tagfile_doxygen); \
	  if test -f $(cmt_final_setup_doxygen); then /bin/rm -f $(cmt_final_setup_doxygen); fi; \
	  $(cmtexe) -tag=$(tags) $(doxygen_extratags) show setup >>$(cmt_final_setup_doxygen)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_doxygen = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_doxygen = /tmp/CMT_INGRIDsetup.make
cmt_local_doxygen_makefile = /tmp/CMT_doxygen$(cmt_lock_pid).make
else
#cmt_local_tagfile_doxygen = $(INGRID_tag).make
cmt_local_tagfile_doxygen = $(bin)$(INGRID_tag).make
cmt_final_setup_doxygen = $(bin)INGRIDsetup.make
cmt_local_doxygen_makefile = $(bin)doxygen.make
endif

endif

ifndef QUICK
$(cmt_local_doxygen_makefile) :: $(doxygen_dependencies) $(cmt_local_tagfile_doxygen) build_library_links dirs
else
$(cmt_local_doxygen_makefile) :: $(cmt_local_tagfile_doxygen)
endif
	$(echo) "(constituents.make) Building doxygen.make"; \
	  $(cmtexe) -tag=$(tags) $(doxygen_extratags) build constituent_makefile -out=$(cmt_local_doxygen_makefile) doxygen

doxygen :: $(doxygen_dependencies) $(cmt_local_doxygen_makefile)
	$(echo) "(constituents.make) Creating doxygen${lock_suffix} and Starting doxygen"
	@${lock_command} doxygen${lock_suffix} || exit $$?; \
	  retval=$$?; \
	  trap '${unlock_command} doxygen${lock_suffix}; exit $${retval}' 1 2 15; \
	  $(MAKE) -f $(cmt_local_doxygen_makefile) cmt_lock_pid=$${cmt_lock_pid} doxygen; \
	  retval=$$?; ${unlock_command} doxygen${lock_suffix}; exit $${retval}
	$(echo) "(constituents.make) doxygen done"

clean :: doxygenclean

doxygenclean :: $(doxygenclean_dependencies) ##$(cmt_local_doxygen_makefile)
	$(echo) "(constituents.make) Starting doxygenclean"
	@-if test -f $(cmt_local_doxygen_makefile); then \
	  $(MAKE) -f $(cmt_local_doxygen_makefile) cmt_lock_pid=$${cmt_lock_pid} doxygenclean; \
	fi

##	  /bin/rm -f $(cmt_local_doxygen_makefile) $(bin)doxygen_dependencies.make

install :: doxygeninstall

doxygeninstall :: $(doxygen_dependencies) $(cmt_local_doxygen_makefile)
	$(echo) "(constituents.make) Starting install doxygen"
	@-$(MAKE) -f $(cmt_local_doxygen_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install doxygen done"

uninstall :: doxygenuninstall

doxygenuninstall :: $(cmt_local_doxygen_makefile)
	$(echo) "(constituents.make) Starting uninstall doxygen"
	@-$(MAKE) -f $(cmt_local_doxygen_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall doxygen done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ doxygen"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ doxygen done"
endif


#-- end of constituent_lock ------
#-- start of constituent ------

cmt_DSTMaker_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_DSTMaker_has_target_tag

ifdef READONLY
cmt_local_tagfile_DSTMaker = /tmp/CMT_$(INGRID_tag)_DSTMaker.make$(cmt_lock_pid)
cmt_final_setup_DSTMaker = /tmp/CMT_INGRID_DSTMakersetup.make
cmt_local_DSTMaker_makefile = /tmp/CMT_DSTMaker$(cmt_lock_pid).make
else
#cmt_local_tagfile_DSTMaker = $(INGRID_tag)_DSTMaker.make
cmt_local_tagfile_DSTMaker = $(bin)$(INGRID_tag)_DSTMaker.make
cmt_final_setup_DSTMaker = $(bin)INGRID_DSTMakersetup.make
cmt_local_DSTMaker_makefile = $(bin)DSTMaker.make
endif

DSTMaker_extratags = -tag_add=target_DSTMaker

#$(cmt_local_tagfile_DSTMaker) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_DSTMaker) ::
else
$(cmt_local_tagfile_DSTMaker) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_DSTMaker)"
	@if test -f $(cmt_local_tagfile_DSTMaker); then /bin/rm -f $(cmt_local_tagfile_DSTMaker); fi ; \
	  $(cmtexe) -tag=$(tags) $(DSTMaker_extratags) build tag_makefile >>$(cmt_local_tagfile_DSTMaker); \
	  if test -f $(cmt_final_setup_DSTMaker); then /bin/rm -f $(cmt_final_setup_DSTMaker); fi; \
	  $(cmtexe) -tag=$(tags) $(DSTMaker_extratags) show setup >>$(cmt_final_setup_DSTMaker)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_DSTMaker = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_DSTMaker = /tmp/CMT_INGRIDsetup.make
cmt_local_DSTMaker_makefile = /tmp/CMT_DSTMaker$(cmt_lock_pid).make
else
#cmt_local_tagfile_DSTMaker = $(INGRID_tag).make
cmt_local_tagfile_DSTMaker = $(bin)$(INGRID_tag).make
cmt_final_setup_DSTMaker = $(bin)INGRIDsetup.make
cmt_local_DSTMaker_makefile = $(bin)DSTMaker.make
endif

endif

ifndef QUICK
$(cmt_local_DSTMaker_makefile) :: $(DSTMaker_dependencies) $(cmt_local_tagfile_DSTMaker) build_library_links dirs
else
$(cmt_local_DSTMaker_makefile) :: $(cmt_local_tagfile_DSTMaker)
endif
	$(echo) "(constituents.make) Building DSTMaker.make"; \
	  $(cmtexe) -tag=$(tags) $(DSTMaker_extratags) build constituent_makefile -out=$(cmt_local_DSTMaker_makefile) DSTMaker

DSTMaker :: $(DSTMaker_dependencies) $(cmt_local_DSTMaker_makefile)
	$(echo) "(constituents.make) Starting DSTMaker"
	@$(MAKE) -f $(cmt_local_DSTMaker_makefile) cmt_lock_pid=$${cmt_lock_pid} DSTMaker
	$(echo) "(constituents.make) DSTMaker done"

clean :: DSTMakerclean

DSTMakerclean :: $(DSTMakerclean_dependencies) ##$(cmt_local_DSTMaker_makefile)
	$(echo) "(constituents.make) Starting DSTMakerclean"
	@-if test -f $(cmt_local_DSTMaker_makefile); then \
	  $(MAKE) -f $(cmt_local_DSTMaker_makefile) cmt_lock_pid=$${cmt_lock_pid} DSTMakerclean; \
	fi

##	  /bin/rm -f $(cmt_local_DSTMaker_makefile) $(bin)DSTMaker_dependencies.make

install :: DSTMakerinstall

DSTMakerinstall :: $(DSTMaker_dependencies) $(cmt_local_DSTMaker_makefile)
	$(echo) "(constituents.make) Starting install DSTMaker"
	@-$(MAKE) -f $(cmt_local_DSTMaker_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install DSTMaker done"

uninstall :: DSTMakeruninstall

DSTMakeruninstall :: $(cmt_local_DSTMaker_makefile)
	$(echo) "(constituents.make) Starting uninstall DSTMaker"
	@-$(MAKE) -f $(cmt_local_DSTMaker_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall DSTMaker done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ DSTMaker"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ DSTMaker done"
endif


#-- end of constituent ------
#-- start of constituent ------

cmt_Calc_MPPC_new_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_Calc_MPPC_new_has_target_tag

ifdef READONLY
cmt_local_tagfile_Calc_MPPC_new = /tmp/CMT_$(INGRID_tag)_Calc_MPPC_new.make$(cmt_lock_pid)
cmt_final_setup_Calc_MPPC_new = /tmp/CMT_INGRID_Calc_MPPC_newsetup.make
cmt_local_Calc_MPPC_new_makefile = /tmp/CMT_Calc_MPPC_new$(cmt_lock_pid).make
else
#cmt_local_tagfile_Calc_MPPC_new = $(INGRID_tag)_Calc_MPPC_new.make
cmt_local_tagfile_Calc_MPPC_new = $(bin)$(INGRID_tag)_Calc_MPPC_new.make
cmt_final_setup_Calc_MPPC_new = $(bin)INGRID_Calc_MPPC_newsetup.make
cmt_local_Calc_MPPC_new_makefile = $(bin)Calc_MPPC_new.make
endif

Calc_MPPC_new_extratags = -tag_add=target_Calc_MPPC_new

#$(cmt_local_tagfile_Calc_MPPC_new) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_Calc_MPPC_new) ::
else
$(cmt_local_tagfile_Calc_MPPC_new) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_Calc_MPPC_new)"
	@if test -f $(cmt_local_tagfile_Calc_MPPC_new); then /bin/rm -f $(cmt_local_tagfile_Calc_MPPC_new); fi ; \
	  $(cmtexe) -tag=$(tags) $(Calc_MPPC_new_extratags) build tag_makefile >>$(cmt_local_tagfile_Calc_MPPC_new); \
	  if test -f $(cmt_final_setup_Calc_MPPC_new); then /bin/rm -f $(cmt_final_setup_Calc_MPPC_new); fi; \
	  $(cmtexe) -tag=$(tags) $(Calc_MPPC_new_extratags) show setup >>$(cmt_final_setup_Calc_MPPC_new)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_Calc_MPPC_new = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_Calc_MPPC_new = /tmp/CMT_INGRIDsetup.make
cmt_local_Calc_MPPC_new_makefile = /tmp/CMT_Calc_MPPC_new$(cmt_lock_pid).make
else
#cmt_local_tagfile_Calc_MPPC_new = $(INGRID_tag).make
cmt_local_tagfile_Calc_MPPC_new = $(bin)$(INGRID_tag).make
cmt_final_setup_Calc_MPPC_new = $(bin)INGRIDsetup.make
cmt_local_Calc_MPPC_new_makefile = $(bin)Calc_MPPC_new.make
endif

endif

ifndef QUICK
$(cmt_local_Calc_MPPC_new_makefile) :: $(Calc_MPPC_new_dependencies) $(cmt_local_tagfile_Calc_MPPC_new) build_library_links dirs
else
$(cmt_local_Calc_MPPC_new_makefile) :: $(cmt_local_tagfile_Calc_MPPC_new)
endif
	$(echo) "(constituents.make) Building Calc_MPPC_new.make"; \
	  $(cmtexe) -tag=$(tags) $(Calc_MPPC_new_extratags) build constituent_makefile -out=$(cmt_local_Calc_MPPC_new_makefile) Calc_MPPC_new

Calc_MPPC_new :: $(Calc_MPPC_new_dependencies) $(cmt_local_Calc_MPPC_new_makefile)
	$(echo) "(constituents.make) Starting Calc_MPPC_new"
	@$(MAKE) -f $(cmt_local_Calc_MPPC_new_makefile) cmt_lock_pid=$${cmt_lock_pid} Calc_MPPC_new
	$(echo) "(constituents.make) Calc_MPPC_new done"

clean :: Calc_MPPC_newclean

Calc_MPPC_newclean :: $(Calc_MPPC_newclean_dependencies) ##$(cmt_local_Calc_MPPC_new_makefile)
	$(echo) "(constituents.make) Starting Calc_MPPC_newclean"
	@-if test -f $(cmt_local_Calc_MPPC_new_makefile); then \
	  $(MAKE) -f $(cmt_local_Calc_MPPC_new_makefile) cmt_lock_pid=$${cmt_lock_pid} Calc_MPPC_newclean; \
	fi

##	  /bin/rm -f $(cmt_local_Calc_MPPC_new_makefile) $(bin)Calc_MPPC_new_dependencies.make

install :: Calc_MPPC_newinstall

Calc_MPPC_newinstall :: $(Calc_MPPC_new_dependencies) $(cmt_local_Calc_MPPC_new_makefile)
	$(echo) "(constituents.make) Starting install Calc_MPPC_new"
	@-$(MAKE) -f $(cmt_local_Calc_MPPC_new_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install Calc_MPPC_new done"

uninstall :: Calc_MPPC_newuninstall

Calc_MPPC_newuninstall :: $(cmt_local_Calc_MPPC_new_makefile)
	$(echo) "(constituents.make) Starting uninstall Calc_MPPC_new"
	@-$(MAKE) -f $(cmt_local_Calc_MPPC_new_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall Calc_MPPC_new done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ Calc_MPPC_new"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ Calc_MPPC_new done"
endif


#-- end of constituent ------
#-- start of constituent_lock ------

cmt_make_has_target_tag = 1

#--------------------------------------------------------

ifdef cmt_make_has_target_tag

ifdef READONLY
cmt_local_tagfile_make = /tmp/CMT_$(INGRID_tag)_make.make$(cmt_lock_pid)
cmt_final_setup_make = /tmp/CMT_INGRID_makesetup.make
cmt_local_make_makefile = /tmp/CMT_make$(cmt_lock_pid).make
else
#cmt_local_tagfile_make = $(INGRID_tag)_make.make
cmt_local_tagfile_make = $(bin)$(INGRID_tag)_make.make
cmt_final_setup_make = $(bin)INGRID_makesetup.make
cmt_local_make_makefile = $(bin)make.make
endif

make_extratags = -tag_add=target_make

#$(cmt_local_tagfile_make) : $(cmt_lock_setup)
ifndef QUICK
$(cmt_local_tagfile_make) ::
else
$(cmt_local_tagfile_make) :
endif
	$(echo) "(constituents.make) Rebuilding setup.make $(cmt_local_tagfile_make)"
	@if test -f $(cmt_local_tagfile_make); then /bin/rm -f $(cmt_local_tagfile_make); fi ; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build tag_makefile >>$(cmt_local_tagfile_make); \
	  if test -f $(cmt_final_setup_make); then /bin/rm -f $(cmt_final_setup_make); fi; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) show setup >>$(cmt_final_setup_make)
	$(echo) setup.make ok

else

ifdef READONLY
cmt_local_tagfile_make = /tmp/CMT_$(INGRID_tag).make$(cmt_lock_pid)
cmt_final_setup_make = /tmp/CMT_INGRIDsetup.make
cmt_local_make_makefile = /tmp/CMT_make$(cmt_lock_pid).make
else
#cmt_local_tagfile_make = $(INGRID_tag).make
cmt_local_tagfile_make = $(bin)$(INGRID_tag).make
cmt_final_setup_make = $(bin)INGRIDsetup.make
cmt_local_make_makefile = $(bin)make.make
endif

endif

ifndef QUICK
$(cmt_local_make_makefile) :: $(make_dependencies) $(cmt_local_tagfile_make) build_library_links dirs
else
$(cmt_local_make_makefile) :: $(cmt_local_tagfile_make)
endif
	$(echo) "(constituents.make) Building make.make"; \
	  $(cmtexe) -tag=$(tags) $(make_extratags) build constituent_makefile -out=$(cmt_local_make_makefile) make

make :: $(make_dependencies) $(cmt_local_make_makefile)
	$(echo) "(constituents.make) Creating make${lock_suffix} and Starting make"
	@${lock_command} make${lock_suffix} || exit $$?; \
	  retval=$$?; \
	  trap '${unlock_command} make${lock_suffix}; exit $${retval}' 1 2 15; \
	  $(MAKE) -f $(cmt_local_make_makefile) cmt_lock_pid=$${cmt_lock_pid} make; \
	  retval=$$?; ${unlock_command} make${lock_suffix}; exit $${retval}
	$(echo) "(constituents.make) make done"

clean :: makeclean

makeclean :: $(makeclean_dependencies) ##$(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting makeclean"
	@-if test -f $(cmt_local_make_makefile); then \
	  $(MAKE) -f $(cmt_local_make_makefile) cmt_lock_pid=$${cmt_lock_pid} makeclean; \
	fi

##	  /bin/rm -f $(cmt_local_make_makefile) $(bin)make_dependencies.make

install :: makeinstall

makeinstall :: $(make_dependencies) $(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting install make"
	@-$(MAKE) -f $(cmt_local_make_makefile) cmt_lock_pid=$${cmt_lock_pid} install
	$(echo) "(constituents.make) install make done"

uninstall :: makeuninstall

makeuninstall :: $(cmt_local_make_makefile)
	$(echo) "(constituents.make) Starting uninstall make"
	@-$(MAKE) -f $(cmt_local_make_makefile) cmt_lock_pid=$${cmt_lock_pid} uninstall
	$(echo) "(constituents.make) uninstall make done"

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(constituents.make) Starting $@ make"
	$(echo) Using default action for $@
	$(echo) "(constituents.make) $@ make done"
endif


#-- end of constituent_lock ------
#-- start of constituents_trailer ------

clean :: remove_library_links

remove_library_links ::
	$(echo) "Removing library links"; \
	  $(remove_library_links)

makefilesclean ::

###	@/bin/rm -f checkuses

###	/bin/rm -f *.make*

clean :: makefilesclean

binclean :: clean
	$(echo) "Removing binary directory $(bin)"
	@if test ! "$(bin)" = "./"; then \
	  /bin/rm -rf $(bin); \
	fi

#-- end of constituents_trailer ------
