PROCESS = Decoder wgMerge_inglib wgMerge_spill wgSpillCorr wgMakeHist
MDIR = process

define maker
$(1):
	@echo "make $1:";	cd $(MDIR)/$1 && make
endef

define clean
	@echo "clean:";	cd $(MDIR)/$1 && rm -f obj/* *.d *.o ; cd ../../
endef

.PHONY: all clean single

all:$(PROCESS)
ARG = test

$(foreach tar,$(PROCESS),$(eval $(call maker,$(tar))))

clean:
	$(foreach tar,$(PROCESS),$(call clean,$(tar)))
	rm -f obj/*	
	rm -f bin/*

single: $(ARG)
$(eval $(call maker,$(ARG)))
