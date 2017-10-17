# Flexible "out of tree" build thanks to Paul D. Smith <psmith@gnu.org>
# http://mad-scientist.net/make/multi-arch.html

.SUFFIXES:

.PHONY: clean $(bldpath)

$(bldpath): $(bldpath)/flags.mk
	+@[ -d $@ ] || mkdir -p $@
#	echo $(MAKEFLAGS)
	+@$(MAKE) -C $@ -f $(CURDIR)/Makefile srcpath=$(CURDIR) $(MAKECMDGOALS)
#	--no-print-directory
#	[ -d $@ ] || mkdir -p $@
#	$(MAKE) --no-print-directory -C $@ -f $(CURDIR)/Makefile srcpath=$(CURDIR) $(MAKECMDGOALS)


Makefile:;


%.mk :: ;

% :: $(bldpath); :

clean:
	rm -rf $(bldpath)
