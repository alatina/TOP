.PHONY: placet placet-htgen utils

placet:
	+make -f Makefile.default

placet-htgen:
	+make -f Makefile.htgen

utils:
	+make -f Makefile.utils

cleanobj:
	@rm -rf placet.obj
	@rm -rf placet-htgen.obj
	@rm -f octave.df

cleanexe:
	@rm -f placet
	@rm -f placet-htgen
	@rm -f ground
	@rm -f grid
	@rm -f mad2gp

clean: cleanobj cleanexe
