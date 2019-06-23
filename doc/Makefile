# contact_trace
### Hooks for the editor to set the default target

current: target
-include target.mk

##################################################################

## uplink-style makestuff

Sources += Makefile
Ignore += .gitignore

msrepo = https://github.com/dushoff
ms = makestuff
makestuff:
	/bin/ln -fs ../makestuff .

-include $(ms)/os.mk

Ignore += $(ms)
Makefile: $(ms) $(ms)/Makefile

######################################################################

## Notes

Sources += notes.md

## Manuscript

Sources += network.tex
network.pdf: network.tex

Sources += network.bib

######################################################################

### Makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/texdeps.mk

