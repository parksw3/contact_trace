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

Sources += network.tex cover.txt
network.pdf: network.tex

network.ld.pdf: network.tex
network.tex.2f572ca5a.oldfile: 
Ignore += *.ld.out *.ld.pdf

Sources += response.tex
## response.pdf: response.tex

## resources/Review.pdf:
resources: dir=~/Dropbox/park_trace
resources:
	$(linkdirname)
resources/%.pdf:
	$(MAKE) resources

Sources += network.bib

## network.tex.2f572ca5a.oldfile.texcount:
## network.tex.texcount:
Ignore += *.texcount
%.texcount: %
	texcount $< > $@

######################################################################

### Makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/texdeps.mk

