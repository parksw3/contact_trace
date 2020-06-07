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

## Sources += network.tex cover.txt
network.pdf: network.tex

network.ld.pdf: network.tex
network.tex.2f572ca5a.oldfile: 
Ignore += *.ld.out *.ld.pdf

## First journal had pretty nice reviews but thought we weren't sexy enough
## Sources += response.tex
## response.pdf: response.tex

## JRSI also nice reviews but one of them _very_ demanding
## A letter to the editor about scoping and timing for revision
# Sources += jscope.txt

## resources/Review.pdf:
Ignore += resources
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

