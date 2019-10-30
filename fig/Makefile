# contact_trace/R
### Hooks for the editor to set the default target

current: target
-include target.mk

##################################################################

## Defs

# stuff

Sources += Makefile

msrepo = https://github.com/dushoff
ms = makestuff
Ignore += $(ms)
Makefile: $(ms)/Makefile
$(ms)/Makefile: 
	cd .. && $(MAKE) makestuff && $(MAKE) makestuff/Makefile
	ln -s ../makestuff .

-include $(ms)/os.mk
# -include $(ms)/perl.def

######################################################################

## Content

## Deleting old rules corresponding to files Daniel already deleted
## 2019 Oct 30 (Wed)

## This is a hybrid (mostly uses Daniel workflow; (source and so on))
## Makes paper figure and talk figures together
compare_methods.Rout: compare_methods.R

Ignore += cmp_reproductive.pdf

######################################################################

### Makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/wrapR.mk
