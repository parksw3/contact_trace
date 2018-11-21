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
$(ms): 
	cd .. && $(MAKE) makestuff && $(MAKE) makestuff/Makefile
	ln -s ../makestuff .

-include $(ms)/os.mk
# -include $(ms)/perl.def

######################################################################

## Content

compare_methods.Rout: compare_methods.R

######################################################################

### Makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/stepR.mk

