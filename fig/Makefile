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

compare_methods.Rout: compare_methods.R

smallpox.Rout: smallpox.R

Sources += looksmall.R
looksmall.Rout: ../data/Smallpox_1951_Tilburg_NL.csv looksmall.R

######################################################################

### Makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/wrapR.mk
