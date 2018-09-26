# contact_trace
### Hooks for the editor to set the default target

current: target
-include target.mk

##################################################################

## Defs

# stuff

Sources += Makefile .ignore .gitignore

Ignore += $(ms)
ms = makestuff
makestuff:
	/bin/ln -fs ../makestuff .
-include $(ms)/os.mk

######################################################################

### Makestuff

-include $(ms)/git.mk
-include $(ms)/visual.mk

-include $(ms)/texdeps.mk

