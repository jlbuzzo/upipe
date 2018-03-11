###############################################################################
#
# A Makefile to rule'em all!
#
###############################################################################
all:





############################## PREAMBLE ####################################### 
# Essential source files to include.
SOURCES ?= $(abspath ./src)
aux := $(strip $(wildcard $(SOURCES)/*.mk))
include $(if $(aux), $(SOURCES)/*.mk, $(error Sources not found!))

# Ancillary configuration file (user's).
CONFIG_FILE ?= ./config.cfg
CONFIG_FILE := $(strip $(wildcard $(CONFIG_FILE)))
include $(if $(CONFIG_FILE), $(CONFIG_FILE), $(error Configuration file not found!))

# Define a CALL message for the log.
CALL = [$(shell /bin/date "+%Y-%m-%d(%H:%M:%S)") $(PIPELINE_NAME)]



############################## MAIN CODE ######################################

# Presentation header.
$(call presentation,                               $(PIPELINE_NAME), With upipe!)


# All processes complete.
all: $(pst_proc)
	$(info )
	$(info $(CALL) All processes complete.)

# Post processing stage.
$(pst_proc): $(mn_proc)

# Main processing stage.
$(mn_proc): $(pre_proc)

# Pre processing stage.
$(pre_proc): $(val_proc)

# Inputs' validation processing stage.
$(val_proc):

# Extra processes for further customizations.
$(ext_proc):


# Including targets definitions files.
include $(strip $(wildcard $(addprefix $(MODULES)/$(MODULE_NAME)/, $(TGT_DEFS))))
#include $(EXT_TGT)


.PHONY:# all pre_proc mn_proc pst_proc val_proc ext_proc
###############################################################################
