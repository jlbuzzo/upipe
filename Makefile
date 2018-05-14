###############################################################################
#
# A Makefile to rule'em all!
#
###############################################################################
all:





############################## PREAMBLE ####################################### 

# Include essential source files.
SOURCES ?= $(abspath ./src)
aux := $(abspath $(strip $(wildcard $(SOURCES)/*.mk)))
include $(if $(aux), $(aux), $(error Sources not found))

# Include user's configuration file (only the first one).
CONFIG_FILE ?= ./config.cfg
aux := $(word 1, $(abspath $(strip $(wildcard $(CONFIG_FILE)))))
include $(if $(aux), $(aux), $(error Configuration file not found))

# Set modules' path (only one module).
MDL := $(abspath $(strip $(wildcard $(addprefix $(MODULES)/, $(MODULE_NAME)))))
$(if $(MDL),, $(error Module '$(notdir $(MDL))' not found))

# Define a CALL message for the log.
CALL = [$(UPIPE): $(shell date "+%Y-%m-%d(%H:%M:%S)")]

# Export all variables.
export
# Don't export duplicated presentation message.
unexport presentation



############################## TARGETS ########################################

# Presentation header.
$(call presentation,                               $(UPIPE), Running: $(PIPELINE_NAME)!)


# All processes complete.
all: $(pst_proc)
	$(info )
	$(info $(CALL) All processes complete!)

# Post process.
$(pst_proc): $(mn_proc)
	$(info )
	$(info $(CALL) Post process: '$(pst_proc)')
ifneq ($(strip $(pst_proc_tgt)),)
	$(MAKE) $(pst_proc_tgt)
endif

# Main process.
$(mn_proc): $(pre_proc)
	$(info )
	$(info $(CALL) Main process: '$(mn_proc)')
ifneq ($(strip $(mn_proc_tgt)),)
	$(MAKE) $(mn_proc_tgt)
endif

# Pre process.
$(pre_proc): $(val_proc)
	$(info )
	$(info $(CALL) Pre process: '$(pre_proc)')
ifneq ($(strip $(pre_proc_tgt)),)
	$(MAKE) $(pre_proc_tgt)
endif

# Validation process.
$(val_proc):
	$(info )
	$(info $(CALL) Validation process: '$(val_proc)')
ifneq ($(strip $(val_proc_tgt)),)
	$(MAKE) $(val_proc_tgt)
endif

# Extra process for further customizations.
$(ext_proc):
	$(info )
	$(info Extra process: '$(ext_proc)')
ifneq ($(strip $(ext_proc_tgt)),)
	$(MAKE) $(ext_proc_tgt)
endif



# Include extra targets' definition files.
aux := $(abspath $(strip $(wildcard $(TGT_EXT))))
ifneq ($(aux),)
$(info $(CALL) Included extra files: $(aux).)
include $(aux)
else
$(info $(CALL) No extra files included.)
endif



# .PHONY targets.
.PHONY: $(TGT_PHONY)
