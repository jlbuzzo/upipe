###############################################################################
#
# Main Makefile of the module.
#
###############################################################################



############################## PREAMBLE ####################################### 

# Essential source files to include.
#SOURCES ?= $(abspath ./src)
#aux := $(strip $(wildcard $(SOURCES)/*.mk))
#include $(if $(aux), $(SOURCES)/*.mk, $(error Sources not found!))

# Include the variable definitions file. Can't overwrite it.
CONFIG_FILE ?= ./config.conf
CONFIG := $(strip $(wildcard $(CONFIG_FILE)))
include $(if $(CONFIG), $(CONFIG), $(error Configuration file not found))

# Define a CALL message for the log.
CALL = [$(shell date "+%Y-%m-%d(%H:%M:%S)") $(PIPELINE_NAME)]

# Shortcut to modules' path and scripts.
MODULES := $(abspath ../)# Dangerous path: to much specific.
MDL := $(MODULES)/$(MODULE_NAME)
SCRIPTS := $(MDL)/scripts



############################## TARGETS ########################################

all: mergeCall

mergeCall: processSample

processSample: validations

validations:



# Include target's specifications code.
include validations.mk processSample.mk mergeCall.mk commons.mk dscasd.mk
