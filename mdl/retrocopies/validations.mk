############################## HEADER #########################################

# Makefile data:
# 		Main target: validations
# 		File: validations.mk
# 		Module: retrocopies
#		Project: upipe
# 		Author: Jose L. L. Buzzo
# 		Organization: RetroTeam
# 		Year: 2018
#
#
# Usage:
# 		Input args: INPUT variable (string list).
# 					The string list must contain, at least, a path to a folder
# 					or a file in the system.
#
# 		Output args: PROCESSED_INPUT variable (string list).
# 					The returned string list caontains all the existent files
# 					specified in the input args with canonical paths and
# 					filtered by the SUFFIXES auxilliary variable.
#
#		Aux args: [SUFFIXES] variable (string list).
#					A list of patterns to filter or search for files in the
#					input args. It's optional.



############################## PREAMBLE #######################################
# procedural pre-processing part.

# Validate INPUT against emptyness.
INPUT_PROCESSED := $(strip $(INPUT))
$(if $(INPUT_PROCESSED),, $(error Variabe INPUT is empty))

# Search and validate INPUT from a given string list, directory, or file, by
# the SUFFIXES list criteria.
#INPUT_PROCESSED := $(shell $(SCRIPTS)/ufinder.sh $(INPUT_PROCESSED) $(SUFFIXES))

# A limited alternative validation process (without 'ufinder' script).
INPUT_PROCESSED := $(if $(shell [ -f $(INPUT_PROCESSED) ] && echo 1), $(wildcard $(filter %$(SUFFIXES), $(shell cat $(INPUT_PROCESSED)))), $(wildcard $(INPUT_PROCESSED)/*$(SUFFIXES)))
INPUT_PROCESSED := $(strip $(abspath $(INPUT_PROCESSED)))

# Verifyng SUFFIXES terminated files' remaining after the filtering of the
# INPUT_PROCESSED list.
$(if $(INPUT_PROCESSED),, $(error No valid $(SUFFIXES) file specified))
$(if $(wildcard $(OUTPUT_DIR)),, $(error Output directory '$(OUTPUT_DIR)' doesnt exist))

# ===> DEBUG CODE <===
#$(error Stop emergency test..)



############################## TARGETS ########################################

validations:
	$(info )
	$(info Target 'validations' complete!)
