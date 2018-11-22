#!/usr/bin/make
############################## HEADER ##########################################
# Makefile for general stuff...
################################################################################





############################## PREAMBLE ########################################

# Command macros.
SHELL				:= bash
MODULE_NAME			:= main
#CALL				:= [$(MODULE_NAME): $(shell date --utc)]
ECHO				:= echo -e
MKDIR				:= mkdir -p
PING				:= ping -c1
DOCKER_RUN			:= docker run --rm -ti



############################## INFRASTRUCTURE ##################################
# This section must contain only general infrastructure variables.


# NOTE: Varable names can end with '_d', '_f', '_r', '_t' and '_l'.
# They represent, respectively, directories, files, references, temporaries and
# lists values inside them.
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 


# Include user's configuration file (only the first one).
SWITCH				:= $(shell [ -f ".ponga_switch" ] && echo "1")
CONFIG_f			?= $(if $(SWITCH),/home/config/config.mk,config.mk)
include				$(if $(wildcard $(CONFIG_f)),$(CONFIG_f),$(error Configuration file is missing))


# There could be an user for the container (ex. lion).
# So, C_BASE_d must be set accordingly:
DEFAULT_USER		?= $(if $(SWITCH),lion,)
POCKET				?= $(if $(prefix),$(prefix),$(POCKET))


# Host's base directory, outside the container: $(PWD).
H_BASE_d			?= $(abspath $(POCKET))
H_CONFIG_d			?= $(H_BASE_d)/config
H_INPUTS_d			?= $(H_BASE_d)/inputs
H_OUTPUTS_d			?= $(H_BASE_d)/outputs
H_ASSETS_d			?= $(H_BASE_d)/assets
H_REFERENCE_d		?= $(H_BASE_d)/reference
H_ANNOTATION_d		?= $(H_BASE_d)/annotation
H_EXTRA_d			?= $(H_BASE_d)/extra
H_TMP_d				?= $(H_BASE_d)/tmp


# Container's base directory: /home.
C_BASE_d			?= /home/$(DEFAULT_USER)
C_CONFIG_d			?= $(C_BASE_d)/config
C_INPUTS_d			?= $(C_BASE_d)/inputs
C_OUTPUTS_d			?= $(C_BASE_d)/outputs
C_ASSETS_d			?= $(C_BASE_d)/assets
C_REFERENCE_d		?= $(C_BASE_d)/reference
C_ANNOTATION_d		?= $(C_BASE_d)/annotation
C_EXTRA_d			?= $(C_BASE_d)/extra
C_TMP_d				?= $(C_BASE_d)/tmp


# Mount some important files in their respective folders. Must use ionly absolute paths!
MP_CONFIG_l			?=
MP_INPUTS_l			?=
MP_OUTPUTS_l		?=
MP_ASSETS_l			?=
MP_REFERENCE_l		?=
MP_ANNOTATION_l		?=
MP_EXTRA_l			?=
MP_TMP_l			?=


# A common internal representation. The best place to overwrite variables.
BASE_d				:= $(if $(SWITCH),$(C_BASE_d),$(H_BASE_d))
CONFIG_d			:= $(if $(SWITCH),$(C_CONFIG_d),$(H_CONFIG_d))
INPUTS_d			:= $(if $(SWITCH),$(C_INPUTS_d),$(H_INPUTS_d))
OUTPUTS_d			:= $(if $(SWITCH),$(C_OUTPUTS_d),$(H_OUTPUTS_d))
ASSETS_d			:= $(if $(SWITCH),$(C_ASSETS_d),$(H_ASSETS_d))
REFERENCE_d			:= $(if $(SWITCH),$(C_REFERENCE_d),$(H_REFERENCE_d))
ANNOTATION_d		:= $(if $(SWITCH),$(C_ANNOTATION_d),$(H_ANNOTATION_d))
EXTRA_d				:= $(if $(SWITCH),$(C_EXTRA_d),$(H_EXTRA_d))
TMP_d				:= $(if $(SWITCH),$(C_TMP_d),$(H_TMP_d))


# Infrastructure directories.
INFRASTRUCTURE_dl	:= $(CONFIG_d) $(INPUTS_d) $(OUTPUTS_d) $(ASSETS_d) $(REFERENCE_d) $(ANNOTATION_d) $(EXTRA_d) $(TMP_d)


# Export all variables or not.
export



############################## VALIDATIONS #####################################

############################## DEBUG ###########################################

# Debug code.
ifeq ($(DBG),yes)
$(info ############################## DEBUG ###########################################)

# General variables.
$(info PIPELINE:$(PIPELINE).)
$(info SUFFIXES:$(SUFFIXES).)
$(info SWITCH:$(SWITCH).)
$(info CONFIG_f:$(CONFIG_f).)
$(info CURDIR:$(CURDIR).)
$(info POCKET:$(POCKET).)
$(info )

$(info H_BASE_d:$(H_BASE_d).)
$(info INFRASTRUCTURE_dl:$(INFRASTRUCTURE_dl).)
$(info )

$(info TARGET0:$(TARGET0).)
$(info TARGET1:$(TARGET1).)
$(info TARGET2:$(TARGET2).)
$(info TARGET3:$(TARGET3).)

$(info ################################################################################)
$(info )
endif


# Emergency estop.
ifeq ($(STP),yes)
$(error Emergency stop)
endif



############################## TARGETS #########################################

# Help must be the first target, the default goal.
help: 
	$(ECHO) "Usage:\n"
	$(ECHO) "\tmake all [options]\n"
	$(ECHO) "\tOptions:"
	$(ECHO) "\t\t-f\tSpecify a Makefile with another name."
	$(ECHO) "\t\t-s\tSilent mode."
	$(ECHO) "\t\t-j\tPropagate parallelism to sub-makes. This is serial only."
	$(ECHO) "\t\tARGS\tSpecify extra arguments (ex. ARGS=\"arg1 arg2 ...\")."
	$(ECHO) "\t\tCONFIG\tSpecify configuration file to parse (ex. CONFIG=\"my_config.mk\")."
	$(ECHO) "\t\tDBG\tDebug messages (ex. DBG=yes)."
	$(ECHO) "\t\tSTP\tEmergency stop (ex. STP=yes)."


# This is a default testing target.
simple_test:
	$(ECHO) "This is a simple test for arguments:$(ARGS)."



# All Makefile is serial, but recursive the ones can be parallelized.
.NOTPARALLEL:

# Main targetis chain.
all: target3
	$(ECHO) "All done.\n"

# Postprocess
target3: target2 $(TARGET3)
	$(ECHO) "Postprocess done.\n"

# Process.
target2: target1 $(TARGET2)
	$(ECHO) "process done.\n"

# Preprocess.
target1: target0 $(TARGET1)
	$(ECHO) "Preprocess done.\n"

# Validations.
target0: $(TARGET0) | $(INFRASTRUCTURE_dl)
	$(ECHO) "Validations done.\n"

# Extra target.
extra:
	$(ECHO) "Extra target done.\n"

$(INFRASTRUCTURE_dl):
	$(ECHO) "Createing directory $@."
	mkdir -p $@

clean:
	$(ECHO) "Removing infrastructure directories:\n$(INFRASTRUCTURE_dl)"
	rm -rf $(INFRASTRUCTURE_dl)


# User targets chain.
$(TARGET3): $(TARGET3_REQUISITES)
	$(TARGET3_RECIPES)
	$(ECHO) "$@ done."

$(TARGET2): $(TARGET2_REQUISITES)
	$(TARGET2_RECIPES)
	$(ECHO) "$@ done."

$(TARGET1): $(TARGET1_REQUISITES)
	$(TARGET1_RECIPES)
	$(ECHO) "$@ done."

$(TARGET0): $(TARGET0_REQUISITES)
	$(TARGET0_RECIPES)
	$(ECHO) "$@ done."



# Obligatory targets.
.PHONY:
