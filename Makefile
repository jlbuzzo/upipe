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



############################## ENVIRONMENT TESTS ###############################

# Define switch to feel environment and take configuration file.
SWITCH			?= $(if $(shell [ -e ".ponga_switch" ] && echo "1"),1,)
CONFIG			?= config.mk
include			$(if $(shell [ -e "$(CONFIG)" ] && echo "1"),$(CONFIG),$(error Configuration file is missing))


# Export all variables or not.
export



############################## DEBUG ###########################################

# Debug code.
ifeq ($(DBG),yes)
$(info ############################## DEBUG ###########################################)
$(info )

# General variables.
$(info CONFIG:$(CONFIG).)
$(info TARGET0:$(TARGET0).)
$(info TARGET1:$(TARGET1).)
$(info )

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
target0: $(TARGET0)
	$(ECHO) "Validations done.\n"

# Extra target.
extra:
	$(ECHO) "Extra target done.\n"


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
