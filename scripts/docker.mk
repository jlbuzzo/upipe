#!/usr/bin/make
############################## HEADER ##########################################
# Makefile for docker stuff...
################################################################################





############################## PREAMBLE ########################################

# Command macros.
SHELL				:= bash
MODULE_NAME			:= docker
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


# Define switch to feel environment.
SWITCH			?= $(if $(shell [ -e ".ponga_switch" ] && echo "1"),1,)


# There could be an user for the container (ex. lion).
# So, C_BASE_d must be set accordingly:
DEFAULT_USER	?= lion


# Host's base directory, outside the container: $(PWD).
prefix			?= $(PWD)
H_BASE_d		?= $(prefix)
H_CONFIG_d		?= $(H_BASE_d)/config
H_INPUTS_d		?= $(H_BASE_d)/inputs
H_OUTPUTS_d		?= $(H_BASE_d)/outputs
H_ASSETS_d		?= $(H_BASE_d)/assets
H_REFERENCE_d	?= $(H_BASE_d)/reference
H_ANNOTATION_d	?= $(H_BASE_d)/annotation
H_EXTRA_d		?= $(H_BASE_d)/extra
H_TMP_d			?= $(H_BASE_d)/tmp


# Container's base directory: /home/$(DEFAULT_USER).
C_BASE_d		?= /home/$(DEFAULT_USER)
C_CONFIG_d		?= $(C_BASE_d)/config
C_INPUTS_d		?= $(C_BASE_d)/inputs
C_OUTPUTS_d		?= $(C_BASE_d)/outputs
C_ASSETS_d		?= $(C_BASE_d)/assets
C_REFERENCE_d	?= $(C_BASE_d)/reference
C_ANNOTATION_d	?= $(C_BASE_d)/annotation
C_EXTRA_d		?= $(C_BASE_d)/extra
C_TMP_d			?= $(C_BASE_d)/tmp


# Locale variables.
# A common internal representation. The best place to overwrite variables.
BASE_d			:= $(if $(SWITCH),$(C_BASE_d),$(H_BASE_d))
CONFIG_d		:= $(if $(SWITCH),$(C_CONFIG_d),$(H_CONFIG_d))
INPUTS_d		:= $(if $(SWITCH),$(C_INPUTS_d),$(H_INPUTS_d))
OUTPUTS_d		:= $(if $(SWITCH),$(C_OUTPUTS_d),$(H_OUTPUTS_d))
ASSETS_d		:= $(if $(SWITCH),$(C_ASSETS_d),$(H_ASSETS_d))
REFERENCE_d		:= $(if $(SWITCH),$(C_REFERENCE_d),$(H_REFERENCE_d))
ANNOTATION_d	:= $(if $(SWITCH),$(C_ANNOTATION_d),$(H_ANNOTATION_d))
EXTRA_d			:= $(if $(SWITCH),$(C_EXTRA_d),$(H_EXTRA_d))
TMP_d			:= $(if $(SWITCH),$(C_TMP_d),$(H_TMP_d))

# All directories.
DIRS_dl			:= $(BASE_d) $(CONFIG_d) $(INPUTS_d) $(OUTPUTS_d) $(ASSETS_d) $(REFERENCE_d) $(ANNOTATION_d) $(EXTRA_d) $(TMP_d)


# Docker specific variables. Tarball v18.06.
CONFIG				?= $(BASE_d)/config.mk
IMAGE				?=
DOCKERFILE			?= $(BASE_d)/Dockerfile
DOCKER_TARBALL		?= $(TMP_d)/docker-18.06.1-ce.tgz
DOCKER_TARBALL_URL	?= https://download.docker.com/linux/static/stable/x86_64/docker-18.06.1-ce.tgz


# Mount some important files in their respective folders. Must use absolute paths!
MP_CONFIG_l		?= $(PWD)/Makefile $(PWD)/config.mk
MP_INPUTS_l		?= $(INPUTS)
MP_OUTPUTS_l	?=
MP_ASSETS_l		?= /home/scratch60/lbuzzo/RTC/neopipe/assets/ref.perldb
MP_REFERENCE_l	?= /home/genomes/Homo_sapiens/hg38/hg38.fa
MP_ANNOTATION_l	?= /home/projects2/databases/gencode/release29/gencode.v29.annotation.gff3.gz
MP_EXTRA_l		?=
MP_TMP_l		?=

# Map files.
STEM_l			:= CONFIG INPUTS OUTPUTS ASSETS REFERENCE ANNOTATION EXTRA TMP



############################## ENVIRONMENT TESTS ###############################

# Essential arguments pre-validations.
CONFIG				?= config.mk
include $(if $(shell [ -e "$(CONFIG)" ] && echo "1"),$(CONFIG),$(error "Configuration file is missing"))
$(if $(strip $(IMAGE)),, $(error "Docker image undefined!"))


# Conditional commands.
HAS_INTERNET		:= read <<< "$$($(PING) www.google.com)" && echo -n $$REPLY
GET_INTERNET		:= $(if $(shell $(HAS_INTERNET)),,internet_configure)

HAS_DOCKER			:= if [ -n "$$(which docker)" ]; then echo "1"; fi
GET_DOCKER			:= $(if $(shell $(HAS_DOCKER)),,docker_install)

HAS_DOCKERFILE		:= if [ -f "$(DOCKERFILE)" ]; then echo "1"; fi
GET_MODE			:= $(if $(shell $(HAS_DOCKERFILE)),docker_build,docker_pull)

HAS_IMAGE			:= echo -n "$$(sed -n '\:^$(IMAGE) :{p}' <<< "$$(docker images 2> /dev/null)")"
GET_IMAGE			:= $(if $(shell $(HAS_IMAGE)),,$(GET_MODE))


# Export all variables or not.
export



############################## DEBUG ###########################################

# Debug code.
ifeq ($(DBG),yes)
$(info "############################## DEBUG ###########################################")
$(info )

# General variables.
$(info PIPELINE:$(PIPELINE).)
$(info CURDIR:$(CURDIR).)
$(info PWD:$(PWD).)
$(info MAKEFILE_LIST:$(MAKEFILE_LIST).)
$(info CONFIG:$(CONFIG).)
$(info SWITCH:$(SWITCH).)
$(info )

# Docker variables.
$(info IMAGE:$(IMAGE).)
$(info DOCKERFILE:$(DOCKERFILE).)
$(info DOCKER_TARBALL:$(DOCKER_TARBALL).)
$(info DOCKER_TARBALL_URL:$(DOCKER_TARBALL_URL).)
$(info )

# Operational variables.
$(info HAS_INTERNET:$(HAS_INTERNET).)
$(info GET_INTERNET:$(GET_INTERNET).)
$(info HAS_DOCKER:$(HAS_DOCKER).)
$(info GET_DOCKER:$(GET_DOCKER).)
$(info HAS_DOCKER_FILE:$(HAS_DOCKER_FILE).)
$(info GET_MODE:$(GET_MODE).)
$(info HAS_IMAGE:$(HAS_IMAGE).)
$(info GET_IMAGE:$(GET_IMAGE).)
$(info )

# Infrastructure variables.
$(info BASE_d:$(BASE_d).)
$(info CONFIG_d:$(CONFIG_d).)
$(info INPUTS_d:$(INPUTS_d).)
$(info OUTPUTS_d:$(OUTPUTS_d).)
$(info ASSETS_d:$(ASSETS_d).)
$(info REFERENCE_d:$(REFERENCE_d).)
$(info ANNOTATION_d:$(ANNOTATION_d).)
$(info EXTRA_d:$(EXTRA_d).)
$(info TMP_d:$(TMP_d).)
$(info )

# Infrastructure files.
$(info MP_CONFIG_l:$(MP_CONFIG_l).)
$(info MP_INPUTS_l:$(MP_INPUTS_l).)
$(info MP_OUTPUTS_l:$(MP_OUTPUTS_l).)
$(info MP_ASSETS_l:$(MP_ASSETS_l).)
$(info MP_REFERENCE_l:$(MP_REFERENCE_l).)
$(info MP_ANNOTATION_l:$(MP_ANNOTATION_l).)
$(info MP_EXTRA_l:$(MP_EXTRA_l).)
$(info MP_TMP_l:$(MP_TMP_l).)
$(info )

$(info STEM_l:$(STEM_l).)
$(info MAPS_l:$(MAPS_l).)
$(info )

$(info "################################################################################")
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
	$(ECHO) "\tmake -f my_Makefile IMAGE=my_image ARGS=\"arg1 arg2\" [options] all"


# This is a default testing target.
simple_test:
	$(ECHO) "This is a simple test for arguments:$(ARGS)."



# Routine targets.
all: docker_run
	$(ECHO) "All done.\n"


# Run in Docker.
docker_run: has_docker_image $(DIRS_dl)#$(CONFIG)
	$(ECHO) "Run docker image $(IMAGE)..."
	$(DOCKER_RUN) \
		-w $(C_BASE_d) \
		-v $(H_CONFIG_d):$(C_CONFIG_d) \
		-v $(H_INPUTS_d):$(C_INPUTS_d) \
		-v $(H_OUTPUTS_d):$(C_OUTPUTS_d) \
		-v $(H_ASSETS_d):$(C_ASSETS_d) \
		-v $(H_REFERENCE_d):$(C_REFERENCE_d) \
		-v $(H_ANNOTATION_d):$(C_ANNOTATION_d) \
		-v $(H_EXTRA_d):$(C_EXTRA_d) \
		-v $(H_TMP_d):$(C_TMP_d) \
		$(MAPS_l) $(IMAGE) $(ARGS)
	$(ECHO) "Done.\n"


has_docker_image: $(GET_IMAGE) # docker_pull or docker_build
	$(ECHO) "This system has the Docker image: $(IMAGE).\n"


docker_pull: has_docker | has_internet
	$(ECHO) "Pulling Docker image: $(IMAGE)."
	docker pull $(IMAGE)
	$(ECHO) "Done.\n"


docker_build: has_docker $(DOCKERFILE) | has_internet
	$(ECHO) "Building Docker image: $(IMAGE)."
	docker build -t $(IMAGE) .
	$(ECHO) "Done.\n"


has_docker: $(GET_DOCKER) # docker_install
	$(ECHO) "This system has Docker.\n"



############################## INSTALLING ######################################

# Base installation targets.
docker_install: $(DOCKER_TARBALL)
	$(ECHO) "Installing Docker...."
	#cd $(DOCKER_TARBALL)
	#./configure
	#$(MAKE) && $(MAKE) install
	sudo cp $(subst .tgz,,$(DOCKER_TARBALL))/* /usr/bin/
	sudo dockerd &
	$(ECHO) "Done.\n"


$(DOCKER_TARBALL): | has_internet $(TMP_d)
	$(ECHO) "Downloading and unpacking Docker tarball..."
	wget $(DOCKER_TARBALL_URL) -P $@
	tar xzvf $@ -C $(TMP_d)
	$(ECHO) "Done.\n"



# Common targets.
has_internet: $(GET_INTERNET) # internet_configure
	$(ECHO) "This system has internet connection.\n"


internet_configure:
	$(ECHO) "Trying to configure system's internet..."
	$(ECHO) "Error: You need to have internet access to proceed."
	#systemctl ...
	$(ECHO) "Done.\n"


$(DIRS_dl):
	$(ECHO) "Creating directory $@."
	$(MKDIR) $@
	$(ECHO) "Done.\n"




# Obligatory targets.
.PHONY:
