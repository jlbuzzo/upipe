#!/usr/bin/make
################################################################################
#
# Make module for Git stuff.
#
################################################################################



############################## PREAMBLE ########################################



############################## CONFIG ##########################################

# NOTE: Varable names can end with '_d', '_f', '_r', '_b', '_t' and '_l', to
# represent, respectively, directories, files, references, booleans, temporaries
# and lists values in them.
#
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 

# Reference directory path.
PREFIX := $(CURDIR)


# Project's data.
PROJECT_TITLE := prj
PROJECT_VERSION := v0.1
PROJECT_DESCRIPTION := A project template.
PROJECT_AUTHORS := José Leonel L. Buzzo
PROJECT_ORGANIZATION := Galantelab
PROJECT_YEAR := $(shell date "+%Y")


# Git user data.
GIT_USERNAME := jlbuzzo
GIT_USER_PASSWD :=
GIT_REMOTE := https://github.com/$(GIT_USERNAME)/$(PROJECT).git 


# Git installation data.
# From distro repositories.
GIT_DEB_PACK :=
GIT_RPM_PACK := 
GIT_AUR_PACK := pacman -S git
# From tarball
GIT_GPG_URL :=
GIT_TARBALL_URL := https://mirrors.edge.kernel.org/pub/software/scm/git/git-2.19.0.tar.xz
GIT_TARBALL_SIGN_URL := https://mirrors.edge.kernel.org/pub/software/scm/git/git-2.19.0.tar.sign
# From inside environments.
GIT_PIP_PACK :=
GIT_CONDA_PACK :=


# Build variables.
REPO := $(PROJECT_TITLE)-$(PROJECT_VERSION)
BASE := $(PREFIX)/$(REPO)
README := $(BASE)/README.md

DOCKERFILE := $(BASE)/Dockerfile
DOCKERIGNORE := $(BASE)/.dockerignore
AUTOTOOLS = yes
GIT_INIT := yes
GITIGNORE := $(BASE)/.gitignore

DATE := $(shell date "+%Y-%m-%d")


############################## DEBUG ###########################################

# Debug code part. Print all variables already set.
ifeq ($(DBG),yes)
$(info PROJECT DATA)
$(info --------------------------------------------------------------------------------)
$(info PROJECT_TITLE: $(PROJECT_TITLE).)
$(info PROJECT_VERSION: $(PROJECT_VERSION).)
$(info PROJECT_DESCRIPTION: $(PROJECT_DESCRIPTION).)
$(info PROJECT_AUTHORS: $(PROJECT_AUTHORS).)
$(info PROJECT_ORGANIZATION: $(PROJECT_ORGANIZATION).)
$(info PROJECT_YEAR: $(PROJECT_YEAR).)
$(info )
$(info GIT INSTALLATION DATA)
$(info --------------------------------------------------------------------------------)
$(info GIT_DEB_PACK: $(GIT_DEB_PACK))
$(info GIT_RPM_PACK: $(GIT_RPM_PACK))
$(info GIT_AUR_PACK: $(GIT_AUR_PACK))
$(info GIT_GPG_URL: $(GIT_GPG_URL))
$(info GIT_TARBALL_URL: $(GIT_TARBALL_URL))
$(info GIT_PIP_PACK: $(GIT_PIP_PACK))
$(info GIT_CONDA_PACK: $(GIT_CONDA_PACK))
$(info )
$(info BUILD VARIABLES)
$(info --------------------------------------------------------------------------------)
$(info PREFIX: $(PREFIX).)
$(info BASE: $(BASE_d).)
$(info DATE: $(DATE).)
$(info )
endif



############################## TARGETS #########################################

all:


push_repo: repo
	@echo -e "Pushing ..."
	@echo -e "Repository: $(REPO)."
	@echo -e "Remote: $(GIT_REMOTE)."

repo: infra
	@echo -e "Make repository infrastructure:"

infra: folder



# Targets for the project's README.md file.
readme: $(README)

$(README): | $(BASE)
	echo "Creating README.md file."
	echo "### NAME ### \n\n\
	### VERSION ### \n\n\
	### PREAMBLE ### \n\n\
	### SYNOPSIS ### \n\n\
	" > $@


# Targets for the base project's directory.
folder: | $(BASE)

$(BASE):
	echo "Creating directory: $@."
	mkdir -p $@


# .PHONY targets.
.PHONY: readme folder
