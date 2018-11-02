#!/usr/bin/make
############################## HEADER ##########################################
#
# This is a template of a generic Makefile.
#
################################################################################



############################## INFRASTRUCTURE ##################################
# This section must contain only general infrastructure variables.


# NOTE: Varable names can end with '_d', '_f', '_r', '_b', '_t' and 'l'.
# They represent, respectively, directories, files, references, booleans,
# temporaries and lists values inside them.
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 


# There could be an user for the container (ex. lion).
# So, C_BASE_d must be set accordingly:
USER ?= lion


# Host's base directory, outside the container: $(PWD).
H_BASE_d		?= $(PWD)
H_CONFIG_d		?= $(H_BASE_d)/config
H_INPUTS_d		?= $(H_BASE_d)/inputs
H_OUTPUTS_d		?= $(H_BASE_d)/outputs
H_ASSETS_d		?= $(H_BASE_d)/assets
H_REFERENCE_d	?= $(H_BASE_d)/reference
H_ANNOTATION_d	?= $(H_BASE_d)/annotation
H_EXTRA_d		?= $(H_BASE_d)/extra
H_TMP_d			?= $(H_BASE_d)/tmp


# Container's base directory: /home.
C_BASE_d		?= /home/$(USER)
C_CONFIG_d		?= $(C_BASE_d)/config
C_INPUTS_d		?= $(C_BASE_d)/inputs
C_OUTPUTS_d		?= $(C_BASE_d)/outputs
C_ASSETS_d		?= $(C_BASE_d)/assets
C_REFERENCE_d	?= $(C_BASE_d)/reference
C_ANNOTATION_d	?= $(C_BASE_d)/annotation
C_EXTRA_d		?= $(C_BASE_d)/extra
C_TMP_d			?= $(C_BASE_d)/tmp


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


# Pay attention here! This variables controls the entire Makefile.
INPUTS			:= $(INPUTS_d)
