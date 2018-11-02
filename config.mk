############################## HEADER ##########################################
#
# This is a template of a configuration file.
#
################################################################################



############################## INFRASTRUCTURE ##################################
# This section must contain only general infrastructure variables.


# NOTE: Varable names can end with '_d', '_f', '_r', '_t' and '_l'.
# They represent, respectively, directories, files, references, temporaries and
# lists values inside them.
# Only names ending with '_l' can have non-single string values, but you can
# compose'em: OUTPUTS_dlt (for a list of temporary directories). 


# There could be an user for the container (ex. lion).
# So, C_BASE_d must be set accordingly:
DEFAULT_USER	:= lion


# Host's base directory, outside the container: $(PWD).
prefix :=		$(PWD)#$(PWD)/pocket
H_BASE_d		:= $(prefix)
H_CONFIG_d		:= $(H_BASE_d)/config
H_INPUTS_d		:= $(H_BASE_d)/inputs
H_OUTPUTS_d		:= $(H_BASE_d)/outputs
H_ASSETS_d		:= $(H_BASE_d)/assets
H_REFERENCE_d	:= $(H_BASE_d)/reference
H_ANNOTATION_d	:= $(H_BASE_d)/annotation
H_EXTRA_d		:= $(H_BASE_d)/extra
H_TMP_d			:= $(H_BASE_d)/tmp


# Container's base directory: /home.
C_BASE_d		:= /home/$(DEFAULT_USER)
C_CONFIG_d		:= $(C_BASE_d)/config
C_INPUTS_d		:= $(C_BASE_d)/inputs
C_OUTPUTS_d		:= $(C_BASE_d)/outputs
C_ASSETS_d		:= $(C_BASE_d)/assets
C_REFERENCE_d	:= $(C_BASE_d)/reference
C_ANNOTATION_d	:= $(C_BASE_d)/annotation
C_EXTRA_d		:= $(C_BASE_d)/extra
C_TMP_d			:= $(C_BASE_d)/tmp


# Mount some important files in their respective folders. Must use absolute paths!
MP_CONFIG_l		:= $(PWD)/Makefile $(PWD)/config.mk
MP_INPUTS_l		:= $(INPUTS)
MP_OUTPUTS_l	:=
MP_ASSETS_l		:= /home/scratch60/lbuzzo/RTC/neopipe/assets/ref.perldb
MP_REFERENCE_l	:= /home/genomes/Homo_sapiens/hg38/hg38.fa
MP_ANNOTATION_l	:= /home/projects2/databases/gencode/release29/gencode.v29.annotation.gff3.gz
MP_EXTRA_l		:=
MP_TMP_l		:=



############################## APP SPECIFIC VARIABLES ##########################
# This section must contain variables specific to the application.


# General pipeline data: name and file suffixes.
PIPELINE 			:= CUSTOM
SUFFIXES 			:=



# Meta-targets.
TARGET3				:=
TARGET3_REQUISITES	:=
TARGET3_RECIPES		:=

TARGET2				:=
TARGET2_REQUISITES	:=
TARGET2_RECIPES		:=

TARGET1				:=
TARGET1_REQUISITES	:=
TARGET1_RECIPES		:=

TARGET0				:=
TARGET0_REQUISITES	:=
TARGET0_RECIPES		:=
