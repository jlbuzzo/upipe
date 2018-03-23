###############################################################################
#
# Validations
#
###############################################################################

# Preamble.





###############################################################################
# Procedural.

$(info )
$(info *** VALIDATION PROCESSES ***)

# Validate INPUT_FILE against emptyness.
INPUT_FILE := $(strip $(INPUT_FILE))
$(if $(INPUT_FILE),, $(error Variabe INPUT_FILE is undefined))

# Specifying and select INPUT_FILEs from a list with a suffix or from a directory.
BAM_LIST := $(if $(shell [ -f $(INPUT_FILE) ] && echo 1), $(wildcard $(filter %$(SUFFIXES), $(shell cat $(INPUT_FILE)))), $(wildcard $(INPUT_FILE)/*$(SUFFIXES)))
BAM_LIST := $(abspath $(strip $(BAM_LIST)))
VERIFY = $(shell [ -s "$i" ] && echo "$i")
BAM_LIST := $(strip $(foreach i, $(BAM_LIST),$(VERIFY)))

# Verifyng .bam files presence in the list.
$(if $(BAM_LIST),, $(error No valid $(SUFFIXES) file specified))

#$(info $(BAM_LIST))
# Store INPUT_DIR list.
INPUT_DIR := $(abspath $(dir $(BAM_LIST)))
OUTPUT_DIR := $(abspath $(OUTPUT_DIR))

# Extract the basename of the files.
BAM_BASENAME_LIST	:= $(notdir $(BAM_LIST))
BAM_BASENAME_LIST	:= $(basename $(BAM_BASENAME_LIST))

# Join the final names list.
BAM_FINAL_LIST	:= $(join $(addsuffix /, $(BAM_BASENAME_LIST)), $(BAM_BASENAME_LIST))
BAM_FINAL_LIST	:= $(addprefix $(OUTPUT_DIR)/, $(addsuffix .bam, $(BAM_FINAL_LIST)))
BAM_FINAL_LIST	:= $(strip $(BAM_FINAL_LIST))

# Compose .abnormal names list.
BAM_ABNORMAL_LIST := $(patsubst %.bam, %.abnormal, $(BAM_FINAL_LIST))

# Define .log files.
LOG_FILE_LIST := $(addprefix $(OUTPUT_DIR)/, $(addsuffix .log, $(BAM_BASENAME_LIST)))

# This is in wrong place.
FILE_NAME := $(notdir $(INPUT_FILE))
ifneq ($(SAMPLE_NAME),NULL)
       SAMPLE_ID := $(FILE_NAME)
endif

LOG_FILE := $(OUTPUT_DIR)/$(SAMPLE_ID).log
timestamp = `/bin/date "+%Y-%m-%d(%H:%M:%S)"`


# Unic target.
validations:
	$(info )
	$(info Finished validations for files $^.)
	$(info *** VALIDATION COMPLETE ***)
