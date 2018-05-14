###############################################################################
#
# Main Makefile of the module.
#
###############################################################################
main:



############################## PREAMBLE ####################################### 

# Some shortcuts.
RES := $(OUTPUT_DIR)/result
REF := $(OUTPUT_DIR)/reference
TEMP_PROCESS_DIR := $(OUTPUT_DIR)/temp/
DUMP_DIR := $(RES)/dump/
SCRIPTS := $(MDL)/scripts
CRIT :=

# Define a CALL message for the log.
CALL = [$(PIPELINE_NAME): $(shell date "+%Y-%m-%d(%H:%M:%S)")]


############################## TARGETS ########################################

main: mergeCall
	$(info Main done!)

count:
	$(info )
	$(info $(CALL) Counting retrocopies.)
	$(SCRIPTS)/counter.sh $(OUTPUT_DIR) > counts.txt



# Include target's specifications code.
include $(addprefix $(MDL)/, validations.mk processSample.mk mergeCall.mk commons.mk)


#.PHONY: validations
