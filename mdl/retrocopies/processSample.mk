############################## HEADER #########################################

# Makefile data:
# 		Main target: processSample
# 		File: processSample.mk
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

# ATTENTION: Variables that are string lists has an '_l' suffix appendedd to
# the end of their names. Their behavior differently in pattern rules.

# Test the succes of the validations process.
$(if $(INPUT_PROCESSED), , $(error "No valid input"))

# Carrefully set the INPUT_l variable after success validations.
INPUT_l := $(INPUT_PROCESSED)
INPUT_FILENAME_l := $(strip $(notdir $(INPUT_l)))
INPUT_BASENAME_l := $(strip $(basename $(notdir $(INPUT_l))))
INPUT_DIR_l := $(strip $(dir $(INPUT_l)))
#INPUT_DIR_BASENAME_l :=


#OUTPUT_l :=
#OUTPUT_FILENAME_l := $(INPUT_FILENAME_l)
#OUTPUT_BASENAME_l := $(INPUT_BASENAME_l)
OUTPUT_DIR := $(abspath $(OUTPUT_DIR))		# This is not a string list!

OUTPUT_DIR_BASENAME_l := $(addprefix $(OUTPUT_DIR)/, $(INPUT_BASENAME_l))
OUTPUT_l := $(strip $(join $(OUTPUT_DIR_BASENAME_l), $(addprefix /, $(INPUT_FILENAME_l))))
OUTPUT_ABNORMAL_l := $(OUTPUT_l:$(SUFFIXES)=.abnormal)

# List of '.abnormal' file names without its path.
ABNORMAL_l := $(INPUT_FILENAME:$(SUFFIX)=.abnormal)


# ===> DEBUG CODE <===
#$(info )
#$(info input: $(INPUT) -- input_dir: $(INPUT_DIR))
#$(info out_basename: $(OUTPUT_DIR_BASENAME))
#$(info )
#$(info abnormal: $(ABNORMAL))
#$(info outputs: $(OUTPUT))
#$(info )
# The pattern must be extended to all, must be constant, not an array of
# different values. See:
#$(info  AA314-$(patsubst /home/leonel/%$(SUFFIX), %.o, $(INPUT)))
#$(info  AA315-$(patsubst $(OUTPUT_DIR)%.abnormal, %.o, $(OUTPUT)))



############################## TARGETS ########################################

#################
#
# MAin target: processSample
#
#################

processSample: $(OUTPUT_ABNORMAL_l)
	$(info )
	$(info $(CALL) Target 'processSample' complete!)



#################
#
# Rule to find abnormal alignments
#
#################

$(OUTPUT_ABNORMAL_l): %.abnormal: %.bai %.bam $(OUTPUT_DIR)/reference/genes.bed
	$(info )
	$(info $(CALL) Selecting abnormal pairs from file $(word 2, $^))
	$(MDL)/findAbnormal.sh $(word 2, $^) $(word 3, $^) $(SEARCH_CRIT) $(@D) > $@



#################
#	
# Rule to index .bam file
#
#################		

%.bai: %.bam
	$(info )
	$(info $(CALL) Indexing file $<)
	if [ ! -s $@ ]; then \
		samtools index -b -@ 8 $< $@ \
	else \
		ln -sf $(shell readlink -f $(REQ)) $@; \
	fi



#################
#
# Rule to sort .bam file
#
#################

# Ancilary variables to this target.
REQ = $(filter %$(*F)$(SUFFIXES), $^)
CMD = `samtools view -H $(REQ) 2> /dev/null | head -n1 | cut -f3`

$(OUTPUT_l): %.bam: $(INPUT_l) | $(OUTPUT_DIR)
	$(info )
	$(info $(CALL) Creating link for file: $(REQ).)
	mkdir -p $*
	if [ "$(CMD)" != "SO:coordinate" ]; then \
		samtools sort -O BAM -m 8G -@ 8 $(REQ) -o $@; \
	else \
		ln -sf $(shell readlink -f $(REQ)) $@; \
	fi
