############################## HEADER #########################################

# Makefile data:
#
# 		Default goal: processSample
# 		File: processSample.mk
# 		Module: retrocopies
#		Project: upipe
# 		Author: Jose L. L. Buzzo
# 		Organization: RetroTeam
# 		Year: 2018


# Usage:
#
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

# Test the success of the validations process.
$(if $(INPUT_PROCESSED),, $(error "No valid input"))

# Carrefully set the INPUT_l variable after success validations.
INPUT_l := $(INPUT_PROCESSED)
INPUT_FILENAME_l := $(strip $(notdir $(INPUT_l)))
INPUT_BASENAME_l := $(strip $(basename $(notdir $(INPUT_l))))
INPUT_DIR_l := $(strip $(dir $(INPUT_l)))
#INPUT_DIR_BASENAME_l :=


#OUTPUT_l :=
#OUTPUT_FILENAME_l := $(INPUT_FILENAME_l)
#OUTPUT_BASENAME_l := $(INPUT_BASENAME_l)
OUTPUT_DIR := $(abspath $(strip $(OUTPUT_DIR)))# This is not a string list!

OUTPUT_DIR_BASENAME_l := $(addprefix $(OUTPUT_DIR)/, $(INPUT_BASENAME_l))
OUTPUT_l := $(strip $(join $(OUTPUT_DIR_BASENAME_l), $(addprefix /, $(INPUT_FILENAME_l))))
OUTPUT_ABNORMAL_l := $(OUTPUT_l:$(SUFFIXES)=.abnormal)

# List of '.abnormal' file names without its path.
ABNORMAL_l := $(INPUT_FILENAME_l:$(SUFFIXES)=.abnormal)

# ===> DEBUG CODE <===

# Enable some prints, if variable DBG="yes".
ifeq ($(DBG),yes)
$(info )
$(info INPUT_l: $(INPUT_l))
$(info INPUT_FILENAME_l: $(INPUT_FILENAME_l))
$(info INPUT_BASENAME_l: $(INPUT_BASENAME_l))
$(info INPUT_DIR: $(INPUT_DIR_l))
$(info )
$(info OUTPUT_DIR: A--$(OUTPUT_DIR)--A)
$(info OUTPUT_DIR_BASENAME_l: $(OUTPUT_DIR_BASENAME_l))
$(info OUTPUT_l: $(OUTPUT_l))
$(info OUTPUT_ABNORMAL_l: $(OUTPUT_ABNORMAL_l))
$(info ABNORMAL_l: $(ABNORMAL_l))
## The pattern must be extended to all, must be constant, not an array of
## different values. See:
$(info )
$(info AA314-$(patsubst /home/leonel/%$(SUFFIX), %.bam.o, $(INPUT)))
$(info AA315-$(patsubst $(OUTPUT_DIR)%.abnormal,%.abnormal.o, $(OUTPUT_ABNORMAL_l)))
$(info )
endif

# Enable emergy stop 2, if variable STP2="yes".
ifdef STP2
$(error Emergency stop 2.)
endif



############################## TARGETS ########################################

#################
#
# Main target: processSample
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

$(OUTPUT_ABNORMAL_l): %.abnormal: %.bai %.bam $(REF)/genes.bed
	$(info )
	$(info $(CALL) Selecting abnormal pairs from file $(word 2, $^))
	$(SCRIPTS)/findAbnormal.sh $(word 2, $^) $(word 3, $^) $(SEARCH_CRIT) $(@D) > $@



#################
#	
# Rule to index .bam file
#
#################

# Look for a *.sorted.bam.bai file, according to the link to 'OUTPUT_l'.
REQ_BAI = $(shell readlink -f $(filter %$(*F)$(SUFFIXES), $(OUTPUT_l)))
aux2 = $(strip $(filter %.sorted.bam, $(REQ_BAI)))
STBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*sorted*.bai')
NSTBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*.bai' -not -regex '.*sorted.*')
BAI = $(if $(aux2),$(STBAI),$(NSTBAI))


%.bai: %.bam
	$(info )
	$(info $(CALL) Indexing file $<)
	if [ -s "$(BAI)" ]; then \
		ln -sf "$(BAI)" $@; \
	else \
		samtools index -b -@ 8 $< $@; \
	fi



#################
#
# Rule to sort .bam file
#
#################

# Ancilary variables to this target.
REQ = $(filter %$(*F)$(SUFFIXES), $^)
CMD = $(shell samtools view -H $(REQ) 2> /dev/null | head -n1 | cut -f3)
BAM = $(shell readlink -f $(REQ))
STBAM = $(shell find $(dir $(BAM)) -type f -name '*$(*F)*.sorted.bam')

$(OUTPUT_l): %.bam: $(INPUT_l) | validations $(OUTPUT_DIR)
	$(info )
	$(info $(CALL) Creating link for file: $(REQ).)
	mkdir -p $(*D)
	if [ "$(CMD)" = "SO:coordinate" ]; then \
		ln -sf "$(BAM)" $@; \
	elif [ -s "$(STBAM)" ]; then \
		ln -sf "$(STBAM)" $@; \
	else \
		samtools sort -O BAM -m 8G -@ 8 $(REQ) -o $@; \
	fi
