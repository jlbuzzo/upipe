#################
#
# Inout validations and variables treatmants.
#
#################

# Validate INPUT_DIR against emptyness.
INPUT_FILE := $(strip $(INPUT_FILE))
$(if $(INPUT_FILE),, $(error Variabe INPUT_FILE is undefined))


# Specifying and select INPUT_FILEs from a list with a suffix or from a directory.
$(if $(shell [ -f $(INPUT_FILE) ] && echo 1), $(eval BAM_LIST := $(wildcard $(filter %$(SUFFIX), $(shell cat $(INPUT_FILE))))), $(eval BAM_LIST := $(wildcard $(INPUT_FILE)/*$(SUFFIX))))
BAM_LIST := $(abspath $(strip $(BAM_LIST)))

# Verifyng .bam files presence in the list.
$(if $(BAM_LIST),, $(error No valid $(SUFFIX) file specified))

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

# Define a CALL message for the log.
CALL = [$(shell /bin/date "+%Y-%m-%d(%H:%M:%S)") $(PIPELINE_NAME)]

# A little information to the user.
$(info *******************************************************************************)
$(info 					BEGIN THE PROCESS!)
$(info *******************************************************************************)
$(info Part1: processSample)
$(info )



#################
#
# Main make target.
#
#################

processSample: $(BAM_ABNORMAL_LIST)
	$(info )
	$(info $(CALL) Finished processing for files: $^.)



#################
#
# Rule to find abnormal alignments
#
#################

$(BAM_ABNORMAL_LIST): %.abnormal: %.bai %.bam $(OUTPUT_DIR)/reference/genes.formated
	$(info )
	$(info $(CALL) Selecting abnormal pairs from file $(word 2, $^))
	$(SRC)/findAbnormal.sh $(@D) $(word 2, $^) $(word 3, $^) > $@



#################
#	
# Rule to index .bam file
#
#################		

%.bai: %.bam
	$(info )
	$(info $(CALL) Indexing file $<)
	samtools index -b -@ 20 $< $@ 



#################
#
# Rule to sort .bam file
#
#################

#$(info BAMS-$(BAM_FINAL_LIST))
$(BAM_FINAL_LIST): $(notdir $(BAM_FINAL_LIST)) 

# Ancilary variables to this goal.
REQ = $(filter %$*$(SUFFIX), $^)
TGT = $(OUTPUT_DIR)/$*/$@ 
CMD = `samtools view -H $(REQ) 2> /dev/null | head -n1 | cut -f3`

$(notdir $(BAM_FINAL_LIST)): %.bam: $(BAM_LIST) | $(OUTPUT_DIR)
	$(info )
	$(info $(CALL) Created link for file: $(REQ).)
	mkdir -p $(OUTPUT_DIR)/$*; \
	if [ "$(CMD)" != "SO:coordinate" ]; then \
		@echo 22222222k; \
		@samtools sort -O BAM -m 10G -@ 20 $(REQ) -o $(TGT); \
	elif [ -L "$(REQ)" ]; then \
		echo copiou!; \
		ln -sf $(shell readlink -f $(REQ)) $(TGT); \
	else \
		echo Não precisou copiar!; \
		ln -sf $(REQ) $(TGT); \
	fi



#################
#
# Rule to format bed file and select protein coding transcripts.
#
#################

$(OUTPUT_DIR)/reference/genes.formated: $(OUTPUT_DIR)/reference/genes.bed
	$(info )
	$(info $(CALL) Created protein coding genes bed at $<)
	grep protein_coding $< | sed 's/\t/\*\*/g' > $@



#################
#
# Rule to process .gtf
#
#################

$(OUTPUT_DIR)/reference/genes.bed: $(REF_ANNOTATION) | $(OUTPUT_DIR)/reference
	$(info )
	$(info $(CALL) Created genes.bed at $@)
	perl $(SRC)/gtf2bed.pl 'gene' $< > $@



#################
#
# Validate if $(SUFFIX) file is not corruped.
#
#################

validation: %.bam
	$(info )
	$(info $(CALL) Verifying integrity of file: $<.)
	samtools flagstat -@ 20 $<



#################
#
# Create output, reference and result directories.
#
#################

$(OUTPUT_DIR):
	$(info )
	$(info $(CALL) Created output dir $@)
	mkdir -p $@

$(OUTPUT_DIR)/reference:
	$(info )
	$(info $(CALL) Created reference dir $@)
	mkdir -p $@


.PHONY: validation
