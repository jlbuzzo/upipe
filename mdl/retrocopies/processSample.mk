###############################################################################
#
# processSample
#
###############################################################################

# Preamble.



#################
#
# MAin target: processSample
#
#################

processSample:: $(BAM_ABNORMAL_LIST)
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
	$(MDL)/findAbnormal.sh $(word 3, $^) $(word 2, $^) $(@D) > $@



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
REQ = $(filter %$*$(SUFFIXES), $^)
TGT = $(OUTPUT_DIR)/$*/$@ 
CMD = `samtools view -H $(REQ) 2> /dev/null | head -n1 | cut -f3`

$(notdir $(BAM_FINAL_LIST)): %.bam: $(BAM_LIST) | $(OUTPUT_DIR)
	$(info )
	$(info $(CALL) Created link for file: $(REQ).)
	mkdir -p $(OUTPUT_DIR)/$*
	if [ "$(CMD)" != "SO:coordinate" ]; then \
		samtools sort -O BAM -m 10G -@ 20 $(REQ) -o $(TGT); \
	else \
		ln -sf $(shell readlink -f $(REQ)) $(TGT); \
	fi

