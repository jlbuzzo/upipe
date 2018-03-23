###############################################################################
#
# A Makefile to rule'em all!
#
###############################################################################


pth := $(MODULES)/$(MODULE_NAME)


#################
#
# Rule to format bed file and select protein coding transcripts
#
#################

$(OUTPUT_DIR)/reference/genes.formated: $(OUTPUT_DIR)/reference/genes.bed | $(OUTPUT_DIR)/reference
	$(info )
	$(info Make genes.formated.)
	cat $(OUTPUT_DIR)/reference/genes.bed | grep protein_coding | sed 's/\t/\*\*/g' > $(OUTPUT_DIR)/reference/genes.formated
	@echo "$(timestamp) $(PIPELINE_NAME): Created protein coding genes bed at: $(OUTPUT_DIR)/reference/exons.bed\n" >> $(LOG_FILE)
 
#################
#
# Rule to process .gtf
#
#################

$(OUTPUT_DIR)/reference/genes.bed: $(REF_ANNOTATION) | $(OUTPUT_DIR)/reference
	$(info )
	$(info Make genes.bed.)
	perl $(pth)/gtf2bed.pl 'gene' $< > $@
	@echo "$(timestamp) $(PIPELINE_NAME): Created genes bed at: $(OUTPUT_DIR)/reference/genes.bed\n" >> $(LOG_FILE)



$(OUTPUT_DIR)/reference/exons.bed: $(REF_ANNOTATION) | $(OUTPUT_DIR)/reference
	$(info )
	$(info Make exons.bed.)
	perl $(pth)/gtf2bed.pl 'exon' $< > $@
	@echo "$(timestamp) $(PIPELINE_NAME): Created exons bed at: $(OUTPUT_DIR)/reference/exons.bed\n" >> $(LOG_FILE)



#################
#
# Create results directory
#
#################

$(OUTPUT_DIR):
	$(info )
	$(info Make output dir.)
	mkdir -p $(OUTPUT_DIR)
	@echo "$(timestamp) $(PIPELINE_NAME): Created output dir: $(OUTPUT_DIR).\n"

$(OUTPUT_DIR)/reference:
	$(info )
	$(info Make reference dir.)
	mkdir -p $(OUTPUT_DIR)/reference
	@echo "$(timestamp) $(PIPELINE_NAME): Created reference dir: $(OUTPUT_DIR)/reference.\n"

$(OUTPUT_DIR)/result:
	$(info )
	$(info Make result dir.)
	mkdir -p $(OUTPUT_DIR)/result
	@echo "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/result.\n"

$(TEMP_PROCESS_DIR):
	$(info )
	$(info Make tmp dir.)
	mkdir -p $(TEMP_PROCESS_DIR)

