###############################################################################
#
# A Makefile to rule'em all!
#
###############################################################################





#################
#
# Rule to format bed file and select protein coding transcripts
#
#################

$(REF)/genes.formated: $(REF)/genes.bed | $(REF)
	$(info )
	$(info $(CALL) Make genes.formated.)
	grep 'protein_coding' $< | sed 's/\t/\*\*/g' > $@
	@echo "$(CALL): Created protein coding genes bed from: $<.\n"
 
#################
#
# Rule to process .gtf
#
#################

$(REF)/genes.bed: $(REF_ANNOTATION) | $(REF)
	$(info )
	$(info Make $(CALL) genes.bed.)
	perl $(MDL)/gtf2bed.pl 'gene' $< > $@
	@echo "$(CALL): Created genes bed at: $@.\n"



$(REF)/exons.bed: $(REF_ANNOTATION) | $(REF)
	$(info )
	$(info Make $(CALL) exons.bed.)
	perl $(MDL)/gtf2bed.pl 'exon' $< > $@
	@echo "$(CALL): Created exons bed at: $@.\n"



#################
#
# Create results directory
#
#################

$(OUTPUT_DIR):
	$(info )
	$(info $(CALL) Make output dir.)
	mkdir -p $@

$(REF):
	$(info )
	$(info $(CALL) Make reference dir.)
	mkdir -p $@

$(RES):
	$(info )
	$(info $(CALL) Make result dir.)
	mkdir -p $@

$(TEMP_PROCESS_DIR):
	$(info )
	$(info $(CALL) Make tmp dir.)
	mkdir -p $@

