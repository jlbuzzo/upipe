# Preamble.

OUTPUT_DIR := save
RES := $(OUTPUT_DIR)/result
$(info $(RES))


#################
#
# mergeCall part.
# 
#################

mergeCall: $(RES)/putativeins.min.norep.exonic.dist.notsimilar.orientation
	$(info )
	$(info $(CALL) Finished processing for files: $^.)



#################
#
# Extract reads from bam files
# 
#################

$(RES)/putativeins.min.norep.exonic.dist.notsimilar.orientation: $(RES)/putativeins.min.norep.exonic.dist.notsimilar
	$(info )
	$(info $(CALL) Recalculating Insertion Point and Support from original BAM files; Checking supporting reads orientation.)
	$(SRC)/extract_reads.sh $@ @< $(RES)



#################
#
# Rule to select non similar IP/PARENTAL sequences 
#
#################

$(RES)/putativeins.min.norep.exonic.dist.notsimilar: $(RES)/putativeins.min.norep.exonic.dist
	$(info )
	$(info $(CALL) Removing clusters which Insertion Point is simmilar to Parental Sequence.)
	$(SRC)/select_non_similar.sh



#################
#
# Neighborhood filter
#
#################

$(RES)/putativeins.min.norep.exonic.dist: $(RES)/putativeins.min.norep.exonic
	$(info )
	$(info $(CALL) neighborhood_filter.)
	$(SRC)/neighborhood_filter.sh



#################
#
# Rule to select exonic (and near exonic) abnormal clusters
#
#################

$(RES)/putativeins.min.norep.exonic: $(RES)/putativeins.min.norep $(OUTPUT_DIR)/reference/exons.bed
	$(info )
	$(info $(CALL) Removing clusters extremities not everlapping exons.)
	$(SRC)/select_exonic.sh



#################
#
# Rule to remove clusters overlapping repetitive elements
#
#################

$(RES)/putativeins.min.norep: $(RES)/putativeins.min $(TEMP_PROCESS_DIR)
	$(info )
	$(info $(CALL) Removing clusters overlapping Repetitive Elements annotated by Repeat Masker.)
	$(SRC)/rm_clusters_overlap.sh



#################
#
# Rule to remove clusters with evidence spanning at least 30bp
#
#################

$(RES)/putativeins.min: $(RES)/putativeins
	$(info )
	$(info $(CALL) Removing small ranged clusters and reformatting putativeins.)
	$(SRC)/rm_clusters_evident.sh



#################
#
# Merge abnormal reads into clusters
#
#################

$(RES)/putativeins: $(OUTPUT_DIR)/reference/genes.formated $(RES)
	$(info )
	$(info $(CALL) Clustering abnormals in $(OUTPUT_DIR)/ into $(OUTPUT_DIR)/result/putativeins.)
	$(SRC)/merge.sh



$(FINAL_BAM_LIST): %.format: %.bam %.bai
	$(info )
	@echo "$(timestamp) $(PIPELINE_NAME): Extracting chr format: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE)\n" >> $(LOG_FILE)
	samtools idxstats $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) | cut -c -3 | sort | uniq | grep chr > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).format
	@echo "$(timestamp) $(PIPELINE_NAME): Alignment chr format file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).format\n" >> $(LOG_FILE)


#################
#	
# Extract (chr) format from reference file 
#
#################		
$(OUTPUT_DIR)/$(SAMPLE_ID)/reference.format: $(REFERENCE_GENOME_FASTA)
	$(info )
	@echo "$(timestamp) $(PIPELINE_NAME): Extracting chr format: $(REFERENCE_GENOME_FASTA)\n" >> $(LOG_FILE)
	grep ">" $(REFERENCE_GENOME_FASTA) | sed 's/^>//g' | cut -c -3 | sort | uniq | grep chr > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/reference.format
	@echo "$(timestamp) $(PIPELINE_NAME): Reference chr format file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/reference.format\n" >> $(LOG_FILE)


#################
#	
# Extract (chr) format from annotation file 
#
#################		
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/annotation.format: $(REF_ANNOTATION)
	$(info )
	@echo "$(timestamp) $(PIPELINE_NAME): Extracting chr format: $(REF_ANNOTATION) \n" >> $(LOG_FILE)
	grep -v "^#" $(REF_ANNOTATION) | cut -c -3 | sort | uniq | grep chr > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/annotation.format
	@echo "$(timestamp) $(PIPELINE_NAME): Annotation chr format file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/annotation.format \n" >> $(LOG_FILE)


#################
#
# Rule to format bed file and select protein coding transcripts
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

$(OUTPUT_DIR)/reference/genes.bed: $(REF_ANNOTATION) $(OUTPUT_DIR)/reference
	$(info )
	$(info $(CALL) Created genes bed at $@)
	perl $(SRC)/gtf2bed.pl 'gene' $< > $@


$(OUTPUT_DIR)/reference/exons.bed: $(REF_ANNOTATION) $(OUTPUT_DIR)/reference
	$(info )
	$(info $(CALL) Created exons bed at $@)
	perl $(SRC)/gtf2bed.pl 'exon' $< > $@

#################
#
# Create result, reference and tmp directories.
#
#################

$(OUTPUT_DIR)/result:
	$(info )
	$(info $(CALL) Created result dir $@)
	mkdir -p $@

$(OUTPUT_DIR)/reference:
	$(info )
	$(info $(CALL) Created reference dir $@)
	mkdir -p $@

$(TEMP_PROCESS_DIR):
	$(info )
	$(info $(CALL) Created tmp dir $@)
	mkdir -p $@



.PHONY: mergeCall
