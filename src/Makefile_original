##
## Main make target
##


$(info ###############################################################################)
$(info                               UFINDER)
$(info ###############################################################################)
$(info processSample)

all : processSample

.PHONY: all processSample
.DEFAULT: all
.INTERMEDIATE: 

include opt.mk

FILE_NAME:=$(notdir $(INPUT_FILE))
ifneq ($(SAMPLE_NAME),NULL)
       SAMPLE_ID := $(FILE_NAME)
endif

LOG_FILE := $(OUTPUT_DIR)/$(SAMPLE_ID).log
timestamp := `/bin/date "+%Y-%m-%d(%H:%M:%S)"`

##
## Create results directory
##

$(OUTPUT_DIR)/$(SAMPLE_ID):
	$(info )
	$(info Make output dir.)
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)
	@echo "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(LOG_FILE)

$(OUTPUT_DIR)/reference:
	$(info )
	$(info Make reference dir.)
	mkdir -p $(OUTPUT_DIR)/reference
	@echo "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(LOG_FILE)

$(OUTPUT_DIR)/result:
	$(info )
	$(info Make result dir.)
	mkdir -p $(OUTPUT_DIR)/result
	@echo "$(timestamp) $(PIPELINE_NAME): Created results dir: $(OUTPUT_DIR)/$(SAMPLE_ID)\n" >> $(LOG_FILE)

$(TEMP_PROCESS_DIR):
	$(info )
	$(info Make tmp dir.)
	mkdir -p $(TEMP_PROCESS_DIR)

#################
#
# Rule to sort .bam file
#
#################
ifeq ($(SORTED),TRUE)
FINAL_BAM_FILE := $(FILE_NAME)
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE): $(INPUT_FILE) | $(OUTPUT_DIR)/$(SAMPLE_ID)
	$(info )
	$(info Make *.bam.)
	ln -s $(INPUT_FILE) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) 
	@echo "$(timestamp) $(PIPELINE_NAME): Created link at $(INPUT_FILE) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE): \n" >> $(LOG_FILE)
else
FINAL_BAM_FILE := $(FILE_NAME).sorted.bam
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE): $(INPUT_FILE) | $(OUTPUT_DIR)/$(SAMPLE_ID)
	$(info )
	$(info Make *.bam.)
	@echo "$(timestamp) $(PIPELINE_NAME): Sorting file: $(INPUT_FILE)\n" >> $(LOG_FILE)
	samtools sort -O BAM -m 10G $(INPUT_FILE) -o $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FILE_NAME).sorted
	@echo "$(timestamp) $(PIPELINE_NAME): Sorted file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE)\n" >> $(LOG_FILE)
endif

#################
#	
# Rule to index .bam file
#
#################		

$(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).bai: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE)
	$(info )
	$(info Make *.bai.)
	@echo "$(timestamp) $(PIPELINE_NAME): Indexing file: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE)\n" >> $(LOG_FILE)
	samtools index -b $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).bai
	@echo "$(timestamp) $(PIPELINE_NAME): Indexed file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).bai\n" >> $(LOG_FILE)


#################
#	
# Extract (chr) format from alignment file 
#
#################		

$(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).format: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).bai
	$(info )
	$(info Make *.format.)
	@echo "$(timestamp) $(PIPELINE_NAME): Extracting chr format: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE)\n" >> $(LOG_FILE)
	samtools idxstats $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) | cut -c -3 | sort | uniq | grep chr > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).format
	@echo "$(timestamp) $(PIPELINE_NAME): Alignment chr format file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).format\n" >> $(LOG_FILE)


#################
#	
# Extract (chr) format from reference file 
#
#################		
$(OUTPUT_DIR)/$(SAMPLE_ID)/reference.format: $(REFERENCE_GENOME_FASTA)
	@echo "$(timestamp) $(PIPELINE_NAME): Extracting chr format: $(REFERENCE_GENOME_FASTA)\n" >> $(LOG_FILE)
	grep ">" $(REFERENCE_GENOME_FASTA) | sed 's/^>//g' | cut -c -3 | sort | uniq | grep chr > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/reference.format
	@echo "$(timestamp) $(PIPELINE_NAME): Reference chr format file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/reference.format\n" >> $(LOG_FILE)


#################
#	
# Extract (chr) format from annotation file 
#
#################		
$(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/annotation.format: $(REF_ANNOTATION)
	@echo "$(timestamp) $(PIPELINE_NAME): Extracting chr format: $(REF_ANNOTATION) \n" >> $(LOG_FILE)
	grep -v "^#" $(REF_ANNOTATION) | cut -c -3 | sort | uniq | grep chr > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/annotation.format
	@echo "$(timestamp) $(PIPELINE_NAME): Annotation chr format file at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(OUTPUT_DIR)/reference/annotation.format \n" >> $(LOG_FILE)


#################
#
# Rule to find abnormal alignments
#
#################

$(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).abnormal: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).bai $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) $(OUTPUT_DIR)/reference/genes.formated
	$(info )
	$(info Make *.abnormal.)
	@echo "$(timestamp) $(PIPELINE_NAME): Selecting abnormal pairs from file: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE)\n" >> $(LOG_FILE)
	#Clean and create temporary path to store temporary files
	rm -Rf $(OUTPUT_DIR)/$(SAMPLE_ID)/genes
	mkdir -p $(OUTPUT_DIR)/$(SAMPLE_ID)/genes
	#For each gene, create a temporary file with abnormal alignments
	for j in $$(cat $(OUTPUT_DIR)/reference/genes.formated ); do \
      GENE=$$(echo $$j | awk -F "[*][*]" '{print $$4}'); \
			CHR=$$(echo $$j | awk -F "[*][*]" '{print $$1}'); \
			START=$$(echo $$j | awk -F "[*][*]" '{print $$2}'); \
      END=$$(echo $$j | awk -F "[*][*]" '{print $$3}'); \
		samtools view $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE) $$CHR:$$START-$$END | \
		awk '{if ( ( $$7 == "=" && ( $$9 > 100000 || $$9 < -100000) ) || ( $$7 != "=" && $$7 != "*" ) ) {print "chr"$$3"\t"$$4"\tchr"$$7"\t"$$8}}' | sed 's/chrchr/chr/g' > $(OUTPUT_DIR)/$(SAMPLE_ID)/genes/$$GENE.abnormal 2>/dev/null; \
      if [[ ! -s $(OUTPUT_DIR)/$(SAMPLE_ID)/genes/$$GENE.abnormal ]]; then rm $(OUTPUT_DIR)/$(SAMPLE_ID)/genes/$$GENE.abnormal; fi; \
    done
	find $(OUTPUT_DIR)/$(SAMPLE_ID)/genes/ -type f -name '*.abnormal' -exec wc -l '{}' ';' >> $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).abnormal
	#wc -l $(OUTPUT_DIR)/$(SAMPLE_ID)/genes/*.abnormal > $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).abnormal
	@echo "$(timestamp) $(PIPELINE_NAME): Abnormal file create at: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).abnormal\n" >> $(LOG_FILE)

#################
#
# Rule to process .gtf
#
#################

$(OUTPUT_DIR)/reference/genes.bed: $(REF_ANNOTATION) | $(OUTPUT_DIR)/reference
	$(info )
	$(info Make genes.bed.)
	perl library/src/gtf2bed.pl 'gene' $< > $@
	@echo "$(timestamp) $(PIPELINE_NAME): Created genes bed at: $(OUTPUT_DIR)/reference/genes.bed\n" >> $(LOG_FILE)

#cat $(REF_ANNOTATION) | awk '{if ( $$3 == "gene" ) {print $$1"\t"$$4"\t"$$5"\t"$$16"@@"$$12"@@"$$18"@@"$$10}}' | sed 's/[";]//g' > $(OUTPUT_DIR)/reference/genes.bed

$(OUTPUT_DIR)/reference/exons.bed: $(REF_ANNOTATION) | $(OUTPUT_DIR)/reference
	$(info )
	$(info Make exons.bed.)
	perl library/src/gtf2bed.pl 'exon' $< > $@
	@echo "$(timestamp) $(PIPELINE_NAME): Created exons bed at: $(OUTPUT_DIR)/reference/exons.bed\n" >> $(LOG_FILE)

#cat $(REF_ANNOTATION) | awk '{if ( $$3 == "exon" ) {print $$1"\t"$$4"\t"$$5"\t"$$18"_"$$14"_"$$26"_"$$12}}' | sed 's/[";]//g' > $(OUTPUT_DIR)/reference/exons.bed

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
# Merge abnormal reads into clusters
#
#################

$(OUTPUT_DIR)/result/putativeins: $(OUTPUT_DIR)/reference/genes.formated | $(OUTPUT_DIR)/result
	$(info )
	$(info Make putativeins.)
	#Foreach gene
	@echo "$(timestamp) $(PIPELINE_NAME): Clustering abnormals in $(OUTPUT_DIR)/ into $(OUTPUT_DIR)/result/putativeins\n" #>> $(LOG_MERGE)
	RAND=$$(echo $$RANDOM);
	SAMPLES=$$(find -L $(OUTPUT_DIR)/ -type d -name genes); \
	for j in $$(cat $(OUTPUT_DIR)/reference/genes.formated ); do \
		GENE=$$(echo $$j | awk -F "[*][*]" '{print $$4}'); \
		echo $$GENE >> $(OUTPUT_DIR)/result/putativeins.processed; \
		CHR=$$(echo $$j | awk -F "[*][*]" '{print $$1}') ; \
		START=$$(echo $$j | awk -F "[*][*]" '{print $$2}'); \
		END=$$(echo $$j | awk -F "[*][*]" '{print $$3}'); \
		for i in $$SAMPLES; do \
			echo $$i >> $(OUTPUT_DIR)/result/putativeins.processed; \
			SAMPLE_NAME=$$(echo $$i | awk -F "/" '{print $$(NF-1)}');\
			cat $$i/"$$GENE".abnormal; \
		done | \
		sort -k3,3V -k4,4n | \
		egrep -v "GL|NC|chrMT|hs|chrM" | \
		perl library/src/cluster_pair.pl -w 4000 -s 5| \
		sort -n -k 11 | \
		awk -v gene="$$j" -v start=$$START -v end=$$END '{ if ( ! ($$6 == "=" && (int($$7) >= int(start) && int($$8) <= int(end)) ) ) {print $$1,$$2,$$3,$$4,$$5,$$6,$$7,$$8,$$9,$$10,$$11,gene} else {print $$1,$$2,$$3,$$4,$$5,"removed"}}'; \
	done > $(OUTPUT_DIR)/result/putativeins;
	@echo "$(timestamp) $(PIPELINE_NAME): Finished clustering abnormals in $(OUTPUT_DIR)/ into $(OUTPUT_DIR)/result/putativeins\n" #>> $(LOG_MERGE)


####################################################################################################
#################
#
# Rule to remove clusters with evidence spanning at least 30bp
#
#################
$(OUTPUT_DIR)/result/putativeins.min: $(OUTPUT_DIR)/result/putativeins
	$(info )
	$(info Make putativeins.min.)
	@echo "$(timestamp) $(PIPELINE_NAME): Removing small ranged clusters and reformatting putativeins\n" #>> $(LOG_MERGE)
	cat $(OUTPUT_DIR)/result/putativeins  | awk '{if ($$6 == "chr=") {print $$1,$$2,$$3,$$4,$$5,$$1,$$7,$$8,$$9,$$10,$$11,$$12} else {print $$_}}' > $(OUTPUT_DIR)/result/temp 
	cat $(OUTPUT_DIR)/result/temp | sed 's/[()]//g' | awk '{if ( $$4 >= 30 && $$8 >= 30 ) {print}}' > $(OUTPUT_DIR)/result/putativeins.min
	rm -r $(OUTPUT_DIR)/result/temp
	@echo "$(timestamp) $(PIPELINE_NAME): Clustering abnormals in $(OUTPUT_DIR)/ into $(OUTPUT_DIR)/result/putativeins\n" #>> $(LOG_MERGE)


#################
#
# Rule to remove clusters overlapping repetitive elements
#
#################
$(OUTPUT_DIR)/result/putativeins.min.norep: $(OUTPUT_DIR)/result/putativeins.min | $(TEMP_PROCESS_DIR)
	$(info )
	$(info Make putativeins.min.norep.)
	@echo "$(timestamp) $(PIPELINE_NAME): Removing clusters overlapping Repetitive Elements annotated by Repeat Masker\n" #>> $(LOG_MERGE)
	perl library/src/remove_rep.pl -p $(TEMP_PROCESS_DIR) -f $(REP_ANNOTATION) -f2 $(OUTPUT_DIR)/result/putativeins.min
	perl library/src/remove_rep.pl -p $(TEMP_PROCESS_DIR) -f $(REP_ANNOTATION_manual) -f2 $(OUTPUT_DIR)/result/putativeins.min.norep
	mv $(OUTPUT_DIR)/result/putativeins.min.norep.norep $(OUTPUT_DIR)/result/putativeins.min.norep
	@echo "$(timestamp) $(PIPELINE_NAME): Removed clusters overlapping Repetitive Elements annotated by Repeat Masker\n" #>> $(LOG_MERGE)
	##This should be very slow.. Is there anyway of parallezing it?

#################
#
# Rule to select exonic (and near exonic) abnormal clusters
#
#################
$(OUTPUT_DIR)/result/putativeins.min.norep.exonic: $(OUTPUT_DIR)/result/putativeins.min.norep $(OUTPUT_DIR)/reference/exons.bed | $(TEMP_PROCESS_DIR)
	$(info )
	$(info Make putativeins.min.norep.exonic.)
	@echo "$(timestamp) $(PIPELINE_NAME): Removing clusters extremities not everlapping exons\n" #>> $(LOG_MERGE)
	RAND=$$$$ && \
	for putativeins in $$(cat $(OUTPUT_DIR)/result/putativeins.min.norep | sed 's/[ ]/#/g' ); do \
		GENE=$$(echo $$putativeins | awk -F "**|@@" '{print $$4}'); \
		grep -P "\t$$GENE\_" $(OUTPUT_DIR)/reference/exons.bed > $(TEMP_PROCESS_DIR)/temp_exons.txt.$$RAND; \
		echo $$putativeins | awk -F "#" '{print $$1"\t"$$2-20"\t"$$2+20"\tUP\n"$$1"\t"$$3-20"\t"$$3+20"\tDOWN"}' > $(TEMP_PROCESS_DIR)/a.$$RAND; \
		BOTH=`intersectBed -a $(TEMP_PROCESS_DIR)/a.$$RAND -b $(TEMP_PROCESS_DIR)/temp_exons.txt.$$RAND | awk '{print $$NF}' | sort | uniq | wc -l`; \
		if [ $$BOTH -eq 2 ]; then \
			echo "$$putativeins" | sed 's/#/ /g'; \
		fi; \
	done > $(OUTPUT_DIR)/result/putativeins.min.norep.exonic
	@echo "$(timestamp) $(PIPELINE_NAME): Removed clusters extremities not everlapping exons." # >> $(LOG_MERGE)
	#rm $(TEMP_PROCESS_DIR)/a.$$RAND; \
	#rm $(TEMP_PROCESS_DIR)/temp_exons.txt.$$RAND;

#################
#
# Neighborhood filter
##
#################

#$(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist: $(OUTPUT_DIR)/result/putativeins.min.norep.exonic
#	cat $(OUTPUT_DIR)/result/putativeins.min.norep.exonic | awk '{if ($$1 != $$6) {print $$_,"IN_DIST"} else { if ($$2 - $$8 >= $(MIN_DIST) && $$7 - $$3 >= $(MIN_DIST) ) {print $$_,"IN_DIST" } else {print $$_,"OUT_DIST"}} }' > $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist_debug
#	grep -v "OUT_DIST" $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist_debug > $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist


#################
#
# Rule to select non similar IP/PARENTAL sequences 
#
#################
## apaguei temporariamente o dist do $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.DIST e na linha do comando perl 254!!!
##

$(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar: $(OUTPUT_DIR)/result/putativeins.min.norep.exonic | $(TEMP_PROCESS_DIR)
	$(info )
	$(info Make putativeins.min.norep.exonic.dist.notsimilar.)
	@echo "$(timestamp) $(PIPELINE_NAME): Removing clusters which Insertion Point is simmilar to Parental Sequence\n" #>> $(LOG_MERGE)
	RAND=$$$$ && \
	perl library/src/similarity_filter.pl -p $(TEMP_PROCESS_DIR) -f $(OUTPUT_DIR)/result/putativeins.min.norep.exonic -g $(REFERENCE_GENOME_FASTA) > $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar_debug
	awk '{if ($$NF == "IN") {print}}' $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar_debug > $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar
	@echo "$(timestamp) $(PIPELINE_NAME): Removed clusters which Insertion Point is simmilar to Parental Sequence\n" #>> $(LOG_MERGE)

#################
#
# Extract reads from bam files
## 
#################

		#PAR_CHR=$$(echo $$putativeins | awk -F "[#]" '{print $$1}'| sed 's/chr//g'); \
		#IP_CHR=$$(echo $$putativeins | awk -F "[#]" '{print $$6}'|sed 's/chr//g'); \
		#PAR_CHR2=$$(echo $$putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$1}}'|sed 's/chr//g'); \
		#IP_CHR2=$$(echo $$putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$6}}')|sed 's/chr//g'; \

$(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation: $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar
	$(info )
	$(info Make putativeins.min.norep.exonic.dist.notsimilar.orientation.)
	@echo "$(timestamp) $(PIPELINE_NAME): Recalculating Insertion Point and Support from original BAM files; Checking supporting reads orientation\n" #>> $(LOG_MERGE)
	rm -Rf $(OUTPUT_DIR)/result/dump/; \
	rm -Rf $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation; \
	mkdir -p $(OUTPUT_DIR)/result/dump/; \
	for putativeins in $$(cat $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar | sed 's/[ ]/#/g' ); do \
		GENE=$$(echo $$putativeins | awk -F "[#]" '{print $$12}'); \
		PAR_CHR=$$(echo $$putativeins | awk -F "[#]" '{print $$1}'); \
		PAR_START=$$(echo $$putativeins | awk -F "[#]" '{print $$2}'); \
		PAR_END=$$(echo $$putativeins | awk -F "[#]" '{print $$3}'); \
		IP_CHR=$$(echo $$putativeins | awk -F "[#]" '{print $$6}'); \
		PAR_CHR2=$$(echo $$putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$1}}'); \
		IP_CHR2=$$(echo $$putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$6}}'); \
		IP_START=$$(echo $$putativeins | awk -F "[#]" '{print $$7-500}'); \
		IP_END=$$(echo $$putativeins | awk -F "[#]" '{print $$8+500}'); \
		DESC=$$(echo $$putativeins | sed 's/#/ /g'); \
		GENE2=$$(echo $$putativeins | awk -F "**|@@" '{print $$4}'); \
		grep -P "\t$$GENE2\_" $(OUTPUT_DIR)/reference/exons.bed > $(OUTPUT_DIR)/result/dump/temp_exons.txt; \
		SAMPLES=$$(find $(OUTPUT_DIR)/ -type d -name "*.bam" -o -name "*.cram"); \
		for i in $$SAMPLES; do \
			SAMPLE=$$(echo $$i | awk -F "/" '{print $$(NF)}'); \
			samtools view $$i/$$SAMPLE $$PAR_CHR\:$$PAR_START\-$$PAR_END | awk -v chr=$$IP_CHR2 -v start=$$IP_START -v end=$$IP_END '{if ($$7 == chr && int($$8) >= int(start) && int($$8) <= int(end) ) {print}}' | sed 's/^/'$$SAMPLE' /'; \
			samtools view $$i/$$SAMPLE $$IP_CHR\:$$IP_START\-$$IP_END | awk -v chr=$$PAR_CHR2 -v start=$$PAR_START -v end=$$PAR_END '{if ($$7 == chr && int($$8) >= int(start) && int($$8) <= int(end) ) {print}}' | sed 's/^/'$$SAMPLE' /'; \
		done >> $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal; \
		cat $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal | perl library/src/or_ip.pl -d "$$DESC" -f $(OUTPUT_DIR)/result/dump/temp_exons.txt > $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal.bed 2>> $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation; \
		sed -i 's/chrchr/chr/g' $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal.bed; \
	done 
	@echo "$(timestamp) $(PIPELINE_NAME): Finished recalculating Insertion Point and Support from original BAM files; Checked supporting reads orientation\n" #>> $(LOG_MERGE)


#################
#
# no evidence - test
##
#################

#$(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation.noevidence: $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar
#	@echo "$(timestamp) $(PIPELINE_NAME): Recalculating Insertion Point and Support from original BAM files; Checking supporting reads orientation\n" >> $(LOG_MERGE)
#	rm -Rf $(OUTPUT_DIR)/result/dump/;
#	rm -Rf $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation; 
#	mkdir -p $(OUTPUT_DIR)/result/dump/;
#	for putativeins in $$(cat $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar | sed 's/[ ]/#/g' ); do \
#		GENE=$$(echo $$putativeins | awk -F "[#]" '{print $$12}'); \
#		PAR_CHR=$$(echo $$putativeins | awk -F "[#]" '{print $$1}'|sed 's/chr//'); \
#		PAR_START=$$(echo $$putativeins | awk -F "[#]" '{print $$2}'); \
#		PAR_END=$$(echo $$putativeins | awk -F "[#]" '{print $$3}'); \
#		IP_CHR=$$(echo $$putativeins | awk -F "[#]" '{print $$6}'|sed 's/chr//'); \
#		PAR_CHR2=$$(echo $$putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$1}}'|sed 's/chr//'); \
#		IP_CHR2=$$(echo $$putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$6}}'|sed 's/chr//'); \
#		IP_START=$$(echo $$putativeins | awk -F "[#]" '{print $$7-500}'); \
#		IP_END=$$(echo $$putativeins | awk -F "[#]" '{print $$8+500}'); \
#		DESC=$$(echo $$putativeins | sed 's/#/ /g'); \
#  	GENE2=$$(echo $$putativeins | awk -F "**|@@" '{print $$4}'); \
#	 	grep -P "\t$$GENE2\_" $(OUTPUT_DIR)/reference/exons.bed > $(OUTPUT_DIR)/result/dump/temp_exons.txt; \
#		SAMPLES=$$(find $(OUTPUT_DIR)/ -type d -name "*.bam"); \
#		for i in $$SAMPLES; do \

processSample: $(OUTPUT_DIR)/$(SAMPLE_ID)/$(FINAL_BAM_FILE).abnormal
	$(info )
	$(info Finished the search for abnormals.)

mergeCall: $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation
	$(info )
	$(info Finished merge.)
