###############################################################################
#
# mergeCall
#
###############################################################################

# Preamble.



#################
#
# Main target: mergeCall
#
#################

mergeCall: $(RES)/putativeins.min.norep.exonic.dist.notsimilar.orientation
	$(info )
	$(info $(CALL) Target 'mergeCall' complete!)



#################
#
# putativeins.min.norep.exonic.dist.notsimilar.orientation
#
#################

$(RES)/putativeins.min.norep.exonic.dist.notsimilar.orientation: $(RES)/putativeins.min.norep.exonic.dist.notsimilar
	$(info )
	$(info $(CALL)  Make $@.)
	@echo "$(CALL): Recalculating Insertion Point and Support from original BAM files; Checking supporting reads orientation\n"
	rm -Rf $(RES)/dump/; \
	rm -Rf $@; \
	mkdir -p $(RES)/dump/; \
	for putativeins in $$(cat $< | sed 's/[ ]/#/g' ); do \
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
		grep -P "\t$$GENE2\_" $(REF)/exons.bed > $(RES)/dump/temp_exons.txt; \
		SAMPLES=$$(find $(OUTPUT_DIR)/ -type d -name "*.bam" -o -name "*.cram"); \
		for i in $$SAMPLES; do \
			SAMPLE=$$(echo $$i | awk -F "/" '{print $$(NF)}'); \
			samtools view $$i/$$SAMPLE $$PAR_CHR\:$$PAR_START\-$$PAR_END | awk -v chr=$$IP_CHR2 -v start=$$IP_START -v end=$$IP_END '{if ($$7 == chr && int($$8) >= int(start) && int($$8) <= int(end) ) {print}}' | sed 's/^/'$$SAMPLE' /'; \
			samtools view $$i/$$SAMPLE $$IP_CHR\:$$IP_START\-$$IP_END | awk -v chr=$$PAR_CHR2 -v start=$$PAR_START -v end=$$PAR_END '{if ($$7 == chr && int($$8) >= int(start) && int($$8) <= int(end) ) {print}}' | sed 's/^/'$$SAMPLE' /'; \
		done >> $(RES)/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal; \
		cat $(RES)/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal | perl $(SCRIPTS)/or_ip.pl -d "$$DESC" -f $(RES)/dump/temp_exons.txt > $(RES)/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal.bed 2>> $@; \
		sed -i 's/chrchr/chr/g' $(RES)/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal.bed; \
	done 
	@echo "$(CALL): Finished recalculating Insertion Point and Support from original BAM files; Checked supporting reads orientation.\n"



#################
#
# Rule to select non similar IP/PARENTAL sequences 
#
#################

# apaguei temporariamente o dist do $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.DIST e na linha do comando perl 254!!!
$(RES)/putativeins.min.norep.exonic.dist.notsimilar: $(RES)/putativeins.min.norep.exonic | $(TEMP_PROCESS_DIR)
	$(info )
	$(info $(CALL) Make $@.)
	@echo "$(CALL): Removing clusters which Insertion Point is simmilar to Parental Sequence.\n"
	RAND=$$$$ && \
	perl $(SCRIPTS)/similarity_filter.pl -p $(TEMP_PROCESS_DIR) -f $< -g $(REFERENCE_GENOME_FASTA) > $@_debug
	awk '{if ($$NF == "IN") {print}}' $@_debug > $@
	@echo "$(CALL): Removed clusters which Insertion Point is simmilar to Parental Sequence.\n"



#################
#
# Neighborhood filter
##
#################

#$(RES)/putativeins.min.norep.exonic.dist: $(RES)/putativeins.min.norep.exonic
	#$(info )
	#$(info $(CALL) Make $@.)
	#cat $< | awk '{if ($$1 != $$6) {print $$_,"IN_DIST"} else { if ($$2 - $$8 >= $(MIN_DIST) && $$7 - $$3 >= $(MIN_DIST) ) {print $$_,"IN_DIST" } else {print $$_,"OUT_DIST"}} }' > $@_debug
	#grep -v "OUT_DIST" $@_debug > $@



#################
#
# Rule to select exonic (and near exonic) abnormal clusters
#
#################
$(RES)/putativeins.min.norep.exonic: $(RES)/putativeins.min.norep $(REF)/exons.bed | $(TEMP_PROCESS_DIR)
	$(info )
	$(info $(CALL) Make $@.)
	@echo "$(CALL): Removing clusters extremities not everlapping exons.\n"
	RAND=$$$$ && \
	for putativeins in $$(cat $< | sed 's/[ ]/#/g' ); do \
		GENE=$$(echo $$putativeins | awk -F "**|@@" '{print $$4}'); \
		grep -P "\t$$GENE\_" $(REF)/exons.bed > $(TEMP_PROCESS_DIR)/temp_exons.txt.$$RAND; \
		echo $$putativeins | awk -F "#" '{print $$1"\t"$$2-20"\t"$$2+20"\tUP\n"$$1"\t"$$3-20"\t"$$3+20"\tDOWN"}' > $(TEMP_PROCESS_DIR)/a.$$RAND; \
		BOTH=`intersectBed -a $(TEMP_PROCESS_DIR)/a.$$RAND -b $(TEMP_PROCESS_DIR)/temp_exons.txt.$$RAND | awk '{print $$NF}' | sort | uniq | wc -l`; \
		if [ $$BOTH -eq 2 ]; then \
			echo "$$putativeins" | sed 's/#/ /g'; \
		fi; \
	done > $@
	@echo "$(CALL): Removed clusters extremities not everlapping exons."
	#rm $(TEMP_PROCESS_DIR)/a.$$RAND; \
	#rm $(TEMP_PROCESS_DIR)/temp_exons.txt.$$RAND;



#################
#
# Rule to remove clusters overlapping repetitive elements
#
#################
$(RES)/putativeins.min.norep: $(RES)/putativeins.min | $(TEMP_PROCESS_DIR)
	$(info )
	$(info $(CALL) Make $@.)
	@echo "$(CALL): Removing clusters overlapping Repetitive Elements annotated by Repeat Masker.\n"
	perl $(SCRIPTS)/remove_rep.pl -p $(TEMP_PROCESS_DIR) -f $(REP_ANNOTATION) -f2 $<
	perl $(SCRIPTS)/remove_rep.pl -p $(TEMP_PROCESS_DIR) -f $(REP_ANNOTATION_manual) -f2 $@
	mv $@.norep $@
	@echo "$(CALL): Removed clusters overlapping Repetitive Elements annotated by Repeat Masker.\n"
	##This should be very slow.. Is there anyway of parallezing it?



#################
#
# Rule to remove clusters with evidence spanning at least 30bp
#
#################

$(RES)/putativeins.min: $(RES)/putativeins
	$(info )
	$(info $(CALL) Make $@.)
	@echo "$(CALL): Removing small ranged clusters and reformatting putativeins.\n"
	exit
	cat $< | awk '{if ($$6 == "chr=") {print $$1,$$2,$$3,$$4,$$5,$$1,$$7,$$8,$$9,$$10,$$11,$$12} else {print $$_}}' > $(RES)/temp 
	cat $(RES)/temp | sed 's/[()]//g' | awk '{if ( $$4 >= 30 && $$8 >= 30 ) {print}}' > $@
	rm -r $(RES)/temp
	@echo "$(CALL): Clustering abnormals in $(OUTPUT_DIR)/ into $<.\n"



#################
#
# Merge abnormal reads into clusters
#
#################

$(RES)/putativeins: processSample $(REF)/genes.bed | $(RES)
	$(info )
	$(info $(CALL) Make $@.)
	@echo "$(CALL): Clustering abnormals in $(OUTPUT_DIR)/ into $@.\n"
	$(SCRIPTS)/merge.sh $(OUTPUT_DIR) $(word 2,$^) $@ $(SCRIPTS) > $@
	@echo "$(CALL): Finished clustering abnormals in $(OUTPUT_DIR)/ into $@.\n"
