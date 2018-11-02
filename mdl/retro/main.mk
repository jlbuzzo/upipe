############################## HEADER #########################################
# Makefile to rulle'em all.
#
# Data:
# 		Main target: validations
# 		File: validations.mk
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

# Default goal
main:



############################## VALIDATIONS ####################################

# procedural pre-processing part.

# Validate INPUT against emptyness.
INPUT_PROCESSED := $(strip $(INPUT))
$(if $(INPUT_PROCESSED),, $(error INPUT variable is empty))

# Search and validate INPUT from a given string list, directory, or file, by
# the SUFFIXES list criteria.
#INPUT_PROCESSED := $(shell $(SCRIPTS)/ufinder.sh $(INPUT_PROCESSED) $(SUFFIXES))

# A limited alternative validation process (without 'ufinder' script).
INPUT_PROCESSED := $(if $(shell [ -f $(INPUT_PROCESSED) ] && echo 1), $(wildcard $(filter %$(SUFFIXES), $(shell cat $(INPUT_PROCESSED)))), $(wildcard $(INPUT_PROCESSED)/*$(SUFFIXES)))
INPUT_PROCESSED := $(abspath $(strip $(INPUT_PROCESSED)))

# Verifying 'SUFFIXES' terminated files' remaining after the filtering of the
# INPUT_PROCESSED list.
$(if $(INPUT_PROCESSED),, $(error No valid $(SUFFIXES) file specified))

# Verifying OUTPUT_DIR previous existence.
ifneq ($(wildcard $(OUTPUT_DIR)),)
$(info )
$(info $(CALL) Directory '$(OUTPUT_DIR)' will be overwritten!)
endif

# ===> DEBUG CODE <===

# Enable emergy stop, if variable STP1="yes".
ifdef STP1
$(error Emergency stop 1)
endif


############################## PROCESSSAMPLE ##################################

# procedural pre-processing part.

# ATTENTION: Variables that are string lists has an '_l' suffix appendedd to
# the end of their names. Their behavior differently in pattern rules.

# Test the succes of the validations process.
$(if $(INPUT_PROCESSED),, $(error "No valid INPUT"))

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

# Look for a *.sorted.bam.bai file, according to the link to 'OUTPUT_l'.
REQ_BAI = $(shell readlink -f $(filter %$(*F)$(SUFFIXES), $(OUTPUT_l)))
aux2 = $(strip $(filter %.sorted.bam, $(REQ_BAI)))
STBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*sorted*.bai')
NSTBAI = $(shell find $(dir $(REQ_BAI)) -type f -name '$(*F)*.bai' -not -regex '.*sorted.*')
BAI = $(if $(aux2),$(STBAI),$(NSTBAI))


############################## OTHERS #########################################

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
	$(info )
	$(info $(CALL) Main done!)



validations:
	$(info )
	$(info $(CALL) Target 'validations' complete!)



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

###############################################################################
#
# mergeCall
#
###############################################################################



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
	$(info $(CALL) Make genes.bed.)
	perl $(SCRIPTS)/gtf2bed.pl 'gene' $< > $@
	@echo "$(CALL): Created genes bed at: $@.\n"



$(REF)/exons.bed: $(REF_ANNOTATION) | $(REF)
	$(info )
	$(info $(CALL) Make exons.bed.)
	perl $(SCRIPTS)/gtf2bed.pl 'exon' $< > $@
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

