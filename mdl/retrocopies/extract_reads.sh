#!/bin/bash/

# Main variable attributions.
goal=$1
req=$2
out_dir=$3


# Removing remnant directories.
rm -Rf $(OUTPUT_DIR)/result/dump/
rm -Rf $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation
mkdir -p $(OUTPUT_DIR)/result/dump/

# wsgweg
for putativeins in $$(cat $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar | sed 's/[ ]/#/g' ); do
	GENE=$(echo $putativeins | awk -F "[#]" '{print $$12}')
	PAR_CHR=$(echo $putativeins | awk -F "[#]" '{print $$1}'| sed 's/chr//g')
	PAR_START=$(echo $putativeins | awk -F "[#]" '{print $$2}')
	PAR_END=$(echo $putativeins | awk -F "[#]" '{print $$3}')
	IP_CHR=$(echo $putativeins | awk -F "[#]" '{print $$6}'|sed 's/chr//g')
	PAR_CHR2=$(echo $putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$1}}'|sed 's/chr//g')
	IP_CHR2=$(echo $putativeins | awk -F "[#]" '{ if ($$6 == $$1) {print "="} else {print $$6}}')|sed 's/chr//g'
	IP_START=$(echo $putativeins | awk -F "[#]" '{print $$7-500}')
	IP_END=$(echo $putativeins | awk -F "[#]" '{print $$8+500}')
	DESC=$(echo $$putativeins | sed 's/#/ /g')
  		GENE2=$(echo $$putativeins | awk -F "**|@@" '{print $$4}')
 	grep -P "\t$GENE2\_" $(OUTPUT_DIR)/reference/exons.bed > $(OUTPUT_DIR)/result/dump/temp_exons.txt
	SAMPLES=$(find $(OUTPUT_DIR)/ -type d -name "*.bam" -o -name "*.cram")


	# fdsfdg
	for i in $$SAMPLES; do
		SAMPLE=$$(echo $$i | awk -F "/" '{print $$(NF)}')
		samtools view $$i/$$SAMPLE.sorted.bam $$PAR_CHR\:$$PAR_START\-$$PAR_END | awk -v chr=$$IP_CHR2 -v start=$$IP_START -v end=$$IP_END '{if ($$7 == chr && int($$8) >= int(start) && int($$8) <= int(end) ) {print}}' | sed 's/^/'$$SAMPLE' /'
		samtools view $$i/$$SAMPLE.sorted.bam $$IP_CHR\:$$IP_START\-$$IP_END | awk -v chr=$$PAR_CHR2 -v start=$$PAR_START -v end=$$PAR_END '{if ($$7 == chr && int($$8) >= int(start) && int($$8) <= int(end) ) {print}}' | sed 's/^/'$$SAMPLE' /'
	done >> $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal


	# sdfgsfgsd
	cat $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal | perl library/src/or_ip.pl -d "$$DESC" -f $(OUTPUT_DIR)/result/dump/temp_exons.txt > $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal.bed 2>> $(OUTPUT_DIR)/result/putativeins.min.norep.exonic.dist.notsimilar.orientation
	sed -i 's/chrchr/chr/g' $(OUTPUT_DIR)/result/dump/$$GENE\_$$PAR_CHR\_$$PAR_START\_$$IP_CHR\_$$IP_START.reads.abnormal.bed
done


exit
