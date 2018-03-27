#!/bin/bash

###############################################################################
# Arguments usage:

# $1: BAM file.
# $2: Formated file.
# $3: Search criteria.
# $4: OUTPUT_DIR.

# Arguments validation.
if (( ${#@} < 3 )); then
  	echo -e "Error in $0: Not enough arguments!" >&2
	exit 1;
fi



###############################################################################
# Preamble:

# 'chr' pattern extracted from BAM and formated files.
#ptn_bam=$(samtools view -H $1 2> /dev/null | head -n 2 | cut -f 2 | sed -rn '2s/SN:([^0-9]*)[0-9]*/\1/p')
ptn_fmt=$(sed -rn '1s/^(chr|Chr|CHR).*/\1/p' $1)

# Foreach gene
RAND=$RANDOM
[[ ! -d $2 ]] && echo -e "Error in $0: Inexistent directory $2!" >&2 && exit 1
SAMPLES=$(find -L $2 -type d -name genes)

while read CHR START END GENE; do
	# Debug code.
	echo "c-$CHR st-$START ed-$END gn-$GENE; do"
	echo "$SAMPLES"
	exit
	
	# Concatenate on 
	#echo $GENE >> $@.processed;
	
	for i in $SAMPLES; do
		echo $i >> $@.processed;
		SAMPLE_NAME=$(echo $i | awk -F "/" '{print $(NF-1)}')
		cat $i/${GENE}.abnormal
	done | sort -k3,3V -k4,4n | \
	egrep -v "GL|NC|chrMT|hs|chrM" | \
	perl $3/cluster_pair.pl -w 4000 -s 5 | \
	sort -n -k 11 | \
	awk -v gene="$j" -v start=$START -v end=$END '{ if ( ! ($6 == "=" && (int($7) >= int(start) && int($8) <= int(end)) ) ) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,gene} else {print $1,$2,$3,$4,$5,"removed"}}'
done < $1

