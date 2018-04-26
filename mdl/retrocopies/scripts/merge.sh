#!/bin/bash

###############################################################################
# Arguments usage:

# $1: Source (dir).
# $2: Index file.
# $3: Search criteria (incorrectly used here as the target name).
# $4: Path to module.

# Arguments validation.
if (( ${#@} < 4 )); then
  	echo -e "Error in $0: Not enough arguments!" >&2
	exit 1;
fi



###############################################################################
# Preamble:

# 'chr' pattern extracted from BAM and formated files.
#ptn_bam=$(samtools view -H $1 2> /dev/null | head -n 2 | cut -f 2 | sed -rn '2s/SN:([^0-9]*)[0-9]*/\1/p')
ptn_fmt=$(sed -rn '1s/^(chr|Chr|CHR).*/\1/p' $2)

# Foreach gene
[[ ! -d $1 ]] && echo -e "Error in $0: Inexistent directory $2!" >&2 && exit 1
SAMPLES=$(find -L $1 -type d -name genes)

# Main loop.
while read CHR START END GENE; do
	# Debug code.
	#echo "c-$CHR st-$START ed-$END gn-$GENE; do"
	#echo "$SAMPLES"
	#exit
	
	# Concatenate on 
	echo $GENE >> $3.processed
	
	for i in $SAMPLES; do
		echo $i >> $3.processed
		#SAMPLE_NAME=$(echo $i | awk -F "/" '{print $(NF-1)}')
		[[ -e $i/${GENE}.abnormal ]] && cat $i/${GENE}.abnormal
	done | sort -k3,3V -k4,4n | \
	egrep -v "GL|NC|chrMT|hs|chrM" | \
	perl $4/cluster_pair.pl -w 4000 -s 5 | \
	sort -n -k 11 | \
	awk -v gene="$j" -v start=$START -v end=$END '{ if ( ! ($6 == "=" && (int($7) >= int(start) && int($8) <= int(end)) ) ) {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,gene} else {print $1,$2,$3,$4,$5,"removed"}}'
done < $2

