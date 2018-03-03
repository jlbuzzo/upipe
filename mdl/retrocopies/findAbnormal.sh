#!/bin/bash

# Arguments:
# $1: OUTPUT_DIR.
# $2: FINAL_BAM_FILE.
# $3: *.formated file.


# Clean and create temporary path to store temporary files.
rm -Rf $1/genes
mkdir -p $1/genes

# For each gene, create a temporary file with abnormal alignments.
while read line; do
	GENE=$(echo $line | awk -F "[*][*]" '{print $4}')
	CHR=$(echo $line | awk -F "[*][*]" '{print $1}' | sed 's/chr//i')
	START=$(echo $line | awk -F "[*][*]" '{print $2}')
	END=$(echo $line | awk -F "[*][*]" '{print $3}')

	samtools view $2 $CHR:$START-$END | awk '
		{
			if ( ( $7 == "=" && ( $9 > 100000 || $9 < -100000) ) || ( $7 != "=" && $7 != "*" ) )
				print "chr"$3"\t"$4"\tchr"$7"\t"$8
		}' | sed 's/chrchr/chr/g' > $1/genes/$GENE.abnormal 2> /dev/null
	
	if [[ ! -s $1/genes/$GENE.abnormal ]]; then
		rm $1/genes/$GENE.abnormal
	fi;
done < $3

# Print.
find $1/genes -name '*.abnormal' -exec wc -l '{}' \;

