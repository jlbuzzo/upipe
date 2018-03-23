#!/bin/bash

# Arguments:
# $1: OUTPUT_DIR.
# $2: FINAL_BAM_FILE.
# $3: *.formated file.

# Counting time.
timeini=$(date +%s.%n)

# Destiny folder.
FDR=$3/genes

# Clean and create temporary path to store temporary files.
[[ -d "$FDR" ]] && rm -Rf $FDR
mkdir -p $FDR


# For each gene, create a temporary file with abnormal alignments.
sed -e 's/\*\*/ /g' $1 | while read CHR START END GENE; do
	#CHR=$(sed 's/chr//i' <<< $CHR)
	#echo "c-$CHR st-$START ed-$END gn-$GENE; $2 do"
	
	samtools view $2 $CHR:$START-$END | awk '
		{
			if ( ( $7 == "=" && ( $9 > 100000 || $9 < -100000) ) || ( $7 != "=" && $7 != "*" ) )
				print "chr"$3"\t"$4"\tchr"$7"\t"$8
		}' | sed 's/chrchr/chr/g' > $FDR/$GENE.abnormal 2> /dev/null
	
	if [[ ! -s $FDR/$GENE.abnormal ]]; then
		rm $FDR/$GENE.abnormal
	fi;
done

# Count time.
timeend=$(date +%s.%n)
echo "timeini-timeend" > rtc.log

# Print and exit.
find $FDR -type f -name '*.abnormal' -exec wc -l '{}' ';'
