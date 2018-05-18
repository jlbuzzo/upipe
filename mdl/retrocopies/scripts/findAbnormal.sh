#!/bin/bash

############################## HEADER #########################################

# Arguments usage:
#
# $1: Source (file).
# $2: Index file.
# $3: Search criteria.
# $4: Output (dir).



############################## PREAMBLE #######################################

# Arguments validation.
if (( ${#@} < 4 )); then
  	echo -e "Error in $0: Not enough arguments!" >&2
	exit 1;
fi

# Counting time and iterations.
#count_time=$(time )
count=1

# Temporary folder.
FDR=$4/genes

# Clean or create the temporary directory to store files.
[[ -d "$FDR" ]] && rm -Rf $FDR
mkdir -p $FDR



############################## MAIN CODE ######################################

# 'chr' pattern extracted from BAM and formated files.
ptn_bam=$(samtools view -H $1 2> /dev/null | sed -rn '1,/@SQ/ s/SN:([^0-9]*)[0-9]*/\1/p' | cut -f 2)
ptn_fmt=$(sed -rn '1s/^(chr|Chr|CHR).*/\1/p' $2)


# For each gene in an annotation file, create a temporary file with abnormal alignments.
while read CHR START END GENE; do
	# Debug code.
	#echo "c-$CHR st-$START ed-$END gn-$GENE; ${ptn_bam:0:3}${CHR:3} $2 do"
	
	aux=$(samtools view $1 ${ptn_bam:0:3}${CHR:3}:$START-$END 2> /dev/null | awk '
		{
			if ( ( $7 == "=" && ( $9 > 100000 || $9 < -100000) ) || ( $7 != "=" && $7 != "*" ) )
			print "chr"$3"\t"$4"\tchr"$7"\t"$8
		}' | sed 's/chrchr/chr/g')
	
	# Write non empty entries to a file .abnormal.
	((${#aux} > 0)) && echo -e "$aux" > $FDR/${GENE}.abnormal
	
	# Debug code.
	#if (( ++count > 300 )); then break; fi
done < $2



###############################################################################

# Print, trimming white spaces and exit.
find $FDR -type f -name '*.abnormal' -print0 | xargs --null wc -l | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
