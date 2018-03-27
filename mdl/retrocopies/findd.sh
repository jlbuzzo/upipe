
OUT=$3/genes


#Clean and create temporary path to store temporary files
rm -Rf $OUT
mkdir -p $OUT

#For each gene, create a temporary file with abnormal alignments
for j in $(cat $1); do
	GENE=$(echo $j | awk -F "[*][*]" '{print $4}')
	CHR=$(echo $j | awk -F "[*][*]" '{print $1}')
	START=$(echo $j | awk -F "[*][*]" '{print $2}')
	END=$(echo $j | awk -F "[*][*]" '{print $3}')
	
	samtools view $2 ${CHR}:${START}-${END} |
	awk '{if ( ( $7 == "=" && ( $9 > 100000 || $9 < -100000) ) || ( $7 != "=" && $7 != "*" ) ) {print "chr"$3"\t"$4"\tchr"$7"\t"$8}}' | sed 's/chrchr/chr/g' > $OUT/${GENE}.abnormal 2>/dev/null
	
	#if [[ ! -s $OUT/${GENE}.abnormal ]]; then
	#	rm $OUT/${GENE}.abnormal
	#fi
done


find $OUT -type f -name '*.abnormal' -exec wc -l '{}' ';'
