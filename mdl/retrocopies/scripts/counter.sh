#!/bin/bash

###############################################################################
# Arguments usage:

# $1: is a folder list to search for retrocopies.


# Arguments validation.
if (( ${#@} < 1 )); then
  	echo -e "Error in $0: Not enough arguments!" >&2
	exit 1;
fi

# The default name of a RTC file.
file=putativeins.min.norep.exonic.dist.notsimilar.orientation


# Search criterium.
#srch=$(sed -rn -e '/IN IN IN IN_/s//>>/p' ${1}result/$file | awk -F ":" '{ print $2 }' | awk -F "_" '{if (($2 != 0) && ($3 != 0)) print }' | wc -l)



###############################################################################
# Preamble:

# Search for retrocopies in results directory of RTC folders.
while [ ! -z $1 ]; do
	# Test for directory existence.
	if [[ -d $1 && -s "$1/result/$file" ]]; then
		# Print directory name and retrocopy number.
		echo -e "Folder $1: $(sed -rn -e '/IN IN IN IN_/s//>>/p' ${1}result/$file | awk -F ":" '{ print $2 }' | awk -F "_" '{if (($2 != 0) && ($3 != 0)) print }' | wc -l)";
	fi
	
	shift
done <<< "$@"
