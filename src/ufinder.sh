#!/bin/bash


# Universal file finder.

# Variable setting.
INPUT=$1

echo $INPUT

# Analyse INPUT variable ($1).
if [ -z "$INPUT" ]; then
	# Can't be a zero length string.
	echo "Empty INPUT variable."
elif [ -s "$INPUT" ]; then
	# Variable INPUT exists as a non-empty file, directory or URL.
	echo Exists.
else
	echo Empty file.
fi


exit
