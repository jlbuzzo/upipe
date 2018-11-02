#!/usr/bin/env bash
############################## HEADER ##########################################
#
# A simple script to simulate a process over time.
#
# Bash script data:
# 		Main target: none
# 		File: sleep.sh
# 		Module: miscellaneous
#		Project: generic
# 		Author: Jose L. L. Buzzo
# 		Organization: GalanteLab
# 		Year: 2018
#
#
# Usage:
# 		Input args:
#			$1: Time for sleep (integer > 0)
#			[$2: An initial message to print (string)]
#			[$3: An final message to print (string)]
#
#		Output args:
#			none
################################################################################





############################## PREAMBLE ########################################

# Arguments validation.
if (( $# < 1 )); then
  	echo -e "Error in $0: Input arguments missing!\n" >&2
	echo -e "Usage:
	
	Input args:
		\$1\tSleeping time (integer > 0).
		[\$2]\tAn initial message to print (string).
		[\$3]\tA final message to print (string).
	
	Output args:
		none" >&2
	exit 1;
fi



############################## MAIN ############################################

# Print messages, if given.
if (( "$1" > 0 )); then
	[[ -n "$2" ]] && echo -e "$2"
	
	sleep $1
	
	[[ -n "$3" ]] && echo -e "$3"
	
	# Must finish with an exit status 0.
	exit 0
else
  	echo -e "Error in $0: Invalid argument value: \"$1\".\n" >&2
  	echo -e "First argument must be an integer > 0." >&2
	exit 1;
fi
