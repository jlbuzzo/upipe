###############################################################################
#
# A Makefile to rule'em all!
#
###############################################################################





############################# PREAMBLE ######################################## 
# Ancillary files.
include variables.mk
include macros.mk
#include functions.mk
#include commands.mk

#$(call presentation, $(PIPELINE_NAME),asfd2)
#export INPUT_FILE


############################# MAIN CODE ######################################## 

all: processSample
	$(call show, All processes finished.)

# Post processing stage.
post_processing: main_processing
	$(call show, Post processing...)
	#$(post_proc)

# Main processing stage.
main_processing: pre_processing
	$(call show, Main processing...)
	#$(m_proc)

# Pre processing stage.
pre_processing: validations
	$(call show, Pre processing...)
	#make -f processSample

# Verify and validate inputs.
validations:
	$(call show, Validating inputs...)
	#$(validate)

# Extra processes for further customizations.
extra:
	$(call show, Extra processing...)
	#$(extra)


# Rule files.
include processSample.mk
#include mergeCall.mk



.PHONY: all
###############################################################################
