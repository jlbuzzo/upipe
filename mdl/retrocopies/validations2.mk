###############################################################################
#
# validations2
#
###############################################################################




# Preamble
define val
	find . -type f -name $1
endef


validations:
	$(info )
	$(info $(CALL) Validation process.)


.PHONY: validations
